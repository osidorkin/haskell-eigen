{-# LANGUAGE RecordWildCards, MultiParamTypeClasses, Rank2Types #-}
module Data.Eigen.Matrix (
    -- * Matrix type
    Matrix(..),
    valid,
    -- * Matrix conversions
    fromList,
    toList,
    generate,
    -- * Standard matrices and special cases
    empty,
    null,
    square,
    zero,
    ones,
    identity,
    constant,
    random,
    -- * Accessing matrix data
    cols,
    rows,
    (!),
    coeff,
    unsafeCoeff,
    col,
    row,
    block,
    topRows,
    bottomRows,
    leftCols,
    rightCols,
    -- * Matrix properties
    sum,
    prod,
    mean,
    minCoeff,
    maxCoeff,
    trace,
    norm,
    squaredNorm,
    blueNorm,
    hypotNorm,
    determinant,
    -- * Boolean reductions
    all,
    any,
    count,
    -- * Matrix operations
    add,
    sub,
    mul,
    -- * Matrix transformations
    diagonal,
    transpose,
    inverse,
    adjoint,
    conjugate,
    normalize,
    modify,
    -- * Mutable matrices
    thaw,
    freeze,
    unsafeThaw,
    unsafeFreeze,
    unsafeWith,
) where

import Prelude hiding (null, sum, all, any)
import Data.List (intercalate)
import Data.Tuple
import Foreign.Ptr
import Foreign.C.Types
import Foreign.C.String
import Text.Printf
import Control.Monad
import Control.Monad.ST
import Control.Monad.Primitive
import Control.Applicative hiding (empty)
import qualified Data.Vector.Storable as VS
import qualified Data.Vector.Storable.Mutable as VSM
import qualified Data.Eigen.Internal as I
import qualified Data.Eigen.Matrix.Mutable as M

-- | Matrix to be used in pure computations, uses column major memory layout
data Matrix = Matrix {
    m_rows :: Int,
    m_cols :: Int,
    m_vals :: VS.Vector CDouble
};

-- | Pretty prints the matrix
instance Show Matrix where
    show m@Matrix{..} = concat [
        "Matrix ", show m_rows, "x", show m_cols, 
        "\n", intercalate "\n" $ map (intercalate "\t" . map show) $ toList m, "\n"]

-- | Nm instance for the matrix
instance Num Matrix where
    (*) = mul
    (+) = add
    (-) = sub
    fromInteger = constant 1 1 . fromInteger
    signum m@Matrix{..} = m { m_vals = VS.map signum m_vals }
    abs m@Matrix{..} = m { m_vals = VS.map abs m_vals }

-- | Empty 0x0 matrix
empty :: Matrix
empty = Matrix 0 0 VS.empty

-- | Is matrix empty?
null :: Matrix -> Bool
null Matrix{..} = m_rows == 0 && m_cols == 0

-- | Is matrix square?
square :: Matrix -> Bool
square Matrix{..} = m_rows == m_cols

-- | Matrix where all coeffs are filled with given value
constant :: Int -> Int -> Double -> Matrix
constant rows cols val = Matrix rows cols $ VS.replicate (rows * cols) (I.cast val)

-- | Matrix where all coeff are 0
zero :: Int -> Int -> Matrix
zero rows cols = constant rows cols 0

-- | Matrix where all coeff are 1
ones :: Int -> Int -> Matrix
ones rows cols = constant rows cols 1

-- | Square matrix with 1 on main diagonal and 0 elsewhere
identity :: Int -> Matrix
identity size = Matrix size size $ VS.create $ do
    vm <- VSM.replicate (size * size) 0
    forM_ [0..pred size] $ \n ->
        VSM.write vm (n * size + n) 1
    return vm

-- | The random matrix of a given size
random :: Int -> Int -> IO Matrix
random rows cols = do
    m <- M.new rows cols
    I.call $ M.unsafeWith m I.c_random
    unsafeFreeze m

-- | Number of rows for the matrix
rows :: Matrix -> Int
rows = m_rows

-- | Number of columns for the matrix
cols :: Matrix -> Int
cols = m_cols

-- | Mtrix size as (rows, cols) pair
dims :: Matrix -> (Int, Int)
dims Matrix{..} = (m_rows, m_cols)

-- | Matrix coefficient at specific row and col
(!) :: Matrix -> (Int,Int) -> Double
(!) m (row,col) = coeff row col m

-- | Matrix coefficient at specific row and col
coeff :: Int -> Int -> Matrix -> Double
coeff row col m@Matrix{..}
    | not (valid m) = error "matrix is not valid"
    | row < 0 || row >= m_rows = error $ printf "Matrix.coeff: row %d is out of bounds [0..%d)" row m_rows
    | col < 0 || col >= m_cols = error $ printf "Matrix.coeff: col %d is out of bounds [0..%d)" col m_cols
    | otherwise = unsafeCoeff row col m

-- | Unsafe version of coeff function. No bounds check performed so SEGFAULT possible
unsafeCoeff :: Int -> Int -> Matrix -> Double
unsafeCoeff row col Matrix{..} = I.cast $ VS.unsafeIndex m_vals $ col * m_rows + row

-- | List of coefficients for the given col
col :: Int -> Matrix -> [Double]
col c m@Matrix{..} = [coeff r c m | r <- [0..pred m_rows]]

-- | List of coefficients for the given row
row :: Int -> Matrix -> [Double]
row r m@Matrix{..} = [coeff r c m | c <- [0..pred m_cols]]

-- | Extract rectangular block from matrix defined by startRow startCol blockRows blockCols
block :: Int -> Int -> Int -> Int -> Matrix -> Matrix
block startRow startCol blockRows blockCols m =
    generate blockRows blockCols $ \row col ->
        coeff (startRow + row) (startCol + col) m

-- | Verify matrix dimensions and memory layout
valid :: Matrix -> Bool
valid Matrix{..} = m_rows >= 0 && m_cols >= 0 && VS.length m_vals == m_rows * m_cols

-- | The maximum of all coefficients of matrix
maxCoeff :: Matrix -> Double
maxCoeff Matrix{..} = I.cast $ VS.maximum m_vals

-- | The minimum of all coefficients of matrix
minCoeff :: Matrix -> Double
minCoeff Matrix{..} = I.cast $ VS.minimum m_vals

-- | Top n rows of matrix
topRows :: Int -> Matrix -> Matrix
topRows rows m@Matrix{..} = block 0 0 rows m_cols m

-- | Bottom n rows of matrix
bottomRows :: Int -> Matrix -> Matrix
bottomRows rows m@Matrix{..} = block (m_rows - rows) 0 rows m_cols m

-- | Left n columns of matrix
leftCols :: Int -> Matrix -> Matrix
leftCols cols m@Matrix{..} = block 0 0 m_rows cols m

-- | Right n columns of matrix
rightCols :: Int -> Matrix -> Matrix
rightCols cols m@Matrix{..} = block 0 (m_cols - cols) m_rows cols m

-- | Construct matrix from a list of rows, column count is detected as maximum row length. Missing values are filled with 0
fromList :: [[Double]] -> Matrix
fromList list = Matrix rows cols vals where
    rows = length list
    cols = maximum $ map length list
    vals = VS.create $ do
        vm <- VSM.replicate (rows * cols) 0
        forM_ (zip [0..] list) $ \(row, vals) ->
            forM_ (zip [0..] vals) $ \(col, val) ->
                VSM.write vm (col * rows + row) (I.cast val)
        return vm

-- | Convert matrix to a list of rows
toList :: Matrix -> [[Double]]
toList Matrix{..} = [[I.cast $ m_vals VS.! (col * m_rows + row) | col <- [0..pred m_cols]] | row <- [0..pred m_rows]]

-- | Create matrix using generator function f :: row -> col -> val
generate :: Int -> Int -> (Int -> Int -> Double) -> Matrix
generate rows cols f = Matrix rows cols $ VS.create $ do
    vals <- VSM.new (rows * cols)
    forM_ [0..pred rows] $ \row ->
        forM_ [0..pred cols] $ \col ->
            VSM.write vals (col * rows + row) (I.cast $ f row col)
    return vals

-- | The sum of all coefficients of the matrix
sum :: Matrix -> Double
sum = _prop I.c_sum

-- | The product of all coefficients of the matrix
prod :: Matrix -> Double
prod = _prop I.c_prod

-- | The mean of all coefficients of the matrix
mean :: Matrix -> Double
mean = _prop I.c_prod

-- | The trace of a matrix is the sum of the diagonal coefficients and can also be computed as sum (diagonal m)
trace :: Matrix -> Double
trace = _prop I.c_trace

-- | Applied to a predicate and a matrix, all determines if all elements of the matrix satisfies the predicate
all :: (Double -> Bool) -> Matrix -> Bool
all f = VS.all (f . I.cast) . m_vals

-- | Applied to a predicate and a matrix, any determines if any element of the matrix satisfies the predicate
any :: (Double -> Bool) -> Matrix -> Bool
any f = VS.any (f . I.cast) . m_vals

-- | Returns the number of coefficients in a given matrix that evaluate to true
count :: (Double -> Bool) -> Matrix -> Int
count f = VS.foldl' (\n x -> if f (I.cast x) then succ n else n) 0 . m_vals

-- | For vectors, the l2 norm, and for matrices the Frobenius norm. In both cases, it consists in the square root of the sum of the square of all the matrix entries. For vectors, this is also equals to the square root of the dot product of this with itself.
norm :: Matrix -> Double
norm = _prop I.c_norm

-- | For vectors, the squared l2 norm, and for matrices the Frobenius norm. In both cases, it consists in the sum of the square of all the matrix entries. For vectors, this is also equals to the dot product of this with itself.
squaredNorm :: Matrix -> Double
squaredNorm = _prop I.c_squaredNorm

-- | The l2 norm of the matrix using the Blue's algorithm. A Portable Fortran Program to Find the Euclidean Norm of a Vector, ACM TOMS, Vol 4, Issue 1, 1978.
blueNorm :: Matrix -> Double
blueNorm = _prop I.c_blueNorm

-- | The l2 norm of the matrix avoiding undeflow and overflow. This version use a concatenation of hypot calls, and it is very slow.
hypotNorm :: Matrix -> Double
hypotNorm = _prop I.c_hypotNorm

-- | The determinant of the matrix
determinant :: Matrix -> Double
determinant m
    | square m = _prop I.c_determinant m
    | otherwise = error "Matrix.determinant: non-square matrix"

-- | Adding two matrices by adding the corresponding entries together. You can use @(+)@ function as well.
add :: Matrix -> Matrix -> Matrix
add m1 m2
    | dims m1 == dims m2 = _binop const I.c_add m1 m2
    | otherwise = error "Matrix.add: matrix should have the same size"

-- | Subtracting two matrices by subtracting the corresponding entries together. You can use @(-)@ function as well.
sub :: Matrix -> Matrix -> Matrix
sub m1 m2
    | dims m1 == dims m2 = _binop const I.c_sub m1 m2
    | otherwise = error "Matrix.add: matrix should have the same size"

-- | Matrix multiplication. You can use @(*)@ function as well.
mul :: Matrix -> Matrix -> Matrix
mul m1 m2
    | m_cols m1 == m_rows m2 = _binop (\(rows, _) (_, cols) -> (rows, cols)) I.c_mul m1 m2
    | otherwise = error "Matrix.mul: number of columns for lhs matrix should be the same as number of rows for rhs matrix"

-- | Diagonal of the matrix
diagonal :: Matrix -> Matrix
diagonal = _unop (\(rows, cols) -> (min rows cols, 1)) I.c_diagonal

{- | Inverse of the matrix

For small fixed sizes up to 4x4, this method uses cofactors. In the general case, this method uses PartialPivLU decomposition
-}
inverse :: Matrix -> Matrix
inverse m@Matrix{..}
    | m_rows == m_cols = _unop id I.c_inverse m
    | otherwise = error "Matrix.inverse: non-square matrix"

-- | Adjoint of the matrix
adjoint :: Matrix -> Matrix
adjoint = _unop swap I.c_adjoint

-- | Transpose of the matrix
transpose :: Matrix -> Matrix
transpose = _unop swap I.c_transpose

-- | Conjugate of the matrix
conjugate :: Matrix -> Matrix
conjugate = _unop id I.c_conjugate

-- | Nomalize the matrix by deviding it on its 'norm'
normalize :: Matrix -> Matrix
normalize Matrix{..} = I.performIO $ do
    vals <- VS.thaw m_vals
    VSM.unsafeWith vals $ \p ->
        I.call $ I.c_normalize p (I.cast m_rows) (I.cast m_cols)
    Matrix m_rows m_cols <$> VS.unsafeFreeze vals

-- | Apply a destructive operation to a matrix. The operation will be performed in place if it is safe to do so and will modify a copy of the matrix otherwise.
modify :: (forall s. M.MMatrix s -> ST s ()) -> Matrix -> Matrix
modify f m@Matrix{..} = m { m_vals = VS.modify (f . M.MMatrix m_rows m_cols) m_vals } where

-- | Yield an immutable copy of the mutable matrix
freeze :: PrimMonad m => M.MMatrix (PrimState m) -> m Matrix
freeze M.MMatrix{..} = VS.freeze mm_vals >>= return . Matrix mm_rows mm_cols

-- | Yield a mutable copy of the immutable matrix
thaw :: PrimMonad m => Matrix -> m (M.MMatrix (PrimState m))
thaw Matrix{..} = VS.thaw m_vals >>= return . M.MMatrix m_rows m_cols

-- | Unsafe convert a mutable matrix to an immutable one without copying. The mutable matrix may not be used after this operation.
unsafeFreeze :: PrimMonad m => M.MMatrix (PrimState m) -> m Matrix
unsafeFreeze M.MMatrix{..} = VS.unsafeFreeze mm_vals >>= return . Matrix mm_rows mm_cols

-- | Unsafely convert an immutable matrix to a mutable one without copying. The immutable matrix may not be used after this operation.
unsafeThaw :: PrimMonad m => Matrix -> m (M.MMatrix (PrimState m))
unsafeThaw Matrix{..} = VS.unsafeThaw m_vals >>= return . M.MMatrix m_rows m_cols

-- | Pass a pointer to the matrix's data to the IO action. The data may not be modified through the pointer.
unsafeWith :: Matrix -> (Ptr CDouble -> CInt -> CInt -> IO a) -> IO a
unsafeWith m@Matrix{..} f
    | not (valid m) = fail "matrix layout is invalid"
    | otherwise = VS.unsafeWith m_vals $ \p -> f p (I.cast m_rows) (I.cast m_cols)

_prop :: (Ptr CDouble -> CInt -> CInt -> IO CDouble) -> Matrix -> Double
_prop f m = I.cast $ I.performIO $ unsafeWith m f

_binop :: ((Int, Int) -> (Int, Int) -> (Int, Int)) -> (Ptr CDouble -> CInt -> CInt -> Ptr CDouble -> CInt -> CInt -> Ptr CDouble -> CInt -> CInt -> IO CString) -> Matrix -> Matrix -> Matrix
_binop f g m1 m2 = I.performIO $ do
    m0 <- uncurry M.new $ f (dims m1) (dims m2)
    M.unsafeWith m0 $ \vals0 rows0 cols0 ->
        unsafeWith m1 $ \vals1 rows1 cols1 ->
            unsafeWith m2 $ \vals2 rows2 cols2 ->
                I.call $ g
                    vals0 rows0 cols0
                    vals1 rows1 cols1
                    vals2 rows2 cols2
    unsafeFreeze m0

_unop :: ((Int,Int) -> (Int,Int)) -> (Ptr CDouble -> CInt -> CInt -> Ptr CDouble -> CInt -> CInt -> IO CString) -> Matrix -> Matrix
_unop f g m1 = I.performIO $ do
    m0 <- uncurry M.new $ f (dims m1)
    M.unsafeWith m0 $ \vals0 rows0 cols0 ->
        unsafeWith m1 $ \vals1 rows1 cols1 ->
            I.call $ g
                vals0 rows0 cols0
                vals1 rows1 cols1
    unsafeFreeze m0
