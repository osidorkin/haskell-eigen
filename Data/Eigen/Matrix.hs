{-# LANGUAGE RecordWildCards, MultiParamTypeClasses, Rank2Types #-}
module Data.Eigen.Matrix (
    -- * Matrix type
    Matrix(..),
    -- * Matrix conversions
    fromList,
    toList,
    generate,
    -- * Standard matrices and special cases
    empty,
    zero,
    ones,
    identity,
    constant,
    -- * Accessing matrix data
    cols,
    rows,
    coeff,
    minCoeff,
    maxCoeff,
    col,
    row,
    block,
    topRows,
    bottomRows,
    leftCols,
    rightCols,
    -- * Matrix properties
    norm,
    squaredNorm,
    determinant,
    -- * Matrix operations
    add,
    sub,
    mul,
    -- * Matrix transformations
    inverse,
    adjoint,
    conjugate,
    transpose,
    normalize,
    modify,
    -- * Mutable matrices
    thaw,
    freeze,
    unsafeThaw,
    unsafeFreeze
) where

import Data.List (intercalate)
import Data.Tuple
import Foreign.Ptr
import Foreign.C.Types
import Foreign.C.String
import Control.Monad
import Control.Monad.ST
import Control.Monad.Primitive
import Control.Applicative hiding (empty)
import qualified Data.Vector.Storable as VS
import qualified Data.Vector.Storable.Mutable as VSM
import Data.Eigen.Internal
import Data.Eigen.Matrix.Mutable

-- | constant Matrix class to be used in pure computations, uses the same column major memory layout as Eigen MatrixXd
data Matrix = Matrix {
    m_rows :: Int,
    m_cols :: Int,
    m_vals :: VS.Vector CDouble --
};

-- | pretty prints the matrix
instance Show Matrix where
    show m@Matrix{..} = concat [
        "Matrix ", show m_rows, "x", show m_cols, "\n", intercalate "\n" $ map (intercalate "\t" . map show) $ toList m, "\n"]

-- | only the following functions are defined for Num instance: (*), (+), (-)
instance Num Matrix where
    (*) = mul
    (+) = add
    (-) = sub
    fromInteger = undefined
    signum = undefined
    abs = undefined

-- | empty 0x0 matrix
empty :: Matrix
empty = Matrix 0 0 VS.empty

-- | matrix where all coeffs are filled with given value
constant :: Int -> Int -> Double -> Matrix
constant rows cols val = Matrix rows cols $ VS.replicate (rows * cols) (cast val)

-- | matrix where all coeff are 0
zero :: Int -> Int -> Matrix
zero rows cols = constant rows cols 0

-- | matrix where all coeff are 1
ones :: Int -> Int -> Matrix
ones rows cols = constant rows cols 1

-- | square matrix with 1 on main diagonal and 0 elsewhere
identity :: Int -> Matrix
identity size = Matrix size size $ VS.create $ do
    vm <- VSM.replicate (size * size) 0
    forM_ [0..pred size] $ \n ->
        VSM.write vm (n * size + n) 1
    return vm

-- | number of rows for the matrix
rows :: Matrix -> Int
rows = m_rows

-- | number of columns for the matrix
cols :: Matrix -> Int
cols = m_cols

-- | matrix coefficient at specific row and col
coeff :: Int -> Int -> Matrix -> Double
coeff row col Matrix{..} = cast $ m_vals VS.! (col * m_rows + row)

-- | list of coefficients for the given col
col :: Int -> Matrix -> [Double]
col c m@Matrix{..} = [coeff r c m | r <- [0..pred m_rows]]

-- | list of coefficients for the given row
row :: Int -> Matrix -> [Double]
row r m@Matrix{..} = [coeff r c m | c <- [0..pred m_cols]]

-- | extract rectangular block from matrix defined by startRow startCol blockRows blockCols
block :: Int -> Int -> Int -> Int -> Matrix -> Matrix
block startRow startCol blockRows blockCols m =
    generate blockRows blockCols $ \row col ->
        coeff (startRow + row) (startCol + col) m

-- | the maximum of all coefficients of matrix
maxCoeff :: Matrix -> Double
maxCoeff Matrix{..} = cast $ VS.maximum m_vals

-- | the minimum of all coefficients of matrix
minCoeff :: Matrix -> Double
minCoeff Matrix{..} = cast $ VS.minimum m_vals

-- | top n rows of matrix
topRows :: Int -> Matrix -> Matrix
topRows rows m@Matrix{..} = block 0 0 rows m_cols m

-- | bottom n rows of matrix
bottomRows :: Int -> Matrix -> Matrix
bottomRows rows m@Matrix{..} = block (m_rows - rows) 0 rows m_cols m

-- | left n columns of matrix
leftCols :: Int -> Matrix -> Matrix
leftCols cols m@Matrix{..} = block 0 0 m_rows cols m

-- | right n columns of matrix
rightCols :: Int -> Matrix -> Matrix
rightCols cols m@Matrix{..} = block 0 (m_cols - cols) m_rows cols m

-- | construct matrix from a list of rows, column count is detected as maximum row length
fromList :: [[Double]] -> Matrix
fromList list = Matrix rows cols vals where
    rows = length list
    cols = maximum $ map length list
    vals = VS.create $ do
        vm <- VSM.replicate (rows * cols) 0
        forM_ (zip [0..] list) $ \(row, vals) ->
            forM_ (zip [0..] vals) $ \(col, val) ->
                VSM.write vm (col * rows + row) (cast val)
        return vm

-- | converts matrix to a list of its rows
toList :: Matrix -> [[Double]]
toList Matrix{..} = [[cast $ m_vals VS.! (col * m_rows + row) | col <- [0..pred m_cols]] | row <- [0..pred m_rows]]

-- | craete matrix using generator function f :: row -> col -> val
generate :: Int -> Int -> (Int -> Int -> Double) -> Matrix
generate rows cols f = Matrix rows cols $ VS.create $ do
    vals <- VSM.new (rows * cols)
    forM_ [0..pred rows] $ \row ->
        forM_ [0..pred cols] $ \col ->
            VSM.write vals (col * rows + row) (cast $ f row col)
    return vals


-- | for vectors, the l2 norm, and for matrices the Frobenius norm. In both cases, it consists in the square root of the sum of the square of all the matrix entries. For vectors, this is also equals to the square root of the dot product of this with itself.
norm :: Matrix -> Double
norm = _unop c_norm

-- | for vectors, the squared l2 norm, and for matrices the Frobenius norm. In both cases, it consists in the sum of the square of all the matrix entries. For vectors, this is also equals to the dot product of this with itself.
squaredNorm :: Matrix -> Double
squaredNorm = _unop c_squaredNorm

-- | the determinant of the matrix
determinant :: Matrix -> Double
determinant m@Matrix{..}
    | m_cols == m_rows = _unop c_determinant m
    | otherwise = error "you tried calling determinant on non-square matrix"

-- | return a - b
add :: Matrix -> Matrix -> Matrix
add = _binop c_add

-- | return a + b
sub :: Matrix -> Matrix -> Matrix
sub = _binop c_sub

-- | return a * b
mul :: Matrix -> Matrix -> Matrix
mul = _binop c_mul


{- | inverse of the matrix

For small fixed sizes up to 4x4, this method uses cofactors. In the general case, this method uses class 'PartialPivLU'
-}
inverse :: Matrix -> Matrix
inverse m@Matrix{..}
    | m_rows == m_cols = _modify id c_inverse m
    | otherwise = error "you tried calling inverse on non-square matrix"

-- | adjoint of the matrix
adjoint :: Matrix -> Matrix
adjoint = _modify swap c_adjoint

-- | transpose of the matrix
transpose :: Matrix -> Matrix
transpose = _modify swap c_transpose

-- | conjugate of the matrix
conjugate :: Matrix -> Matrix
conjugate = _modify id c_conjugate

-- | nomalize the matrix by deviding it on its 'norm'
normalize :: Matrix -> Matrix
normalize Matrix{..} = performIO $ do
    vals <- VS.thaw m_vals
    VSM.unsafeWith vals $ \p ->
        call $ c_normalize p (cast m_rows) (cast m_cols)
    Matrix m_rows m_cols <$> VS.unsafeFreeze vals

-- | Apply a destructive operation to a matrix. The operation will be performed in place if it is safe to do so and will modify a copy of the matrix otherwise.
modify :: (forall s. MMatrix s -> ST s ()) -> Matrix -> Matrix
modify f m@Matrix{..} = m { m_vals = VS.modify f' m_vals } where
    f' vals = f (MMatrix m_rows m_cols vals)

-- | Yield an immutable copy of the mutable matrix
freeze :: PrimMonad m => MMatrix (PrimState m) -> m Matrix
freeze MMatrix{..} = VS.freeze mm_vals >>= \vals -> return $ Matrix mm_rows mm_cols vals

-- | Yield a mutable copy of the immutable matrix
thaw :: PrimMonad m => Matrix -> m (MMatrix (PrimState m))
thaw Matrix{..} = VS.thaw m_vals >>= \vals -> return $ MMatrix m_rows m_cols vals

-- | Unsafe convert a mutable matrix to an immutable one without copying. The mutable matrix may not be used after this operation.
unsafeFreeze :: PrimMonad m => MMatrix (PrimState m) -> m Matrix
unsafeFreeze MMatrix{..} = VS.unsafeFreeze mm_vals >>= \vals -> return $ Matrix mm_rows mm_cols vals

-- | Unsafely convert an immutable matrix to a mutable one without copying. The immutable matrix may not be used after this operation.
unsafeThaw :: PrimMonad m => Matrix -> m (MMatrix (PrimState m))
unsafeThaw Matrix{..} = VS.unsafeThaw m_vals >>= \vals -> return $ MMatrix m_rows m_cols vals

_unop :: (Ptr CDouble -> CInt -> CInt -> IO CDouble) -> Matrix -> Double
_unop f Matrix{..} = performIO $ VS.unsafeWith m_vals $ \p ->
    cast <$> f p (cast m_rows) (cast m_cols)

_binop :: (Ptr CDouble -> CInt -> CInt -> Ptr CDouble -> CInt -> CInt -> IO CString) -> Matrix -> Matrix -> Matrix
_binop f m1 m2 = performIO $ do
    vals <- VS.thaw (m_vals m1)
    VSM.unsafeWith vals $ \lhs ->
        VS.unsafeWith (m_vals m2) $ \rhs ->
            call $ f
                lhs (cast $ m_rows m1) (cast $ m_cols m1)
                rhs (cast $ m_rows m2) (cast $ m_cols m2)
    Matrix (m_rows m1) (m_cols m1) <$> VS.unsafeFreeze vals

_modify :: ((Int,Int) -> (Int,Int)) -> (Ptr CDouble -> CInt -> CInt -> Ptr CDouble -> CInt -> CInt -> IO CString) -> Matrix -> Matrix
_modify f g Matrix{..} = performIO $ do
    let (rows, cols) = f (m_rows, m_cols)
    vals <- VSM.new (rows * cols)
    VS.unsafeWith m_vals $ \src ->
        VSM.unsafeWith vals $ \dst ->
            call $ g
                dst (cast rows) (cast cols)
                src (cast m_rows) (cast m_cols)
    Matrix rows cols <$> VS.unsafeFreeze vals
