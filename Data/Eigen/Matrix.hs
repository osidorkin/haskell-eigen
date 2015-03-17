{-# LANGUAGE RecordWildCards, MultiParamTypeClasses, Rank2Types, ScopedTypeVariables, FunctionalDependencies, GADTs #-}
module Data.Eigen.Matrix (
    -- * Matrix type
    -- | Matrix aliases follows Eigen naming convention
    Matrix(..),
    MatrixXf,
    MatrixXd,
    MatrixXcf,
    MatrixXcd,
    I.Elem,
    I.CComplex,
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
    dims,
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
    -- * Generic reductions
    fold,
    fold',
    ifold,
    ifold',
    fold1,
    fold1',
    -- * Boolean reductions
    all,
    any,
    count,
    -- * Basic matrix algebra
    add,
    sub,
    mul,
    -- * Mapping over elements
    map,
    imap,
    filter,
    ifilter,
    -- * Matrix transformations
    diagonal,
    transpose,
    inverse,
    adjoint,
    conjugate,
    normalize,
    modify,
    convert,
    upperTriangle,
    lowerTriangle,
    -- * Matrix serialization
    encode,
    decode,
    -- * Mutable matrices
    thaw,
    freeze,
    unsafeThaw,
    unsafeFreeze,
    -- * Raw pointers
    unsafeWith,
) where

import qualified Prelude as P
import qualified Data.List as L
import Prelude hiding (null, sum, all, any, map, filter)
import Data.Tuple
import Data.STRef
import Data.Complex hiding (conjugate)
import Foreign.Ptr
import Foreign.ForeignPtr
import Foreign.C.Types
import Foreign.C.String
import Foreign.Storable
import Foreign.Marshal.Alloc
import Text.Printf
import Control.Monad
import Control.Monad.ST
import Control.Monad.Primitive
import Control.Applicative hiding (empty)
import qualified Data.Vector.Storable as VS
import qualified Data.Vector.Storable.Mutable as VSM
import qualified Data.Eigen.Internal as I
import qualified Data.Eigen.Matrix.Mutable as M
import qualified Data.ByteString as BS
import qualified Data.ByteString.Lazy as BSL
import qualified Data.ByteString.Internal as BSI

-- | Matrix to be used in pure computations, uses column major memory layout, features copy-free FFI with C++ <http://eigen.tuxfamily.org Eigen> library.

data Matrix a b where
    Matrix :: I.Elem a b => !Int -> !Int -> !(VS.Vector b) -> Matrix a b

-- | Alias for single precision matrix
type MatrixXf = Matrix Float CFloat
-- | Alias for double precision matrix
type MatrixXd = Matrix Double CDouble
-- | Alias for single previsiom matrix of complex numbers
type MatrixXcf = Matrix (Complex Float) (I.CComplex CFloat)
-- | Alias for double prevision matrix of complex numbers
type MatrixXcd = Matrix (Complex Double) (I.CComplex CDouble)

-- | Pretty prints the matrix
instance (I.Elem a b, Show a) => Show (Matrix a b) where
    show m@(Matrix rows cols _) = concat [
        "Matrix ", show rows, "x", show cols,
        "\n", L.intercalate "\n" $ P.map (L.intercalate "\t" . P.map show) $ toList m, "\n"]


-- | Shortcuts for basic matrix math
instance I.Elem a b => Num (Matrix a b) where
    (*) = mul
    (+) = add
    (-) = sub
    fromInteger = constant 1 1 . fromInteger
    signum = map signum
    abs = map abs
    negate = map negate

-- | Empty 0x0 matrix
empty :: I.Elem a b => Matrix a b
empty = Matrix 0 0 VS.empty

-- | Is matrix empty?
null :: I.Elem a b => Matrix a b -> Bool
null (Matrix rows cols _) = rows == 0 && cols == 0

-- | Is matrix square?
square :: I.Elem a b => Matrix a b -> Bool
square (Matrix rows cols _) = rows == cols

-- | Matrix where all coeffs are filled with given value
constant :: I.Elem a b => Int -> Int -> a -> Matrix a b
constant rows cols val = Matrix rows cols $ VS.replicate (rows * cols) (I.cast val)

-- | Matrix where all coeff are 0
zero :: I.Elem a b => Int -> Int -> Matrix a b
zero rows cols = constant rows cols 0

-- | Matrix where all coeff are 1
ones :: I.Elem a b => Int -> Int -> Matrix a b
ones rows cols = constant rows cols 1

-- | The identity matrix (not necessarily square).
identity :: I.Elem a b => Int -> Int -> Matrix a b
identity rows cols = I.performIO $ do
    m <- M.new rows cols
    I.call $ M.unsafeWith m I.identity
    unsafeFreeze m

-- | The random matrix of a given size
random :: I.Elem a b => Int -> Int -> IO (Matrix a b)
random rows cols = do
    m <- M.new rows cols
    I.call $ M.unsafeWith m I.random
    unsafeFreeze m

-- | Number of rows for the matrix
rows :: I.Elem a b => Matrix a b -> Int
rows (Matrix rows _ _) = rows

-- | Number of columns for the matrix
cols :: I.Elem a b => Matrix a b -> Int
cols (Matrix _ cols _) = cols

vals :: I.Elem a b => Matrix a b -> VS.Vector b
vals (Matrix _ _ vals) = vals

-- | Mtrix size as (rows, cols) pair
dims :: I.Elem a b => Matrix a b -> (Int, Int)
dims (Matrix rows cols _) = (rows, cols)

-- | Matrix coefficient at specific row and col
(!) :: I.Elem a b => Matrix a b -> (Int, Int) -> a
(!) m (row,col) = coeff row col m

-- | Matrix coefficient at specific row and col
coeff :: I.Elem a b => Int -> Int -> Matrix a b -> a
coeff row col m@(Matrix rows cols _)
    | not (valid m) = error "matrix is not valid"
    | row < 0 || row >= rows = error $ printf "Matrix.coeff: row %d is out of bounds [0..%d)" row rows
    | col < 0 || col >= cols = error $ printf "Matrix.coeff: col %d is out of bounds [0..%d)" col cols
    | otherwise = unsafeCoeff row col m

-- | Unsafe version of coeff function. No bounds check performed so SEGFAULT possible
unsafeCoeff :: I.Elem a b => Int -> Int -> Matrix a b -> a
unsafeCoeff row col (Matrix rows _ vals) = I.cast $ VS.unsafeIndex vals $ col * rows + row

-- | List of coefficients for the given col
col :: I.Elem a b => Int -> Matrix a b -> [a]
col c m@(Matrix rows _ _) = [coeff r c m | r <- [0..pred rows]]

-- | List of coefficients for the given row
row :: I.Elem a b => Int -> Matrix a b -> [a]
row r m@(Matrix _ cols _) = [coeff r c m | c <- [0..pred cols]]

-- | Extract rectangular block from matrix defined by startRow startCol blockRows blockCols
block :: I.Elem a b => Int -> Int -> Int -> Int -> Matrix a b -> Matrix a b
block startRow startCol blockRows blockCols m =
    generate blockRows blockCols $ \row col ->
        coeff (startRow + row) (startCol + col) m

-- | Verify matrix dimensions and values memory layout
valid :: I.Elem a b => Matrix a b -> Bool
valid (Matrix rows cols vals) = rows >= 0 && cols >= 0 && VS.length vals == rows * cols

-- | The maximum coefficient of the matrix
maxCoeff :: (I.Elem a b, Ord a) => Matrix a b -> a
maxCoeff = fold1' max

-- | The minimum coefficient of the matrix
minCoeff :: (I.Elem a b, Ord a) => Matrix a b -> a
minCoeff = fold1' min

-- | Top @N@ rows of matrix
topRows :: I.Elem a b => Int -> Matrix a b -> Matrix a b
topRows n m@(Matrix _ cols _) = block 0 0 n cols m

-- | Bottom @N@ rows of matrix
bottomRows :: I.Elem a b => Int -> Matrix a b -> Matrix a b
bottomRows n m@(Matrix rows cols _) = block (rows - n) 0 n cols m

-- | Left @N@ columns of matrix
leftCols :: I.Elem a b => Int -> Matrix a b -> Matrix a b
leftCols n m@(Matrix rows _ _) = block 0 0 rows n m

-- | Right @N@ columns of matrix
rightCols :: I.Elem a b => Int -> Matrix a b -> Matrix a b
rightCols n m@(Matrix rows cols _) = block 0 (cols - n) rows n m

-- | Construct matrix from a list of rows, column count is detected as maximum row length. Missing values are filled with 0
fromList :: I.Elem a b => [[a]] -> Matrix a b
fromList list = Matrix rows cols vals where
    rows = length list
    cols = L.foldl' max 0 $ P.map length list
    vals = VS.create $ do
        vm <- VSM.replicate (rows * cols) (I.cast (0 `asTypeOf` (head (head list))))
        forM_ (zip [0..] list) $ \(row, vals) ->
            forM_ (zip [0..] vals) $ \(col, val) ->
                VSM.write vm (col * rows + row) (I.cast val)
        return vm

-- | Convert matrix to a list of rows
toList :: I.Elem a b => Matrix a b -> [[a]]
toList (Matrix rows cols vals) = [[I.cast $ vals VS.! (col * rows + row) | col <- [0..pred cols]] | row <- [0..pred rows]]

-- | [generate rows cols (λ row col -> val)]
--
-- Create matrix using generator function @f row col val@
--
generate :: I.Elem a b => Int -> Int -> (Int -> Int -> a) -> Matrix a b
generate rows cols f = Matrix rows cols $ VS.create $ do
    vals <- VSM.new (rows * cols)
    forM_ [0..pred rows] $ \row ->
        forM_ [0..pred cols] $ \col ->
            VSM.write vals (col * rows + row) (I.cast $ f row col)
    return vals

-- | The sum of all coefficients of the matrix
sum :: I.Elem a b => Matrix a b -> a
sum = _prop I.sum

-- | The product of all coefficients of the matrix
prod :: I.Elem a b => Matrix a b -> a
prod = _prop I.prod

-- | The mean of all coefficients of the matrix
mean :: I.Elem a b => Matrix a b -> a
mean = _prop I.prod

-- | The trace of a matrix is the sum of the diagonal coefficients and can also be computed as sum (diagonal m)
trace :: I.Elem a b => Matrix a b -> a
trace = _prop I.trace

-- | Applied to a predicate and a matrix, all determines if all elements of the matrix satisfies the predicate
all :: I.Elem a b => (a -> Bool) -> Matrix a b -> Bool
all f = VS.all (f . I.cast) . vals

-- | Applied to a predicate and a matrix, any determines if any element of the matrix satisfies the predicate
any :: I.Elem a b => (a -> Bool) -> Matrix a b -> Bool
any f = VS.any (f . I.cast) . vals

-- | Returns the number of coefficients in a given matrix that evaluate to true
count :: I.Elem a b => (a -> Bool) -> Matrix a b -> Int
count f = VS.foldl' (\n x -> if f (I.cast x) then succ n else n) 0 . vals

{-| For vectors, the l2 norm, and for matrices the Frobenius norm.
    In both cases, it consists in the square root of the sum of the square of all the matrix entries.
    For vectors, this is also equals to the square root of the dot product of this with itself.
-}
norm :: I.Elem a b => Matrix a b -> a
norm = _prop I.norm

-- | For vectors, the squared l2 norm, and for matrices the Frobenius norm. In both cases, it consists in the sum of the square of all the matrix entries. For vectors, this is also equals to the dot product of this with itself.
squaredNorm :: I.Elem a b => Matrix a b -> a
squaredNorm = _prop I.squaredNorm

-- | The l2 norm of the matrix using the Blue's algorithm. A Portable Fortran Program to Find the Euclidean Norm of a Vector, ACM TOMS, Vol 4, Issue 1, 1978.
blueNorm :: I.Elem a b => Matrix a b -> a
blueNorm = _prop I.blueNorm

-- | The l2 norm of the matrix avoiding undeflow and overflow. This version use a concatenation of hypot calls, and it is very slow.
hypotNorm :: I.Elem a b => Matrix a b -> a
hypotNorm = _prop I.hypotNorm

-- | The determinant of the matrix
determinant :: I.Elem a b => Matrix a b -> a
determinant m
    | square m = _prop I.determinant m
    | otherwise = error "Matrix.determinant: non-square matrix"

-- | Adding two matrices by adding the corresponding entries together. You can use @(+)@ function as well.
add :: I.Elem a b => Matrix a b -> Matrix a b -> Matrix a b
add m1 m2
    | dims m1 == dims m2 = _binop const I.add m1 m2
    | otherwise = error "Matrix.add: matrices should have the same size"

-- | Subtracting two matrices by subtracting the corresponding entries together. You can use @(-)@ function as well.
sub :: I.Elem a b => Matrix a b -> Matrix a b -> Matrix a b
sub m1 m2
    | dims m1 == dims m2 = _binop const I.sub m1 m2
    | otherwise = error "Matrix.add: matrices should have the same size"

-- | Matrix multiplication. You can use @(*)@ function as well.
mul :: I.Elem a b => Matrix a b -> Matrix a b -> Matrix a b
mul m1 m2
    | cols m1 == rows m2 = _binop (\(rows, _) (_, cols) -> (rows, cols)) I.mul m1 m2
    | otherwise = error "Matrix.mul: number of columns for lhs matrix should be the same as number of rows for rhs matrix"

{- | Apply a given function to each element of the matrix.

Here is an example how to implement scalar matrix multiplication:

>>> let a = fromList [[1,2],[3,4]] :: MatrixXf

>>> a
Matrix 2x2
1.0 2.0
3.0 4.0

>>> map (*10) a
Matrix 2x2
10.0    20.0
30.0    40.0

-}
map :: I.Elem a b => (a -> a) -> Matrix a b -> Matrix a b
map f (Matrix rows cols vals) = Matrix rows cols (VS.map (I.cast . f . I.cast) vals)

{- | Apply a given function to each element of the matrix.

Here is an example how getting upper triangular matrix can be implemented:

>>> let a = fromList [[1,2,3],[4,5,6],[7,8,9]] :: MatrixXf

>>> a
Matrix 3x3
1.0 2.0 3.0
4.0 5.0 6.0
7.0 8.0 9.0

>>> imap (\row col val -> if row <= col then val else 0) a
Matrix 3x3
1.0 2.0 3.0
0.0 5.0 6.0
0.0 0.0 9.0

-}

imap :: I.Elem a b => (Int -> Int -> a -> a) -> Matrix a b -> Matrix a b
imap f (Matrix rows cols vals) = Matrix rows cols (VS.imap (\n -> let (c, r) = divMod n rows in I.cast . f r c . I.cast) vals)

-- | Upper trinagle of the matrix
upperTriangle :: I.Elem a b => Matrix a b -> Matrix a b
upperTriangle = imap (\row col val -> if row <= col then val else 0)

-- | Lower trinagle of the matrix
lowerTriangle :: I.Elem a b => Matrix a b -> Matrix a b
lowerTriangle = imap (\row col val -> if row >= col then val else 0)

-- | Filter elements in the matrix. Filtered elements will be replaced by 0
filter :: I.Elem a b => (a -> Bool) -> Matrix a b -> Matrix a b
filter f = map (\x -> if f x then x else 0)

-- | Filter elements in the matrix. Filtered elements will be replaced by 0
ifilter :: I.Elem a b => (Int -> Int -> a -> Bool) -> Matrix a b -> Matrix a b
ifilter f = imap (\r c x -> if f r c x then x else 0)

-- | Reduce matrix using user provided function applied to each element.
fold :: I.Elem a b => (c -> a -> c) -> c -> Matrix a b -> c
fold f a (Matrix _ _ vals) = VS.foldl (\a x -> f a (I.cast x)) a vals

-- | Reduce matrix using user provided function applied to each element. This is strict version of 'fold'
fold' :: I.Elem a b => (c -> a -> c) -> c -> Matrix a b -> c
fold' f a (Matrix _ _ vals) = VS.foldl' (\a x -> f a (I.cast x)) a vals

-- | Reduce matrix using user provided function applied to each element and it's index
ifold :: I.Elem a b => (Int -> Int -> c -> a -> c) -> c -> Matrix a b -> c
ifold f a (Matrix rows _ vals) = VS.ifoldl (\a n x -> let (c,r) = divMod n rows in f r c a (I.cast x)) a vals

-- | Reduce matrix using user provided function applied to each element and it's index. This is strict version of 'ifold'
ifold' :: I.Elem a b => (Int -> Int -> c -> a -> c) -> c -> Matrix a b -> c
ifold' f a (Matrix rows _ vals) = VS.ifoldl' (\a n x -> let (c,r) = divMod n rows in f r c a (I.cast x)) a vals

-- | Reduce matrix using user provided function applied to each element.
fold1 :: I.Elem a b => (a -> a -> a) -> Matrix a b -> a
fold1 f = foldl1 f . P.map I.cast . VS.toList . vals

-- | Reduce matrix using user provided function applied to each element. This is strict version of 'fold'
fold1' :: I.Elem a b => (a -> a -> a) -> Matrix a b -> a
fold1' f = L.foldl1' f . P.map I.cast . VS.toList . vals

-- | Diagonal of the matrix
diagonal :: I.Elem a b => Matrix a b -> Matrix a b
diagonal = _unop (\(rows, cols) -> (min rows cols, 1)) I.diagonal

{- | Inverse of the matrix

For small fixed sizes up to 4x4, this method uses cofactors. In the general case, this method uses PartialPivLU decomposition
-}
inverse :: I.Elem a b => Matrix a b -> Matrix a b
inverse m
    | square m = _unop id I.inverse m
    | otherwise = error "Matrix.inverse: non-square matrix"

-- | Adjoint of the matrix
adjoint :: I.Elem a b => Matrix a b -> Matrix a b
adjoint = _unop swap I.adjoint

-- | Transpose of the matrix
transpose :: I.Elem a b => Matrix a b -> Matrix a b
transpose = _unop swap I.transpose

-- | Conjugate of the matrix
conjugate :: I.Elem a b => Matrix a b -> Matrix a b
conjugate = _unop id I.conjugate

-- | Nomalize the matrix by deviding it on its 'norm'
normalize :: I.Elem a b => Matrix a b -> Matrix a b
normalize (Matrix rows cols vals) = I.performIO $ do
    vals <- VS.thaw vals
    VSM.unsafeWith vals $ \p ->
        I.call $ I.normalize p (I.cast rows) (I.cast cols)
    Matrix rows cols <$> VS.unsafeFreeze vals

-- | Apply a destructive operation to a matrix. The operation will be performed in place if it is safe to do so and will modify a copy of the matrix otherwise.
modify :: I.Elem a b => (forall s. M.MMatrix a b s -> ST s ()) -> Matrix a b -> Matrix a b
modify f (Matrix rows cols vals) = Matrix rows cols (VS.modify (f . M.MMatrix rows cols) vals)

-- | Convert matrix to different type using user provided element converter
convert :: (I.Elem a b, I.Elem c d) => (a -> c) -> Matrix a b -> Matrix c d
convert f (Matrix rows cols vals) = Matrix rows cols $ VS.map (I.cast . f . I.cast) vals

-- | Yield an immutable copy of the mutable matrix
freeze :: I.Elem a b => PrimMonad m => M.MMatrix a b (PrimState m) -> m (Matrix a b)
freeze (M.MMatrix mrows mcols mvals) = VS.freeze mvals >>= return . Matrix mrows mcols

-- | Yield a mutable copy of the immutable matrix
thaw :: I.Elem a b => PrimMonad m => Matrix a b -> m (M.MMatrix a b (PrimState m))
thaw (Matrix rows cols vals) = VS.thaw vals >>= return . M.MMatrix rows cols

-- | Unsafe convert a mutable matrix to an immutable one without copying. The mutable matrix may not be used after this operation.
unsafeFreeze :: I.Elem a b => PrimMonad m => M.MMatrix a b (PrimState m) -> m (Matrix a b)
unsafeFreeze (M.MMatrix mrows mcols mvals) = VS.unsafeFreeze mvals >>= return . Matrix mrows mcols

-- | Unsafely convert an immutable matrix to a mutable one without copying. The immutable matrix may not be used after this operation.
unsafeThaw :: I.Elem a b => PrimMonad m => Matrix a b -> m (M.MMatrix a b (PrimState m))
unsafeThaw (Matrix rows cols vals) = VS.unsafeThaw vals >>= return . M.MMatrix rows cols

-- | Pass a pointer to the matrix's data to the IO action. The data may not be modified through the pointer.
unsafeWith :: I.Elem a b => Matrix a b -> (Ptr b -> CInt -> CInt -> IO c) -> IO c
unsafeWith m@(Matrix rows cols vals) f
    | not (valid m) = fail "matrix layout is invalid"
    | otherwise = VS.unsafeWith vals $ \p -> f p (I.cast rows) (I.cast cols)

-- | Encode the matrix as a lazy byte string
encode :: I.Elem a b => Matrix a b -> BSL.ByteString
encode (Matrix rows cols vals) = BSL.fromChunks [
        I.encodeInt (I.code (VS.head vals)),
        I.encodeInt (I.cast rows),
        I.encodeInt (I.cast cols),
        let (fp, fs) = VS.unsafeToForeignPtr0 vals in BSI.PS (castForeignPtr fp) 0 (fs * sizeOf (VS.head vals))]

-- | Decode matrix from the lazy byte string
decode :: I.Elem a b => BSL.ByteString -> Matrix a b
decode st = Matrix rows cols vals where
    (rows, cols, vals) = runST $ do
        ref <- newSTRef st
        let next size = readSTRef ref >>= \a ->
                let (b,c) = BSL.splitAt (fromIntegral size) a
                in if BSL.length b /= fromIntegral size
                    then error "Matrix.decode: stream exhausted"
                    else do
                        writeSTRef ref c
                        return . BS.concat . BSL.toChunks $ b

        code <- I.decodeInt <$> next I.intSize
        when (code /= I.code (VS.head vals)) $
            fail "Matrix.decode: wrong matrix type"
        rows <- I.cast . I.decodeInt <$> next I.intSize
        cols <- I.cast . I.decodeInt <$> next I.intSize
        BSI.PS fp fo _ <- next (rows * cols * sizeOf (undefined :: CInt))
        BSL.null <$> readSTRef ref >>= (`unless` fail "Matrix.decode: stream underrun")
        return (rows, cols, VS.unsafeFromForeignPtr0 (I.plusForeignPtr fp fo) (rows * cols))


_prop :: I.Elem a b => (Ptr b -> Ptr b -> CInt -> CInt -> IO CString) -> Matrix a b -> a
_prop f m = I.cast $ I.performIO $ alloca $ \p -> do
    I.call $ unsafeWith m (f p)
    peek p

_binop :: I.Elem a b => ((Int, Int) -> (Int, Int) -> (Int, Int)) -> (Ptr b -> CInt -> CInt -> Ptr b -> CInt -> CInt -> Ptr b -> CInt -> CInt -> IO CString) -> Matrix a b -> Matrix a b -> Matrix a b
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

_unop :: I.Elem a b => ((Int,Int) -> (Int,Int)) -> (Ptr b -> CInt -> CInt -> Ptr b -> CInt -> CInt -> IO CString) -> Matrix a b -> Matrix a b
_unop f g m1 = I.performIO $ do
    m0 <- uncurry M.new $ f (dims m1)
    M.unsafeWith m0 $ \vals0 rows0 cols0 ->
        unsafeWith m1 $ \vals1 rows1 cols1 ->
            I.call $ g
                vals0 rows0 cols0
                vals1 rows1 cols1
    unsafeFreeze m0
