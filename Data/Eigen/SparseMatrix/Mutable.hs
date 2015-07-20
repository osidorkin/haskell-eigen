{-# LANGUAGE GADTs, RecordWildCards, ScopedTypeVariables #-}
module Data.Eigen.SparseMatrix.Mutable (
    -- * Mutable SparseMatrix
    IOSparseMatrix(..),
    IOSparseMatrixXf,
    IOSparseMatrixXd,
    IOSparseMatrixXcf,
    IOSparseMatrixXcd,
    new,
    reserve,
    -- * Matrix properties
    rows,
    cols,
    innerSize,
    outerSize,
    nonZeros,
    -- * Matrix compression
    compressed,
    compress,
    uncompress,
    -- * Accessing matrix data
    read,
    write,
    setZero,
    setIdentity,
    -- * Changing matrix shape
    resize,
    conservativeResize
) where

import Prelude hiding (read)
import Data.Complex
import Foreign.C.String
import Foreign.C.Types
import Foreign.ForeignPtr
import Foreign.Marshal.Alloc
import Foreign.Ptr
import Foreign.Storable
import qualified Foreign.Concurrent as FC
import qualified Data.Eigen.Internal as I

-- | Mutable version of sparse matrix. See `Data.Eigen.SparseMatrix.SparseMatrix` for details about matrix layout.
data IOSparseMatrix a b where
    IOSparseMatrix :: I.Elem a b => !(ForeignPtr (I.CSparseMatrix a b)) -> IOSparseMatrix a b

-- | Alias for single precision mutable matrix
type IOSparseMatrixXf = IOSparseMatrix Float CFloat
-- | Alias for double precision mutable matrix
type IOSparseMatrixXd = IOSparseMatrix Double CDouble
-- | Alias for single previsiom mutable matrix of complex numbers
type IOSparseMatrixXcf = IOSparseMatrix (Complex Float) (I.CComplex CFloat)
-- | Alias for double prevision mutable matrix of complex numbers
type IOSparseMatrixXcd = IOSparseMatrix (Complex Double) (I.CComplex CDouble)


-- | Creates new matrix with the given size @rows x cols@
new :: I.Elem a b => Int -> Int -> IO (IOSparseMatrix a b)
new rows cols = alloca $ \pm -> do
    I.call $ I.sparse_new (I.cast rows) (I.cast cols) pm
    m <- peek pm
    fm <- FC.newForeignPtr m $ I.call $ I.sparse_free m
    return $! IOSparseMatrix fm

-- | Returns the number of rows of the matrix
rows :: I.Elem a b => IOSparseMatrix a b -> IO Int
rows = _prop I.sparse_rows (return . I.cast)

-- | Returns the number of columns of the matrix
cols :: I.Elem a b => IOSparseMatrix a b  -> IO Int
cols = _prop I.sparse_cols (return . I.cast)

-- | Returns the number of rows (resp. columns) of the matrix if the storage order column major (resp. row major)
innerSize :: I.Elem a b => IOSparseMatrix a b  -> IO Int
innerSize = _prop I.sparse_innerSize (return . I.cast)

-- | Returns the number of columns (resp. rows) of the matrix if the storage order column major (resp. row major)
outerSize :: I.Elem a b => IOSparseMatrix a b  -> IO Int
outerSize = _prop I.sparse_outerSize (return . I.cast)

-- | Returns whether this matrix is in compressed form.
compressed :: I.Elem a b => IOSparseMatrix a b -> IO Bool
compressed = _prop I.sparse_isCompressed (return . (==1))

-- | Turns the matrix into the compressed format.
compress :: I.Elem a b => IOSparseMatrix a b -> IO ()
compress = _inplace I.sparse_compressInplace

-- | Turns the matrix into the uncompressed mode.
uncompress :: I.Elem a b => IOSparseMatrix a b -> IO ()
uncompress = _inplace I.sparse_uncompressInplace

-- | Reads the value of the matrix at position @i@, @j@.
-- This function returns @Scalar(0)@ if the element is an explicit zero.
read :: I.Elem a b => IOSparseMatrix a b -> Int -> Int -> IO a
read (IOSparseMatrix fm) row col = withForeignPtr fm $ \m -> alloca $ \px -> do
    I.call $ I.sparse_coeff m (I.cast row) (I.cast col) px
    I.cast <$> peek px

{- | Writes the value of the matrix at position @i@, @j@.
    This function turns the matrix into a non compressed form if that was not the case.

    This is a @O(log(nnz_j))@ operation (binary search) plus the cost of element insertion if the element does not already exist.
        
    Cost of element insertion is sorted insertion in O(1) if the elements of each inner vector are inserted in increasing inner index order, and in @O(nnz_j)@ for a random insertion.
-}
write :: I.Elem a b => IOSparseMatrix a b -> Int -> Int -> a -> IO ()
write (IOSparseMatrix fm) row col x = withForeignPtr fm $ \m -> alloca $ \px -> do
    I.call $ I.sparse_coeffRef m (I.cast row) (I.cast col) px
    peek px >>= (`poke` I.cast x)

-- | Sets the matrix to the identity matrix
setIdentity :: I.Elem a b => IOSparseMatrix a b -> IO ()
setIdentity = _inplace I.sparse_setIdentity

-- | Removes all non zeros but keep allocated memory
setZero :: I.Elem a b => IOSparseMatrix a b -> IO ()
setZero = _inplace I.sparse_setZero

-- | The number of non zero coefficients
nonZeros :: I.Elem a b => IOSparseMatrix a b -> IO Int
nonZeros = _prop I.sparse_nonZeros (return . I.cast)

-- | Preallocates space for non zeros. The matrix must be in compressed mode.
reserve :: I.Elem a b => IOSparseMatrix a b -> Int -> IO ()
reserve m s = _inplace (\p -> I.sparse_reserve p (I.cast s)) m

-- | Resizes the matrix to a rows x cols matrix and initializes it to zero.
resize :: I.Elem a b => IOSparseMatrix a b -> Int -> Int -> IO ()
resize m rows cols = _inplace (\p -> I.sparse_resize p (I.cast rows) (I.cast cols)) m

-- | Resizes the matrix to a rows x cols matrix leaving old values untouched.
conservativeResize :: I.Elem a b => IOSparseMatrix a b -> Int -> Int -> IO ()
conservativeResize m rows cols = _inplace (\p -> I.sparse_conservativeResize p (I.cast rows) (I.cast cols)) m

_inplace :: I.Elem a b => (Ptr (I.CSparseMatrix a b) -> IO CString) -> IOSparseMatrix a b -> IO ()
_inplace f (IOSparseMatrix fm) = withForeignPtr fm $ \m -> I.call $ f m

_prop :: Storable c => (I.CSparseMatrixPtr a b -> Ptr c -> IO CString) -> (c -> IO d) -> IOSparseMatrix a b -> IO d
_prop f g (IOSparseMatrix fp) =
    withForeignPtr fp $ \p ->
        alloca $ \pq -> do
            I.call (f p pq)
            peek pq >>= g

