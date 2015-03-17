{-# LANGUAGE RecordWildCards, ScopedTypeVariables #-}
module Data.Eigen.Matrix.Mutable (
    MMatrix(..),
    MMatrixXf,
    MMatrixXd,
    MMatrixXcf,
    MMatrixXcd,
    IOMatrix,
    STMatrix,
    -- * Construction
    new,
    replicate,
    -- * Consistency check
    valid,
    -- * Accessing individual elements
    read,
    write,
    unsafeRead,
    unsafeWrite,
    -- * Modifying matrices
    set,
    copy,
    unsafeCopy,
    -- * Raw pointers
    unsafeWith
) where

import Prelude hiding (read, replicate)
import Control.Monad.Primitive
import Foreign.Ptr
import Foreign.C.Types
import Data.Complex
import Text.Printf
import qualified Data.Vector.Storable.Mutable as VSM
import qualified Data.Eigen.Internal as I

-- | Mutable matrix. You can modify elements
data MMatrix a b s = MMatrix {
    mm_rows :: Int,
    mm_cols :: Int,
    mm_vals :: VSM.MVector s b
}

-- | Alias for single precision mutable matrix
type MMatrixXf = MMatrix Float CFloat
-- | Alias for double precision mutable matrix
type MMatrixXd = MMatrix Double CDouble
-- | Alias for single previsiom mutable matrix of complex numbers
type MMatrixXcf = MMatrix (Complex Float) (I.CComplex CFloat)
-- | Alias for double prevision mutable matrix of complex numbers
type MMatrixXcd = MMatrix (Complex Double) (I.CComplex CDouble)

type IOMatrix a b = MMatrix a b RealWorld
type STMatrix a b s = MMatrix a b s

-- | Verify matrix dimensions and memory layout
valid :: I.Elem a b => MMatrix a b s -> Bool
valid MMatrix{..} = mm_rows >= 0 && mm_cols >= 0 && VSM.length mm_vals == mm_rows * mm_cols

-- | Create a mutable matrix of the given size and fill it with 0 as an initial value.
new :: (PrimMonad m, I.Elem a b) => Int -> Int -> m (MMatrix a b (PrimState m))
new rows cols = replicate rows cols 0

-- | Create a mutable matrix of the given size and fill it with as an initial value.
replicate :: (PrimMonad m, I.Elem a b) => Int -> Int -> a -> m (MMatrix a b (PrimState m))
replicate rows cols val = do
    vals <- VSM.replicate (rows * cols) (I.cast val)
    return $ MMatrix rows cols vals

-- | Set all elements of the matrix to the given value
set :: (PrimMonad m, I.Elem a b) => (MMatrix a b (PrimState m)) -> a -> m ()
set MMatrix{..} val = VSM.set mm_vals (I.cast val)

-- | Copy a matrix. The two matrices must have the same size and may not overlap.
copy :: (PrimMonad m, I.Elem a b) => (MMatrix a b (PrimState m)) -> (MMatrix a b (PrimState m)) -> m ()
copy m1 m2
    | not (valid m1) = fail "MMatrix.copy: lhs matrix layout is invalid"
    | not (valid m2) = fail "MMatrix.copy: rhs matrix layout is invalid"
    | mm_rows m1 /= mm_rows m2 = fail "MMatrix.copy: matrices have different number of cols"
    | mm_cols m1 /= mm_cols m2 = fail "MMatrix.copy: matrices have different number of rows"
    | otherwise = VSM.copy (mm_vals m1) (mm_vals m2)

-- | Yield the element at the given position.
read :: (PrimMonad m, I.Elem a b) => MMatrix a b (PrimState m) -> Int -> Int -> m a
read mm@MMatrix{..} row col
    | not (valid mm) = fail "MMatrix.read: matrix layout is invalid"
    | row < 0 || row >= mm_rows = fail $ printf "MMatrix.read: row %d is out of bounds [0..%d)" row mm_rows
    | col < 0 || col >= mm_cols = fail $ printf "MMatrix.read: col %d is out of bounds [0..%d)" col mm_cols
    | otherwise = unsafeRead mm row col

-- | Replace the element at the given position.
write :: (PrimMonad m, I.Elem a b) => MMatrix a b (PrimState m) -> Int -> Int -> a -> m ()
write mm@MMatrix{..} row col val
    | not (valid mm) = fail "MMatrix.write: matrix layout is invalid"
    | row < 0 || row >= mm_rows = fail $ printf "MMatrix.write: row %d is out of bounds [0..%d)" row mm_rows
    | col < 0 || col >= mm_cols = fail $ printf "MMatrix.write: col %d is out of bounds [0..%d)" col mm_cols
    | otherwise = unsafeWrite mm row col val

-- | Copy a matrix. The two matrices must have the same size and may not overlap however no bounds check performaned to it may SEGFAULT for incorrect input.
unsafeCopy :: (PrimMonad m, I.Elem a b) => (MMatrix a b (PrimState m)) -> (MMatrix a b (PrimState m)) -> m ()
unsafeCopy m1 m2 = VSM.unsafeCopy (mm_vals m1) (mm_vals m2)

-- | Yield the element at the given position. No bounds checks are performed.
unsafeRead :: (PrimMonad m, I.Elem a b) => MMatrix a b (PrimState m) -> Int -> Int -> m a
unsafeRead MMatrix{..} row col = VSM.unsafeRead mm_vals (col * mm_rows + row) >>= \val -> return (I.cast val)

-- | Replace the element at the given position. No bounds checks are performed.
unsafeWrite :: (PrimMonad m, I.Elem a b) => MMatrix a b (PrimState m) -> Int -> Int -> a -> m ()
unsafeWrite MMatrix{..} row col val = VSM.unsafeWrite mm_vals (col * mm_rows + row) (I.cast val)

-- | Pass a pointer to the matrix's data to the IO action. Modifying data through the pointer is unsafe if the matrix could have been frozen before the modification.
unsafeWith :: I.Elem a b => IOMatrix a b -> (Ptr b -> CInt -> CInt -> IO c) -> IO c
unsafeWith mm@MMatrix{..} f
    | not (valid mm) = fail "mutable matrix layout is invalid"
    | otherwise = VSM.unsafeWith mm_vals $ \p -> f p (I.cast mm_rows) (I.cast mm_cols)

