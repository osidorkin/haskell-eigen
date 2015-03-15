{-# LANGUAGE RecordWildCards #-}
module Data.Eigen.Matrix.Mutable (
    MMatrix(..),
    IOMatrix,
    STMatrix,
    -- * Construction
    new,
    replicate,
    -- * Accessing individual elements
    read,
    write,
    unsafeRead,
    unsafeWrite,
    -- * Modifying matrices
    set,
    copy,
    -- * Raw pointers
    unsafeWith
) where

import Prelude hiding (read, replicate)
import Control.Monad.Primitive
import Foreign.Ptr
import Foreign.C.Types
import Text.Printf
import qualified Data.Vector.Storable.Mutable as VSM
import qualified Data.Eigen.Internal as I

-- | Mutable matrix. You can modify elements
data MMatrix s = MMatrix {
    mm_rows :: Int,
    mm_cols :: Int,
    mm_vals :: VSM.MVector s CDouble
}

type IOMatrix = MMatrix RealWorld
type STMatrix s = MMatrix s

-- | Verify matrix dimensions and memory layout
valid :: MMatrix s -> Bool
valid MMatrix{..} = mm_rows >= 0 && mm_cols >= 0 && VSM.length mm_vals == mm_rows * mm_cols

-- | Create a mutable matrix of the given size and fill it with 0 as an initial value.
new :: PrimMonad m => Int -> Int -> m (MMatrix (PrimState m))
new rows cols = do
    vals <- VSM.replicate (rows * cols) 0
    return $ MMatrix rows cols vals

-- | Create a mutable matrix of the given size and fill it with as an initial value.
replicate :: PrimMonad m => Int -> Int -> m (MMatrix (PrimState m))
replicate rows cols = do
    vals <- VSM.replicate (rows * cols) 0
    return $ MMatrix rows cols vals

-- | Set all elements of the matrix to the given value
set :: PrimMonad m => (MMatrix (PrimState m)) -> Double -> m ()
set MMatrix{..} val = VSM.set mm_vals (I.cast val)

-- | Copy a matrix. The two matrices must have the same length and may not overlap.
copy :: PrimMonad m => (MMatrix (PrimState m)) -> (MMatrix (PrimState m)) -> m ()
copy m1 m2 = VSM.copy (mm_vals m1) (mm_vals m2)

-- | Yield the element at the given position.
read :: PrimMonad m => MMatrix (PrimState m) -> Int -> Int -> m Double
read mm@MMatrix{..} row col
    | not (valid mm) = fail "MMatrix.read: matrix layout is invalid"
    | row < 0 || row >= mm_rows = fail $ printf "MMatrix.read: row %d is out of bounds [0..%d)" row mm_rows
    | col < 0 || col >= mm_cols = fail $ printf "MMatrix.read: col %d is out of bounds [0..%d)" col mm_cols
    | otherwise = unsafeRead mm row col

-- | Replace the element at the given position.
write :: PrimMonad m => MMatrix (PrimState m) -> Int -> Int -> Double -> m ()
write mm@MMatrix{..} row col val
    | not (valid mm) = fail "MMatrix.write: matrix layout is invalid"
    | row < 0 || row >= mm_rows = fail $ printf "MMatrix.write: row %d is out of bounds [0..%d)" row mm_rows
    | col < 0 || col >= mm_cols = fail $ printf "MMatrix.write: col %d is out of bounds [0..%d)" col mm_cols
    | otherwise = unsafeWrite mm row col val

-- | Yield the element at the given position. No bounds checks are performed.
unsafeRead :: PrimMonad m => MMatrix (PrimState m) -> Int -> Int -> m Double
unsafeRead MMatrix{..} row col = VSM.unsafeRead mm_vals (col * mm_rows + row) >>= \val -> return (I.cast val)

-- | Replace the element at the given position. No bounds checks are performed.
unsafeWrite :: PrimMonad m => MMatrix (PrimState m) -> Int -> Int -> Double -> m ()
unsafeWrite MMatrix{..} row col val = VSM.unsafeWrite mm_vals (col * mm_rows + row) (I.cast val)

-- | Pass a pointer to the matrix's data to the IO action. Modifying data through the pointer is unsafe if the matrix could have been frozen before the modification.
unsafeWith :: IOMatrix -> (CInt -> CInt -> Ptr CDouble -> IO a) -> IO a
unsafeWith mm@MMatrix{..} f
    | not (valid mm) = fail "mutable matrix layout is invalid"
    | otherwise = VSM.unsafeWith mm_vals $ \p -> f (I.cast mm_rows) (I.cast mm_cols) p

