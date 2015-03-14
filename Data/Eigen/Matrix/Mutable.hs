{-# LANGUAGE RecordWildCards #-}
module Data.Eigen.Matrix.Mutable where

import Control.Monad.Primitive
import Foreign.C.Types
import qualified Data.Vector.Storable.Mutable as VSM
import Data.Eigen.Internal

-- | Mutable matrix. You can modify elements
data MMatrix s = MMatrix {
    mm_rows :: Int,
    mm_cols :: Int,
    mm_vals :: VSM.MVector s CDouble
}

type IOMatrix = MMatrix RealWorld
type STMatrix s = MMatrix s

-- | Creates a mutable matrix of the given dimension. Elements are initialized with 0.
new :: PrimMonad m => Int -> Int -> m (MMatrix (PrimState m))
new rows cols = do
    vals <- VSM.replicate (rows * cols) 0
    return $ MMatrix rows cols vals

-- | Set all elements of the matrix to the given value
set :: PrimMonad m => (MMatrix (PrimState m)) -> Double -> m ()
set MMatrix{..} val = VSM.set mm_vals (cast val)

-- | Yield the element at the given position.
read :: PrimMonad m => MMatrix (PrimState m) -> Int -> Int -> m Double
read MMatrix{..} row col
    | row >= 0 && row < mm_rows && 
      col >= 0 && col < mm_cols = do
        val <- VSM.read mm_vals (col * mm_rows + row)
        return $ cast val
    | otherwise = fail "MMatrix.read: wrong index"

-- | Replace the element at the given position.
write :: PrimMonad m => MMatrix (PrimState m) -> Int -> Int -> Double -> m ()
write MMatrix{..} row col val
    | row >= 0 && row < mm_rows && 
      col >= 0 && col < mm_cols = 
        VSM.write mm_vals (col * mm_rows + row) (cast val)
    | otherwise = fail "MMatrix.write: wrong index"


