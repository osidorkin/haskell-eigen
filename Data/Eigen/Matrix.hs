{-# LANGUAGE RecordWildCards, MultiParamTypeClasses #-}
module Data.Eigen.Matrix (
    -- * Matrix type
    Matrix(..),
    -- * Matrix conversions
    fromList,
    toList,
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
    -- * Matrix transformations
    inverse,
    adjoint,
    transpose,
    normalize,
    -- * Mutable operations
    freeze,
    thaw,
    modify,
    with
) where

import Data.List (intercalate)
import Text.Printf
import Foreign.Ptr
import Foreign.Storable
import Foreign.ForeignPtr
import Foreign.C.Types
import Foreign.Marshal.Array
import Control.Monad
import Control.Monad.ST
import System.IO.Unsafe
import qualified Data.Vector.Storable as VS
import qualified Data.Vector.Storable.Mutable as VSM
import qualified Data.Eigen.Matrix.Mutable as MM
import Data.Eigen.Internal

-- | constant Matrix class to be used in pure computations, uses the same column major memory layout as Eigen MatrixXd
data Matrix = Matrix {
    m_rows :: Int,
    m_cols :: Int,
    m_vals :: VS.Vector CDouble --
};

-- | pretty prints the matrix
instance Show Matrix where
    show m@Matrix{..} = printf "Matrix %dx%d\n%s\n" m_rows m_cols $
        intercalate "\n" $ map (intercalate "\t" . map show) $ toList m

-- | only the following functions are defined for Num instance: (*), (+), (-)
instance Num Matrix where
    (*) = binop MM.mul
    (+) = binop MM.add
    (-) = binop MM.sub
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
identity size = Matrix size size $ runST $ do
    vm <- VSM.replicate (size * size) 0
    forM_ [0..pred size] $ \n ->
        VSM.write vm (n * size + n) 1
    VS.unsafeFreeze vm

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
block startRow startCol blockRows blockCols m = fromList $
    [[coeff row col m | col <- take blockCols [startCol..]] | row <- take blockRows [startRow..]]

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

binop :: (MM.MMatrix -> MM.MMatrix -> MM.MMatrix -> IO a) -> Matrix -> Matrix -> Matrix
binop f lhs rhs = unsafePerformIO $ do
    ret <- MM.new 0 0
    lhs <- thaw lhs
    rhs <- thaw rhs
    _ <- f ret lhs rhs
    freeze ret

-- | construct matrix from a list of rows, column count is detected as maximum row length
fromList :: [[Double]] -> Matrix
fromList list = Matrix rows cols vals where
    rows = length list
    cols = maximum $ map length list
    vals = runST $ do
        vm <- VSM.replicate (rows * cols) 0
        zipWithM_ (\row vals ->
            zipWithM_ (\col val ->
                VSM.write vm (col * rows + row) (cast val)) [0..] vals) [0..] list
        VS.unsafeFreeze vm

-- | converts matrix to a list of its rows
toList :: Matrix -> [[Double]]
toList Matrix{..} = [[cast $ m_vals VS.! (col * m_rows + row) | col <- [0..pred m_cols]] | row <- [0..pred m_rows]]

-- | for vectors, the l2 norm, and for matrices the Frobenius norm. In both cases, it consists in the square root of the sum of the square of all the matrix entries. For vectors, this is also equals to the square root of the dot product of this with itself.
norm :: Matrix -> Double
norm m = unsafePerformIO $ thaw m >>= MM.norm

-- | for vectors, the squared l2 norm, and for matrices the Frobenius norm. In both cases, it consists in the sum of the square of all the matrix entries. For vectors, this is also equals to the dot product of this with itself.
squaredNorm :: Matrix -> Double
squaredNorm m = unsafePerformIO $ thaw m >>= MM.squaredNorm

-- | the determinant of the matrix
determinant :: Matrix -> Double
determinant m = unsafePerformIO $ thaw m >>= MM.determinant

{- | inverse of the matrix

For small fixed sizes up to 4x4, this method uses cofactors. In the general case, this method uses class 'PartialPivLU'
-}
inverse :: Matrix -> Matrix
inverse = modify MM.inverse

-- | adjoint of the matrix
adjoint :: Matrix -> Matrix
adjoint = modify MM.adjoint

-- | transpose of the matrix
transpose :: Matrix -> Matrix
transpose = modify MM.transpose

-- | nomalize the matrix by deviding it on its 'norm'
normalize :: Matrix -> Matrix
normalize = modify MM.normalize

-- | create a snapshot of mutable matrix
freeze :: MM.MMatrix -> IO Matrix
freeze mm = MM.with mm $ \pm -> do
    rows <- fmap cast $ c_rows pm
    cols <- fmap cast $ c_cols pm
    let len = rows * cols
    src <- c_data pm
    fp <- mallocForeignPtrArray len
    withForeignPtr fp $ \dst -> copyArray dst src len
    Matrix {
        m_rows = rows,
        m_cols = cols,
        m_vals = VS.unsafeFromForeignPtr fp 0 len
    }

-- | create mutable copy of the matrix
thaw :: Matrix -> IO MM.MMatrix
thaw Matrix{..} = withForeignPtr fp $ \src -> do
    let len = m_rows * m_cols
    pm <- c_create (cast m_rows) (cast m_cols)
    dst <- c_data pm
    copyArray dst (plusPtr src $ off * sizeOf (undefined :: CDouble)) len
    fmap MM.MMatrix $ newForeignPtr c_destroy pm
    where (fp, off, _) = VS.unsafeToForeignPtr m_vals

-- | apply mutable operation to the mutable copy of the matrix and snapshot of this copy
modify :: (MM.MMatrix -> IO ()) -> Matrix -> Matrix
modify f m = unsafePerformIO $ thaw m >>= \mm -> f mm >> freeze mm

-- | apply foreign operation to the mutable copy of the matrix and operation result
with :: Matrix -> (Ptr C_MatrixXd -> IO a) -> IO a
with m f = thaw m >>= \mm -> MM.with mm f
