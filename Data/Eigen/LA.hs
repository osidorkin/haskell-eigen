{-# LANGUAGE RecordWildCards, ForeignFunctionInterface #-}

{- |

The problem: You have a system of equations, that you have written as a single matrix equation

@Ax = b@

Where A and b are matrices (b could be a vector, as a special case). You want to find a solution x.

The solution: You can choose between various decompositions, depending on what your matrix A looks like, and depending on whether you favor speed or accuracy. However, let's start with an example that works in all cases, and is a good compromise:

@
import Data.Eigen.Matrix
import Data.Eigen.LA

main = do
    let
        a = fromList [[1,2,3], [4,5,6], [7,8,10]]
        b = fromList [[3],[3],[4]]
        x = solve ColPivHouseholderQR a b
    putStrLn \"Here is the matrix A:\" >> print a
    putStrLn \"Here is the vector b:\" >> print b
    putStrLn \"The solution is:\" >> print x
@

produces the following output

@
Here is the matrix A:
Matrix 3x3
1.0 2.0 3.0
4.0 5.0 6.0
7.0 8.0 10.0

Here is the vector b:
Matrix 3x1
3.0
3.0
4.0

The solution is:
Matrix 3x1
-2.0000000000000004
1.0000000000000018
0.9999999999999989
@

Checking if a solution really exists: Only you know what error margin you want to allow for a solution to be considered valid.

You can compute relative error using @norm (ax - b) / norm b@ formula or use 'relativeError' function which provides the same calculation implemented slightly more efficient.

-}

module Data.Eigen.LA (
    -- * Basic linear solving
    Decomposition(..),
    solve,
    relativeError,
    -- * Multiple linear regression
    linearRegression
) where

import Foreign.Ptr
import Foreign.C.Types
import Foreign.C.String
import Foreign.Storable
import Foreign.Marshal.Alloc
import Control.Applicative
import Data.Eigen.Matrix
import Data.Eigen.Internal
import qualified Data.Vector.Storable as VS
import qualified Data.Vector.Storable.Mutable as VSM

foreign import ccall "eigen-proxy.h eigen_solve" c_solve :: CInt -> Ptr CDouble -> CInt -> CInt -> Ptr CDouble -> CInt -> CInt -> Ptr CDouble -> CInt -> CInt -> IO CString
foreign import ccall "eigen-proxy.h eigen_relativeError" c_relativeError :: Ptr CDouble -> Ptr CDouble -> CInt -> CInt -> Ptr CDouble -> CInt -> CInt -> Ptr CDouble -> CInt -> CInt -> IO CString



{- |
@
Decomposition           Requirements on the matrix          Speed   Accuracy

PartialPivLU            Invertible                          ++      +
FullPivLU               None                                -       +++
HouseholderQR           None                                ++      +
ColPivHouseholderQR     None                                +       ++
FullPivHouseholderQR    None                                -       +++
LLT                     Positive definite                   +++     +
LDLT                    Positive or negative semidefinite   +++     ++
JacobiSVD               None                                -       +++

The best way to do least squares solving for square matrices is with a SVD decomposition (JacobiSVD)
@
-}

data Decomposition
    -- | LU decomposition of a matrix with partial pivoting.
    = PartialPivLU
    -- | LU decomposition of a matrix with complete pivoting.
    | FullPivLU
    -- | Householder QR decomposition of a matrix.
    | HouseholderQR
    -- | Householder rank-revealing QR decomposition of a matrix with column-pivoting.
    | ColPivHouseholderQR
    -- | Householder rank-revealing QR decomposition of a matrix with full pivoting.
    | FullPivHouseholderQR
    -- | Standard Cholesky decomposition (LL^T) of a matrix.
    | LLT
    -- | Robust Cholesky decomposition of a matrix with pivoting.
    | LDLT
    -- | Two-sided Jacobi SVD decomposition of a rectangular matrix.
    | JacobiSVD deriving (Show, Enum)

-- | [x = solve d a b] finds a solution @x@ of @ax = b@ equation using decomposition @d@
solve :: Decomposition -> Matrix -> Matrix -> Matrix
solve d a b = performIO $ do
    let
        cols = 1
        rows = m_cols a
    vals <- VSM.new (rows * cols)
    VSM.unsafeWith vals $ \px ->
        VS.unsafeWith (m_vals a) $ \pa ->
            VS.unsafeWith (m_vals b) $ \pb ->
                call $ c_solve (cast $ fromEnum d)
                    px (cast rows) (cast cols)
                    pa (cast $ m_rows a) (cast $ m_cols a)
                    pb (cast $ m_rows b) (cast $ m_cols b)
    Matrix rows cols <$> VS.unsafeFreeze vals


-- | [e = relativeError x a b] computes @norm (ax - b) / norm b@ where @norm@ is L2 norm
relativeError :: Matrix -> Matrix -> Matrix -> Double
relativeError x a b = performIO $ do
    VS.unsafeWith (m_vals x) $ \px ->
        VS.unsafeWith (m_vals a) $ \pa ->
            VS.unsafeWith (m_vals b) $ \pb ->
                alloca $ \pe -> do
                    call $ c_relativeError pe
                        px (cast $ m_rows x) (cast $ m_cols x)
                        pa (cast $ m_rows a) (cast $ m_cols a)
                        pb (cast $ m_rows b) (cast $ m_cols b)
                    cast <$> peek pe

{- |
[(coeffs, error) = linearRegression points] computes multiple linear regression @y = a1 x1 + a2 x2 + ... + an xn + b@ using 'ColPivHouseholderQR' decomposition

* point format is @[y, x1..xn]@

* coeffs format is @[b, a1..an]@

* error is calculated using 'relativeError'

@
import Data.Eigen.LA
main = print $ linearRegression [
    [-4.32, 3.02, 6.89],
    [-3.79, 2.01, 5.39],
    [-4.01, 2.41, 6.01],
    [-3.86, 2.09, 5.55],
    [-4.10, 2.58, 6.32]]
 @

 produces the following output

 @
 ([-2.3466569233817127,-0.2534897541434826,-0.1749653335680988],1.8905965120153139e-3)
 @

-}
linearRegression :: [[Double]] -> ([Double], Double)
linearRegression points = (coeffs, e) where
    a = fromList $ map ((1:).tail) points
    b = fromList $ map ((:[]).head) points
    x = solve ColPivHouseholderQR a b
    e = relativeError x a b
    coeffs = map head $ toList x





