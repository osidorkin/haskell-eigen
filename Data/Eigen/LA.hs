{-# LANGUAGE RecordWildCards, ForeignFunctionInterface #-}

{- |

/This description is a based on/ <http://eigen.tuxfamily.org/dox/TutorialLinearAlgebra.html> /tutorial/

Basic linear solving

The problem: You have a system of equations, that you have written as a single matrix equation

@Ax = b@

Where A and b are matrices (b could be a vector, as a special case). You want to find a solution x.

The solution: You can choose between various decompositions, depending on what your matrix A looks like, and depending on whether you favor speed or accuracy. However, let's start with an example that works in all cases, and is a good compromise:

@
import "Data.Eigen.Matrix"
import "Data.Eigen.LA"

main = do
    let
        a = fromList [[1,2,3], [4,5,6], [7,8,10]]
        b = fromList [[3],[3],[4]]
        x = solve ColPivHouseholderQR a b
    putStrLn \"Here is the matrix A:\"
    print a

    putStrLn \"Here is the vector b:\"
    print b

    putStrLn \"The solution is:\"
    print x
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

You can compute relative error using @norm (ax - b) / norm b@ formula or use 'relativeError' function which provides the same but slightly more efficient calculation.

-}

module Data.Eigen.LA (
    Decomposition(..),
    solve,
    relativeError,
    linearRegression
) where

import Foreign.Ptr
import Foreign.C.Types
import Foreign.C.String
import Foreign.Storable
import Foreign.Marshal.Alloc
import Data.Eigen.Matrix
import qualified Data.Eigen.Matrix.Mutable as MM
import Data.Eigen.Internal
import System.IO.Unsafe

foreign import ccall "eigen-proxy.h eigen_solve" c_solve :: CInt -> Ptr C_MatrixXd -> Ptr C_MatrixXd -> Ptr C_MatrixXd -> IO CString
foreign import ccall "eigen-proxy.h eigen_relativeError" c_relativeError :: Ptr CDouble -> Ptr C_MatrixXd -> Ptr C_MatrixXd -> Ptr C_MatrixXd -> IO CString



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

The best way to do least squares solving for square matrixes is with a SVD decomposition (JacobiSVD)
@
-}

data Decomposition
    -- | LU decomposition of a matrix with partial pivoting, and related features
    = PartialPivLU
    -- | LU decomposition of a matrix with complete pivoting, and related features
    | FullPivLU
    -- | Householder QR decomposition of a matrix.
    | HouseholderQR
    -- | Householder rank-revealing QR decomposition of a matrix with column-pivoting.
    | ColPivHouseholderQR
    -- | Householder rank-revealing QR decomposition of a matrix with full pivoting.
    | FullPivHouseholderQR
    -- | Standard Cholesky decomposition (LL^T) of a matrix and associated features.
    | LLT
    -- | Robust Cholesky decomposition of a matrix with pivoting.
    | LDLT
    -- | Two-sided Jacobi SVD decomposition of a rectangular matrix.
    | JacobiSVD deriving (Show, Enum)

-- | /solve d a b/ finds a solution x of @ax = b@ equation using decomposition @d@
solve :: Decomposition -> Matrix -> Matrix -> Matrix
solve d a b = (`modify` empty) $ \x ->
    with a $ \pa ->
    with b $ \pb ->
    MM.with x $ \px ->
        call $ c_solve (fromIntegral $ fromEnum d) px pa pb

-- | /relativeError x a b/ computes @norm (ax - b) / norm b@ where norm is L2 norm
relativeError :: Matrix -> Matrix -> Matrix -> Double
relativeError x a b = unsafePerformIO $
    with x $ \px ->
    with a $ \pa ->
    with b $ \pb ->
    alloca $ \pr -> do
        call $ c_relativeError pr px pa pb
        fmap cast $ peek pr

{- | compute multiple linear regression @y = a1 x1 + a2 x2 + ... + an xn + b@ using 'ColPivHouseholderQR' decomposition

    * argument is a list of point in format @[y, x1..xn]@

    * return value is (coeffs, relative error)

    * coeffs format is @[b, a1..an]@
-}
linearRegression :: [[Double]] -> ([Double], Double)
linearRegression points = (coeffs, e) where
    a = fromList $ map ((1:).tail) points
    b = fromList $ map ((:[]).head) points
    x = solve ColPivHouseholderQR a b
    e = relativeError x a b
    coeffs = map head $ toList x





