{-# LANGUAGE CPP #-}
{-# LANGUAGE ForeignFunctionInterface #-}
{-# LANGUAGE RecordWildCards #-}

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
        a :: MatrixXd
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

You can compute relative error using @'norm' (ax - b) / 'norm' b@ formula or use 'relativeError' function which provides the same calculation implemented slightly more efficient.

-}

module Data.Eigen.LA (
    -- * Basic linear solving
    Decomposition(..),
    solve,
    relativeError,
    -- * Rank-revealing decompositions
    {- |
Certain decompositions are rank-revealing, i.e. are able to compute the 'rank' of a matrix. These are typically also the decompositions that behave best in the face of a non-full-rank matrix (which in the 'square' case means a singular matrix).

@
import Data.Eigen.Matrix
import Data.Eigen.LA

main = do
    let a = fromList [[1,2,5],[2,1,4],[3,0,3]] :: MatrixXd
    putStrLn "Here is the matrix A:" >> print a
    putStrLn "The rank of A is:" >> print (rank FullPivLU a)
    putStrLn "Here is a matrix whose columns form a basis of the null-space of A:" >> print (kernel FullPivLU a)
    putStrLn "Here is a matrix whose columns form a basis of the column-space of A:" >> print (image FullPivLU a)
@

produces the following output

@
Here is the matrix A:
Matrix 3x3
1.0 2.0 5.0
2.0 1.0 4.0
3.0 0.0 3.0

The rank of A is:
2
Here is a matrix whose columns form a basis of the null-space of A:
Matrix 3x1
0.5000000000000001
1.0
-0.5

Here is a matrix whose columns form a basis of the column-space of A:
Matrix 3x2
5.0 1.0
4.0 2.0
3.0 3.0
@
    -}
    rank,
    kernel,
    image,
    -- * Multiple linear regression
    {- | A linear regression model that contains more than one predictor variable. -}
    linearRegression
) where

import Prelude as P
import Foreign.Storable
import Foreign.Marshal.Alloc
import qualified Foreign.Concurrent as FC
#if __GLASGOW_HASKELL__ >= 710
#else
import Control.Applicative
#endif
import Data.Eigen.Matrix
import qualified Data.Eigen.Internal as I
import qualified Data.Eigen.Matrix.Mutable as M
import qualified Data.Vector.Storable as VS

{- |
@
Decomposition           Requirements on the matrix          Speed   Accuracy  Rank  Kernel  Image

PartialPivLU            Invertible                          ++      +         -     -       -
FullPivLU               None                                -       +++       +     +       +
HouseholderQR           None                                ++      +         -     -       -
ColPivHouseholderQR     None                                +       ++        +     -       -
FullPivHouseholderQR    None                                -       +++       +     -       -
LLT                     Positive definite                   +++     +         -     -       -
LDLT                    Positive or negative semidefinite   +++     ++        -     -       -
JacobiSVD               None                                -       +++       +     -       -
@
The best way to do least squares solving for square matrices is with a SVD decomposition ('JacobiSVD')
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
    | JacobiSVD deriving (Eq, Enum, Show, Read)


-- | [x = solve d a b] finds a solution @x@ of @ax = b@ equation using decomposition @d@
solve :: I.Elem a b => Decomposition -> Matrix a b -> Matrix a b -> Matrix a b
solve d a b = I.performIO $ do
    x <- M.new (cols a) 1
    M.unsafeWith x $ \x_vals x_rows x_cols ->
        unsafeWith a $ \a_vals a_rows a_cols ->
            unsafeWith b $ \b_vals b_rows b_cols ->
                I.call $ I.solve (I.cast $ fromEnum d)
                    x_vals x_rows x_cols
                    a_vals a_rows a_cols
                    b_vals b_rows b_cols
    unsafeFreeze x

-- | [e = relativeError x a b] computes @norm (ax - b) / norm b@ where @norm@ is L2 norm
relativeError :: I.Elem a b => Matrix a b -> Matrix a b -> Matrix a b -> a
relativeError x a b = I.performIO $
    unsafeWith x $ \x_vals x_rows x_cols ->
        unsafeWith a $ \a_vals a_rows a_cols ->
            unsafeWith b $ \b_vals b_rows b_cols ->
                alloca $ \pe -> do
                    I.call $ I.relativeError pe
                        x_vals x_rows x_cols
                        a_vals a_rows a_cols
                        b_vals b_rows b_cols
                    I.cast <$> peek pe

-- | The rank of the matrix
rank :: I.Elem a b => Decomposition -> Matrix a b -> Int
rank d m = I.performIO $ alloca $ \pr -> do
    I.call $ unsafeWith m $ I.rank (I.cast $ fromEnum d) pr
    I.cast <$> peek pr

-- | Return matrix whose columns form a basis of the null-space of @A@
kernel :: I.Elem a b => Decomposition -> Matrix a b -> Matrix a b
kernel d m1 = I.performIO $
    alloca $ \pvals ->
    alloca $ \prows ->
    alloca $ \pcols ->
        unsafeWith m1 $ \vals1 rows1 cols1 -> do
            I.call $ I.kernel (I.cast $ fromEnum d)
                pvals prows pcols
                vals1 rows1 cols1
            vals <- peek pvals
            rows <- I.cast <$> peek prows
            cols <- I.cast <$> peek pcols
            fp <- FC.newForeignPtr vals $ I.free vals
            return $ Matrix rows cols $ VS.unsafeFromForeignPtr0 fp $ rows * cols


-- | Return a matrix whose columns form a basis of the column-space of @A@
image :: I.Elem a b => Decomposition -> Matrix a b -> Matrix a b
image d m1 = I.performIO $
    alloca $ \pvals ->
    alloca $ \prows ->
    alloca $ \pcols ->
        unsafeWith m1 $ \vals1 rows1 cols1 -> do
            I.call $ I.image (I.cast $ fromEnum d)
                pvals prows pcols
                vals1 rows1 cols1
            vals <- peek pvals
            rows <- I.cast <$> peek prows
            cols <- I.cast <$> peek pcols
            fp <- FC.newForeignPtr vals $ I.free vals
            return $ Matrix rows cols $ VS.unsafeFromForeignPtr0 fp $ rows * cols


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
    a = fromList $ P.map ((1:).tail) points
    b = fromList $ P.map ((:[]).head) points
    x = solve ColPivHouseholderQR a b
    e = relativeError x a b
    coeffs = P.map head $ toList x





