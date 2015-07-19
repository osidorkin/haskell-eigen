{-# LANGUAGE RecordWildCards #-}
{-# LANGUAGE ForeignFunctionInterface #-}
{-# LANGUAGE MultiParamTypeClasses #-}
{-# LANGUAGE FunctionalDependencies #-}
{-# LANGUAGE FlexibleInstances #-}

{- |

This documentation is based on original Eigen page <http://eigen.tuxfamily.org/dox/group__TopicSparseSystems.html Solving Sparse Linear Systems>

Eigen currently provides a limited set of built-in MPL2 compatible solvers.
They are summarized in the following table:

@
Sparse solver       Solver kind             Matrix kind         Notes

ConjugateGradient   Classic iterative CG    SPD                 Recommended for large symmetric
                                                                problems (e.g., 3D Poisson eq.)
BiCGSTAB            Iterative stabilized    Square
                    bi-conjugate gradient
SparseLU            LU factorization        Square              Optimized for small and large problems
                                                                with irregular patterns
SparseQR            QR factorization        Any, rectangular    Recommended for least-square problems,
                                                                has a basic rank-revealing feature
@

All these solvers follow the same general concept. Here is a typical and general example:

@
let
    a :: SparseMatrixXd
    a = ... -- fill a

    b :: SparseMatrixXd
    b = ... -- fill b

    validate msg = info >>= (`when` fail msg) . (/= Success)

// solve Ax = b
runSolverT solver $ do
    compute a
    validate "decomposition failed"

    x <- solve b
    validate "solving failed"

    // solve for another right hand side
    x1 <- solve b1
@

In the case where multiple problems with the same sparsity pattern have to be solved, then the "compute" step can be decomposed as follow:

@
runSolverT solver $ do
    analyzePattern a1
    factorize a1
    x1 <- solve b1
    x2 <- solve b2

    factorize a2
    x1 <- solve b1
    x2 <- solve b2
@

Finally, each solver provides some specific features, such as determinant, access to the factors, controls of the iterations, and so on.

= The Compute Step

In the `compute` function, the matrix is generally factorized: LLT for self-adjoint matrices, LDLT for general hermitian matrices,
LU for non hermitian matrices and QR for rectangular matrices. These are the results of using direct solvers.
For this class of solvers precisely, the compute step is further subdivided into `analyzePattern` and `factorize`.

The goal of `analyzePattern` is to reorder the nonzero elements of the matrix, such that the factorization step creates less fill-in.
This step exploits only the structure of the matrix. Hence, the results of this step can be used for other linear systems where the
matrix has the same structure.

In `factorize`, the factors of the coefficient matrix are computed. This step should be called each time the values of the matrix change.
However, the structural pattern of the matrix should not change between multiple calls.

For iterative solvers, the `compute` step is used to eventually setup a preconditioner.
Remember that, basically, the goal of the preconditioner is to speedup the convergence of an iterative method by solving a modified linear
system where the coefficient matrix has more clustered eigenvalues.
For real problems, an iterative solver should always be used with a preconditioner.

= The Solve Step

The `solve` function computes the solution of the linear systems with one or many right hand sides.

@
x <- solve b
@

Here, @b@ can be a vector or a matrix where the columns form the different right hand sides.
The `solve` function can be called several times as well, for instance when all the right hand sides are not available at once.

@
x1 <- solve b1
-- Get the second right hand side b2
x2 <- solve b2
--  ...
@

For direct methods, the solution are computed at the machine precision. Sometimes, the solution need not be too accurate.
In this case, the iterative methods are more suitable and the desired accuracy can be set before the solve step using `setTolerance`.
-}

module Data.Eigen.SparseLA (
    SolverInfo(..),
    ComputationInfo(..),
    SolverT,
    runSolverT,
    -- * The Compute step
    analyzePattern,
    factorize,
    compute,
    -- * The Solve step
    solve,
    tolerance,
    setTolerance,
    maxIterations,
    setMaxIterations,
    --solveWithGuess,
    info,
    Data.Eigen.SparseLA.error,
    iterations,
) where

import Prelude as P
import Foreign.Ptr
import Foreign.ForeignPtr
import Foreign.Storable
import Foreign.C.String
import Foreign.Marshal.Alloc
import Control.Monad.IO.Class
import Control.Monad.Trans.Reader
import qualified Foreign.Concurrent as FC
import Control.Applicative
import qualified Data.Eigen.Internal as I
import qualified Data.Eigen.SparseMatrix as SM

data SolverInfo
    {- | A conjugate gradient solver for sparse self-adjoint problems.

        This class allows to solve for @A.x = b@ sparse linear problems using a conjugate gradient algorithm. The sparse matrix @A@ must be selfadjoint.

        The maximal number of iterations and tolerance value can be controlled via the `setMaxIterations` and `setTolerance` methods.
        The defaults are the size of the problem for the maximal number of iterations and @epsilon@ for the tolerance
    -}
    = ConjugateGradient
    {- | A bi conjugate gradient stabilized solver for sparse square problems.

        This class allows to solve for @A.x = b@ sparse linear problems using a bi conjugate gradient stabilized algorithm.
        The vectors @x@ and @b@ can be either dense or sparse.

        The maximal number of iterations and tolerance value can be controlled via the `setMaxIterations` and `setTolerance` methods.
        The defaults are the size of the problem for the maximal number of iterations and @epsilon@ for the tolerance
    -}
    | BiCGSTAB
    {- | Sparse supernodal LU factorization for general matrices.

        This class implements the supernodal LU factorization for general matrices. It uses the main techniques from the sequential
        <http://crd-legacy.lbl.gov/~xiaoye/SuperLU/ SuperLU package>. It handles transparently real and complex arithmetics with
        single and double precision, depending on the scalar type of your input matrix. The code has been optimized to provide BLAS-3
        operations during supernode-panel updates. It benefits directly from the built-in high-performant Eigen BLAS routines.
        Moreover, when the size of a supernode is very small, the BLAS calls are avoided to enable a better optimization from the compiler.
        For best performance, you should compile it with NDEBUG flag to avoid the numerous bounds checking on vectors.
    -}
    | SparseLU
    {- | Sparse left-looking rank-revealing QR factorization.

        This class implements a left-looking rank-revealing QR decomposition of sparse matrices. When a column has a norm less than a given
        tolerance it is implicitly permuted to the end. The QR factorization thus obtained is given by @A*P = Q*R@ where @R@ is upper triangular or trapezoidal.

        @P@ is the column permutation which is the product of the fill-reducing and the rank-revealing permutations.

        @Q@ is the orthogonal matrix represented as products of Householder reflectors.

        @R@ is the sparse triangular or trapezoidal matrix. The later occurs when @A@ is rank-deficient.
    -}
    | SparseQR
    deriving (Eq, Enum, Show, Read)

data ComputationInfo
    -- | Computation was successful.
    = Success
    -- | The provided data did not satisfy the prerequisites.
    | NumericalIssue
    -- | Iterative procedure did not converge.
    | NoConvergence
    -- | The inputs are invalid, or the algorithm has been improperly called. When assertions are enabled, such errors trigger an error.
    | InvalidInput
    deriving (Eq, Enum, Show, Read)

type SolverT a b m = ReaderT (SolverInfo, ForeignPtr (I.CSolver a b)) m

runSolverT :: (MonadIO m, I.Elem a b) => SolverInfo -> SolverT a b m c -> m c
runSolverT i f = do
    fs <- liftIO $ alloca $ \ps -> do
        I.call $ I.sparse_la_newSolver i ps
        s <- peek ps
        FC.newForeignPtr s (I.call $ I.sparse_la_freeSolver i s)
    runReaderT f (i,fs)

-- | Initializes the iterative solver for the sparsity pattern of the matrix @A@ for further solving @Ax=b@ problems.
analyzePattern :: (MonadIO m, I.Elem a b) => SM.SparseMatrix a b -> SolverT a b m ()
analyzePattern (SM.SparseMatrix fa) = ask >>= \(i,fs) -> liftIO $
    withForeignPtr fs $ \s ->
    withForeignPtr fa $ \a ->
        I.call $ I.sparse_la_analyzePattern i s a

-- | nitializes the iterative solver with the numerical values of the matrix @A@ for further solving @Ax=b@ problems.
factorize :: (MonadIO m, I.Elem a b) => SM.SparseMatrix a b -> SolverT a b m ()
factorize (SM.SparseMatrix fa) = ask >>= \(i,fs) -> liftIO $
    withForeignPtr fs $ \s ->
    withForeignPtr fa $ \a ->
        I.call $ I.sparse_la_factorize i s a

-- | Initializes the iterative solver with the matrix @A@ for further solving @Ax=b@ problems.
--
-- The `compute` method is equivalent to calling both `analyzePattern` and `factorize`.
compute :: (MonadIO m, I.Elem a b) => SM.SparseMatrix a b -> SolverT a b m ()
compute (SM.SparseMatrix fa) = ask >>= \(i,fs) -> liftIO $
    withForeignPtr fs $ \s ->
    withForeignPtr fa $ \a ->
        I.call $ I.sparse_la_compute i s a

-- | The tolerance threshold used by the stopping criteria.
tolerance :: (MonadIO m, I.Elem a b) => SolverT a b m Double
tolerance = _get_prop I.sparse_la_tolerance

-- | Sets the tolerance threshold used by the stopping criteria.
-- | This value is used as an upper bound to the relative residual error: @|Ax-b|/|b|@. The default value is the machine precision given by epsilon
setTolerance :: (MonadIO m, I.Elem a b) => Double -> SolverT a b m ()
setTolerance = _set_prop I.sparse_la_setTolerance

-- | The max number of iterations. It is either the value setted by setMaxIterations or, by default, twice the number of columns of the matrix.
maxIterations :: (MonadIO m, I.Elem a b) => SolverT a b m Int
maxIterations = _get_prop I.sparse_la_maxIterations

-- | Sets the max number of iterations. Default is twice the number of columns of the matrix.
setMaxIterations :: (MonadIO m, I.Elem a b) => Int -> SolverT a b m ()
setMaxIterations = _set_prop I.sparse_la_setMaxIterations

-- | An expression of the solution x of @A x = b@ using the current decomposition of @A@.
solve :: (MonadIO m, I.Elem a b) => SM.SparseMatrix a b -> SolverT a b m (SM.SparseMatrix a b)
solve (SM.SparseMatrix fb) = ask >>= \(i,fs) -> liftIO $
    withForeignPtr fs $ \s ->
    withForeignPtr fb $ \b ->
    alloca $ \px -> do
        I.call $ I.sparse_la_solve i s b px
        x <- peek px
        SM.SparseMatrix <$> FC.newForeignPtr x (I.call $ I.sparse_free x)

{-
-- | The solution @x@ of @A x = b@ using the current decomposition of @A@ and @x0@ as an initial solution.
solveWithGuess :: (MonadIO m, I.Elem a b) => SM.SparseMatrix a b -> SM.SparseMatrix a b -> SolverT a b m (SM.SparseMatrix a b)
solveWithGuess (SM.SparseMatrix fb) (SM.SparseMatrix fx0) = ask >>= \(i,fs) -> liftIO $
    withForeignPtr fs $ \s ->
    withForeignPtr fb $ \b ->
    withForeignPtr fx0 $ \x0 ->
    alloca $ \px -> do
        I.call $ I.sparse_la_solveWithGuess i s b x0 px
        x <- peek px
        SM.SparseMatrix <$> FC.newForeignPtr x (I.call $ I.sparse_free x)
-}

-- | Success if the iterations converged, and NoConvergence otherwise. 
info :: (MonadIO m, I.Elem a b) => SolverT a b m ComputationInfo
info = _get_prop I.sparse_la_info >>= \x -> return (toEnum x)

-- | The tolerance error reached during the last solve. It is a close approximation of the true relative residual error @|Ax-b|/|b|@.
error :: (MonadIO m, I.Elem a b) => SolverT a b m Double
error = _get_prop I.sparse_la_error

-- | The number of iterations performed during the last solve
iterations :: (MonadIO m, I.Elem a b) => SolverT a b m Int
iterations = _get_prop I.sparse_la_iterations

_get_prop :: (I.Cast c d, MonadIO m, Storable c) => (SolverInfo -> I.CSolverPtr a b -> Ptr c -> IO CString) -> SolverT a b m d
_get_prop f = ask >>= \(i,fs) -> liftIO $
    withForeignPtr fs $ \s -> alloca $ \px -> do
        I.call $ f i s px
        I.cast <$> peek px

_set_prop :: (I.Cast c d, MonadIO m, Storable c) => (SolverInfo -> I.CSolverPtr a b -> d -> IO CString) -> c -> SolverT a b m ()
_set_prop f x = ask >>= \(i,fs) -> liftIO $
    withForeignPtr fs $ \s -> I.call $ f i s (I.cast x)
