{-# LANGUAGE CPP #-}
{-# LANGUAGE FlexibleInstances #-}
{-# LANGUAGE ForeignFunctionInterface #-}
{-# LANGUAGE FunctionalDependencies #-}
{-# LANGUAGE MultiParamTypeClasses #-}
{-# LANGUAGE RecordWildCards #-}
{-# LANGUAGE GeneralizedNewtypeDeriving #-}

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

-}

module Data.Eigen.SparseLA (
    -- * Sparse Solvers
    Solver,
    DirectSolver,
    IterativeSolver,
    OrderingMethod(..),
    Preconditioner(..),
    ConjugateGradient(ConjugateGradient),
    BiCGSTAB(BiCGSTAB),
    SparseLU(SparseLU),
    SparseQR(SparseQR),
    ComputationInfo(..),
    SolverT,
    runSolverT,
    -- * The Compute step
    {- |
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
    -}
    analyzePattern,
    factorize,
    compute,
    -- * The Solve step
    {- |
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
    -}
    solve,
    --solveWithGuess,
    info,
    -- * Iterative Solvers
    tolerance,
    setTolerance,
    maxIterations,
    setMaxIterations,
    Data.Eigen.SparseLA.error,
    iterations,
    -- * SparseQR Solver
    matrixR,
    matrixQ,
    rank,
    setPivotThreshold,
    -- * SparseLU Solver
    setSymmetric,
    matrixL,
    matrixU,
    determinant,
    absDeterminant,
    signDeterminant,
    logAbsDeterminant,
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
#if __GLASGOW_HASKELL__ >= 710
#else
import Control.Applicative
#endif
import qualified Data.Eigen.Internal as I
import qualified Data.Eigen.SparseMatrix as SM

{- | Ordering methods for sparse matrices. They are typically used to reduce the number of elements during the sparse matrix
    decomposition (@LLT@, @LU@, @QR@). Precisely, in a preprocessing step, a permutation matrix @P@ is computed using those ordering methods
    and applied to the columns of the matrix. Using for instance the sparse Cholesky decomposition, it is expected that the nonzeros
    elements in @LLT(A*P)@ will be much smaller than that in @LLT(A)@.
-}
data OrderingMethod
    -- | The column approximate minimum degree ordering The matrix should be in column-major and compressed format
    = COLAMDOrdering
    -- | The natural ordering (identity)
    | NaturalOrdering deriving (Show, Read)

data Preconditioner
    {- | A preconditioner based on the digonal entries

        It allows to approximately solve for A.x = b problems assuming A is a diagonal matrix.
        In other words, this preconditioner neglects all off diagonal entries and, in Eigen's language, solves for:
        @
        A.diagonal().asDiagonal() . x = b
        @
        This preconditioner is suitable for both selfadjoint and general problems.
        The diagonal entries are pre-inverted and stored into a dense vector.

        A variant that has yet to be implemented would attempt to preserve the norm of each column.
    -}
    = DiagonalPreconditioner
    -- | A naive preconditioner which approximates any matrix as the identity matrix
    | IdentityPreconditioner deriving (Show, Read)


class I.Code s => Solver s where
-- | For direct methods, the solution is computed at the machine precision.
class Solver s => DirectSolver s where
-- | Sometimes, the solution need not be too accurate.
-- In this case, the iterative methods are more suitable and the desired accuracy can be set before the solve step using `setTolerance`.
class Solver s => IterativeSolver s where

{- | A conjugate gradient solver for sparse self-adjoint problems.

    This class allows to solve for @A.x = b@ sparse linear problems using a conjugate gradient algorithm. The sparse matrix @A@ must be selfadjoint.

    The maximal number of iterations and tolerance value can be controlled via the `setMaxIterations` and `setTolerance` methods.
    The defaults are the size of the problem for the maximal number of iterations and @epsilon@ for the tolerance
-}
data ConjugateGradient = ConjugateGradient Preconditioner deriving (Show, Read)
instance Solver ConjugateGradient
instance IterativeSolver ConjugateGradient
instance I.Code ConjugateGradient where
    code (ConjugateGradient DiagonalPreconditioner) = 0
    code (ConjugateGradient IdentityPreconditioner) = 1

{- | A bi conjugate gradient stabilized solver for sparse square problems.

    This class allows to solve for @A.x = b@ sparse linear problems using a bi conjugate gradient stabilized algorithm.
    The vectors @x@ and @b@ can be either dense or sparse.

    The maximal number of iterations and tolerance value can be controlled via the `setMaxIterations` and `setTolerance` methods.
    The defaults are the size of the problem for the maximal number of iterations and @epsilon@ for the tolerance
-}
data BiCGSTAB = BiCGSTAB Preconditioner deriving (Show, Read)
instance Solver BiCGSTAB
instance IterativeSolver BiCGSTAB
instance I.Code BiCGSTAB where
    code (BiCGSTAB DiagonalPreconditioner) = 2
    code (BiCGSTAB IdentityPreconditioner) = 3

{- | Sparse supernodal LU factorization for general matrices.

    This class implements the supernodal LU factorization for general matrices. It uses the main techniques from the sequential
    <http://crd-legacy.lbl.gov/~xiaoye/SuperLU/ SuperLU package>. It handles transparently real and complex arithmetics with
    single and double precision, depending on the scalar type of your input matrix. The code has been optimized to provide BLAS-3
    operations during supernode-panel updates. It benefits directly from the built-in high-performant Eigen BLAS routines.
    Moreover, when the size of a supernode is very small, the BLAS calls are avoided to enable a better optimization from the compiler.
    For best performance, you should compile it with NDEBUG flag to avoid the numerous bounds checking on vectors.

    An important parameter of this class is the ordering method. It is used to reorder the columns
    (and eventually the rows) of the matrix to reduce the number of new elements that are created during
    numerical factorization. The cheapest method available is COLAMD.
    See <http://eigen.tuxfamily.org/dox/group__OrderingMethods__Module.html OrderingMethods module> for the list of
    built-in and external ordering methods.
-}
data SparseLU = SparseLU OrderingMethod deriving (Show, Read)
instance Solver SparseLU
instance DirectSolver SparseLU
instance I.Code SparseLU where
    code (SparseLU NaturalOrdering) = 4
    code (SparseLU COLAMDOrdering) = 5

{- | Sparse left-looking rank-revealing QR factorization.

    This class implements a left-looking rank-revealing QR decomposition of sparse matrices. When a column has a norm less than a given
    tolerance it is implicitly permuted to the end. The QR factorization thus obtained is given by @A*P = Q*R@ where @R@ is upper triangular or trapezoidal.

    @P@ is the column permutation which is the product of the fill-reducing and the rank-revealing permutations.

    @Q@ is the orthogonal matrix represented as products of Householder reflectors.

    @R@ is the sparse triangular or trapezoidal matrix. The later occurs when @A@ is rank-deficient.
-}
data SparseQR = SparseQR OrderingMethod deriving (Show, Read)
instance Solver SparseQR
instance DirectSolver SparseQR
instance I.Code SparseQR where
    code (SparseQR NaturalOrdering) = 6
    code (SparseQR COLAMDOrdering) = 7


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

type SolverT s a b m = ReaderT (s, ForeignPtr (I.CSolver a b)) m

runSolverT :: (Solver s, MonadIO m, I.Elem a b) => s -> SolverT s a b m c -> m c
runSolverT i f = do
    fs <- liftIO $ alloca $ \ps -> do
        I.call $ I.sparse_la_newSolver i ps
        s <- peek ps
        FC.newForeignPtr s (I.call $ I.sparse_la_freeSolver i s)
    runReaderT f (i,fs)

-- | Initializes the iterative solver for the sparsity pattern of the matrix @A@ for further solving @Ax=b@ problems.
analyzePattern :: (Solver s, MonadIO m, I.Elem a b) => SM.SparseMatrix a b -> SolverT s a b m ()
analyzePattern (SM.SparseMatrix fa) = ask >>= \(i,fs) -> liftIO $
    withForeignPtr fs $ \s ->
    withForeignPtr fa $ \a ->
        I.call $ I.sparse_la_analyzePattern i s a

-- | Initializes the iterative solver with the numerical values of the matrix @A@ for further solving @Ax=b@ problems.
factorize :: (Solver s, MonadIO m, I.Elem a b) => SM.SparseMatrix a b -> SolverT s a b m ()
factorize (SM.SparseMatrix fa) = ask >>= \(i,fs) -> liftIO $
    withForeignPtr fs $ \s ->
    withForeignPtr fa $ \a ->
        I.call $ I.sparse_la_factorize i s a

-- | Initializes the iterative solver with the matrix @A@ for further solving @Ax=b@ problems.
--
-- The `compute` method is equivalent to calling both `analyzePattern` and `factorize`.
compute :: (Solver s, MonadIO m, I.Elem a b) => SM.SparseMatrix a b -> SolverT s a b m ()
compute (SM.SparseMatrix fa) = ask >>= \(i,fs) -> liftIO $
    withForeignPtr fs $ \s ->
    withForeignPtr fa $ \a ->
        I.call $ I.sparse_la_compute i s a

-- | An expression of the solution @x@ of @Ax=b@ using the current decomposition of @A@.
solve :: (Solver s, MonadIO m, I.Elem a b) => SM.SparseMatrix a b -> SolverT s a b m (SM.SparseMatrix a b)
solve (SM.SparseMatrix fb) = ask >>= \(i,fs) -> liftIO $
    withForeignPtr fs $ \s ->
    withForeignPtr fb $ \b ->
    alloca $ \px -> do
        I.call $ I.sparse_la_solve i s b px
        x <- peek px
        SM.SparseMatrix <$> FC.newForeignPtr x (I.call $ I.sparse_free x)

{-
-- | The solution @x@ of @Ax=b@ using the current decomposition of @A@ and @x0@ as an initial solution.
solveWithGuess :: (MonadIO m, I.Elem a b) => SM.SparseMatrix a b -> SM.SparseMatrix a b -> SolverT s a b m (SM.SparseMatrix a b)
solveWithGuess (SM.SparseMatrix fb) (SM.SparseMatrix fx0) = ask >>= \(i,fs) -> liftIO $
    withForeignPtr fs $ \s ->
    withForeignPtr fb $ \b ->
    withForeignPtr fx0 $ \x0 ->
    alloca $ \px -> do
        I.call $ I.sparse_la_solveWithGuess i s b x0 px
        x <- peek px
        SM.SparseMatrix <$> FC.newForeignPtr x (I.call $ I.sparse_free x)
-}

-- |
-- * `Success` if the iterations converged or computation was succesful
-- * `NumericalIssue` if the factorization reports a numerical problem
-- * `NoConvergence` if the iterations are not converged
-- * `InvalidInput` if the input matrix is invalid
info :: (Solver s, MonadIO m, I.Elem a b) => SolverT s a b m ComputationInfo
info = _get_prop I.sparse_la_info >>= \x -> return (toEnum x)

-- | The tolerance threshold used by the stopping criteria.
tolerance :: (IterativeSolver s, MonadIO m, I.Elem a b) => SolverT s a b m Double
tolerance = _get_prop I.sparse_la_tolerance

-- | Sets the tolerance threshold used by the stopping criteria.
--
--   This value is used as an upper bound to the relative residual error: @|Ax-b|/|b|@. The default value is the machine precision given by @epsilon@
setTolerance :: (IterativeSolver s, MonadIO m, I.Elem a b) => Double -> SolverT s a b m ()
setTolerance = _set_prop I.sparse_la_setTolerance

-- | The max number of iterations. It is either the value setted by setMaxIterations or, by default, twice the number of columns of the matrix.
maxIterations :: (IterativeSolver s, MonadIO m, I.Elem a b) => SolverT s a b m Int
maxIterations = _get_prop I.sparse_la_maxIterations

-- | Sets the max number of iterations. Default is twice the number of columns of the matrix.
setMaxIterations :: (IterativeSolver s, MonadIO m, I.Elem a b) => Int -> SolverT s a b m ()
setMaxIterations = _set_prop I.sparse_la_setMaxIterations

-- | The tolerance error reached during the last solve. It is a close approximation of the true relative residual error @|Ax-b|/|b|@.
error :: (IterativeSolver s, MonadIO m, I.Elem a b) => SolverT s a b m Double
error = _get_prop I.sparse_la_error

-- | The number of iterations performed during the last solve
iterations :: (IterativeSolver s, MonadIO m, I.Elem a b) => SolverT s a b m Int
iterations = _get_prop I.sparse_la_iterations

-- | Returns the @b@ sparse upper triangular matrix @R@ of the QR factorization.
matrixR :: (MonadIO m, I.Elem a b) => SolverT SparseQR a b m (SM.SparseMatrix a b)
matrixR = _get_matrix I.sparse_la_matrixR

-- | Returns the matrix @Q@ as products of sparse Householder reflectors.
matrixQ :: (MonadIO m, I.Elem a b) => SolverT SparseQR a b m (SM.SparseMatrix a b)
matrixQ = _get_matrix I.sparse_la_matrixQ

-- | Sets the threshold that is used to determine linearly dependent columns during the factorization.
--
-- In practice, if during the factorization the norm of the column that has to be eliminated is below
-- this threshold, then the entire column is treated as zero, and it is moved at the end.
setPivotThreshold :: (MonadIO m, I.Elem a b) => Double -> SolverT SparseQR a b m ()
setPivotThreshold = _set_prop I.sparse_la_setPivotThreshold

-- | Returns the number of non linearly dependent columns as determined by the pivoting threshold.
rank :: (MonadIO m, I.Elem a b) => SolverT SparseQR a b m Int
rank = _get_prop I.sparse_la_rank

-- | Indicate that the pattern of the input matrix is symmetric
setSymmetric :: (MonadIO m, I.Elem a b) => Bool -> SolverT SparseLU a b m ()
setSymmetric = _set_prop I.sparse_la_setSymmetric . fromEnum

-- | Returns the matrix @L@
matrixL :: (MonadIO m, I.Elem a b) => SolverT SparseLU a b m (SM.SparseMatrix a b)
matrixL = _get_matrix I.sparse_la_matrixL

-- | Returns the matrix @U@
matrixU :: (MonadIO m, I.Elem a b) => SolverT SparseLU a b m (SM.SparseMatrix a b)
matrixU = _get_matrix I.sparse_la_matrixU

-- | The determinant of the matrix.
determinant :: (MonadIO m, I.Elem a b) => SolverT SparseLU a b m a
determinant = _get_prop I.sparse_la_determinant

-- | The natural log of the absolute value of the determinant of the matrix of which this is the QR decomposition
--
-- This method is useful to work around the risk of overflow/underflow that's inherent to the determinant computation.
logAbsDeterminant :: (MonadIO m, I.Elem a b) => SolverT SparseLU a b m a
logAbsDeterminant = _get_prop I.sparse_la_logAbsDeterminant

-- | The absolute value of the determinant of the matrix of which *this is the QR decomposition.
--
-- A determinant can be very big or small, so for matrices of large enough dimension, there is a risk of overflow/underflow.
-- One way to work around that is to use `logAbsDeterminant` instead.
absDeterminant :: (MonadIO m, I.Elem a b) => SolverT SparseLU a b m a
absDeterminant = _get_prop I.sparse_la_absDeterminant

-- | A number representing the sign of the determinant
signDeterminant :: (MonadIO m, I.Elem a b) => SolverT SparseLU a b m a
signDeterminant = _get_prop I.sparse_la_signDeterminant

_get_prop :: (I.Cast c d, Solver s, MonadIO m, Storable c) => (s -> I.CSolverPtr a b -> Ptr c -> IO CString) -> SolverT s a b m d
_get_prop f = ask >>= \(i,fs) -> liftIO $
    withForeignPtr fs $ \s -> alloca $ \px -> do
        I.call $ f i s px
        I.cast <$> peek px

_get_matrix :: (Solver s, MonadIO m, I.Elem a b) => (s -> I.CSolverPtr a b -> Ptr (I.CSparseMatrixPtr a b) -> IO CString) -> SolverT s a b m (SM.SparseMatrix a b)
_get_matrix f = ask >>= \(i,fs) -> liftIO $
    withForeignPtr fs $ \s -> alloca $ \px -> do
        I.call $ f i s px
        x <- peek px
        SM.SparseMatrix <$> FC.newForeignPtr x (I.call $ I.sparse_free x)

_set_prop :: (I.Cast c d, Solver s, MonadIO m, Storable c) => (s -> I.CSolverPtr a b -> d -> IO CString) -> c -> SolverT s a b m ()
_set_prop f x = ask >>= \(i,fs) -> liftIO $
    withForeignPtr fs $ \s -> I.call $ f i s (I.cast x)
