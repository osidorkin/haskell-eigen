import qualified Data.Eigen.Matrix as M
import Data.Eigen.SparseMatrix
import Data.Eigen.SparseLA as LA
import Control.Monad.Trans

main = do
    let
        a :: SparseMatrixXd
        b :: SparseMatrixXd
        a = fromDenseList [[1,2,3], [4,5,6], [7,8,10]]
        b = fromDenseList [[3],[3],[4]]
    putStrLn "Here is the matrix A:"
    print $ a

    putStrLn "Here is the vector b:"
    print $ b

    runSolverT (SparseLU COLAMDOrdering) $ do
        compute a
        x <- solve b
        info >>= lift.print
        determinant >>= lift.print
        lift $ putStrLn "The solution is:"
        lift $ print x
