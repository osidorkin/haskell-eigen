import Data.Eigen.Matrix
import Data.Eigen.LA

main = do
    let
        a :: MatrixXd
        a = fromList [[1,2,3], [4,5,6], [7,8,10]]
        b = fromList [[3],[3],[4]]
        x = solve ColPivHouseholderQR a b
    putStrLn "Here is the matrix A:"
    print a

    putStrLn "Here is the vector b:"
    print b

    putStrLn "The solution is:"
    print x
