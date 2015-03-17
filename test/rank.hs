import Data.Eigen.Matrix
import Data.Eigen.LA

main = do
    let a = fromList [[1,2,5],[2,1,4],[3,0,3]] :: MatrixXf
    putStrLn "Here is the matrix A:"
    print a

    putStrLn "The rank of A is:"
    print $ rank FullPivLU a

    putStrLn "Here is a matrix whose columns form a basis of the null-space of A:"
    print $ kernel FullPivLU a

    putStrLn "Here is a matrix whose columns form a basis of the column-space of A:"
    print $ image FullPivLU a
