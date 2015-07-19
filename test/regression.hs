{-# LANGUAGE RecordWildCards #-}
import Data.Eigen.Matrix as M
import Data.Eigen.LA
import Data.List as L
import Control.Monad

main = do
    let
        a :: MatrixXd
        a = fromList [
            [1,3.02, 6.89],
            [1,2.01, 5.39],
            [1,2.41, 6.01],
            [1,2.09, 5.55],
            [1,2.58, 6.32]]

        b = fromList $ Prelude.map return [-4.32,-3.79,-4.01,-3.86,-4.10]

    print a
    print b
    forM_ [FullPivLU, HouseholderQR, ColPivHouseholderQR, FullPivHouseholderQR, JacobiSVD] $ \d -> do
        let x = solve d a b
            e = relativeError x a b
            e' = norm (a*x - b) / norm b
        putStrLn $ replicate 20 '*'
        print d
        print x
        print e
        print e'
    putStrLn "\n-2.34666 - 0.25349 x1 - 0.174965 x2"
    putStrLn "done"

    print $ (identity 4 4 :: MatrixXd)
    print $ M.normalize a
    print $ M.transpose a

    let
        a :: MatrixXd
        a = M.fromList [[0.68,  0.597,  -0.33],[-0.211,  0.823,  0.536],[ 0.566, -0.605, -0.444]]
        b = M.inverse a
    print a
    print b
    print $ a * b
    print $ linearRegression [
            [-4.32, 3.02, 6.89],
            [-3.79, 2.01, 5.39],
            [-4.01, 2.41, 6.01],
            [-3.86, 2.09, 5.55],
            [-4.10, 2.58, 6.32]]


