{- |
Some Eigen's algorithms can exploit the multiple cores present in your hardware. To this end, it is enough to enable OpenMP on your compiler, for instance: GCC: -fopenmp ICC: -openmp MSVC: check the respective option in the build properties. You can control the number of thread that will be used using 'setNbThreads'
-}

module Data.Eigen.Parallel where

import Data.Eigen.Internal

-- | Must be call first when calling Eigen from multiple threads
initParallel :: IO ()
initParallel = c_initParallel

-- | Sets the max number of threads reserved for Eigen
setNbThreads :: Int -> IO ()
setNbThreads = c_setNbThreads . cast
