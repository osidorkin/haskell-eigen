{-# LANGUAGE MultiParamTypeClasses, ForeignFunctionInterface, EmptyDataDecls #-}
module Data.Eigen.Internal where

import Foreign.Ptr
import Foreign.C.Types
import Foreign.C.String
import Control.Monad
import System.IO.Unsafe

class Cast a b where
    cast :: a -> b

instance Cast CDouble Double where
    cast (CDouble x) = x

instance Cast Double CDouble where
    cast = CDouble

instance Cast CInt Int where
    cast = fromIntegral

instance Cast Int CInt where
    cast = fromIntegral

performIO :: IO a -> a
performIO = unsafeDupablePerformIO

foreign import ccall "eigen-proxy.h free" c_freeString :: CString -> IO ()

call :: IO CString -> IO ()
call func = func >>= \c_str -> when (c_str /= nullPtr) $
    peekCString c_str >>= \str -> c_freeString c_str >> fail str

foreign import ccall "eigen-proxy.h eigen_setNbThreads" c_setNbThreads :: CInt -> IO ()
foreign import ccall "eigen-proxy.h eigen_getNbThreads" c_getNbThreads :: IO CInt

foreign import ccall "eigen-proxy.h eigen_random"      c_random      :: Ptr CDouble -> CInt -> CInt -> IO CString
foreign import ccall "eigen-proxy.h eigen_add"         c_add         :: Ptr CDouble -> CInt -> CInt -> Ptr CDouble -> CInt -> CInt -> Ptr CDouble -> CInt -> CInt -> IO CString
foreign import ccall "eigen-proxy.h eigen_sub"         c_sub         :: Ptr CDouble -> CInt -> CInt -> Ptr CDouble -> CInt -> CInt -> Ptr CDouble -> CInt -> CInt -> IO CString
foreign import ccall "eigen-proxy.h eigen_mul"         c_mul         :: Ptr CDouble -> CInt -> CInt -> Ptr CDouble -> CInt -> CInt -> Ptr CDouble -> CInt -> CInt -> IO CString
foreign import ccall "eigen-proxy.h eigen_diagonal"    c_diagonal    :: Ptr CDouble -> CInt -> CInt -> Ptr CDouble -> CInt -> CInt -> IO CString
foreign import ccall "eigen-proxy.h eigen_transpose"   c_transpose   :: Ptr CDouble -> CInt -> CInt -> Ptr CDouble -> CInt -> CInt -> IO CString
foreign import ccall "eigen-proxy.h eigen_inverse"     c_inverse     :: Ptr CDouble -> CInt -> CInt -> Ptr CDouble -> CInt -> CInt -> IO CString
foreign import ccall "eigen-proxy.h eigen_adjoint"     c_adjoint     :: Ptr CDouble -> CInt -> CInt -> Ptr CDouble -> CInt -> CInt -> IO CString
foreign import ccall "eigen-proxy.h eigen_conjugate"   c_conjugate   :: Ptr CDouble -> CInt -> CInt -> Ptr CDouble -> CInt -> CInt -> IO CString
foreign import ccall "eigen-proxy.h eigen_normalize"   c_normalize   :: Ptr CDouble -> CInt -> CInt -> IO CString
foreign import ccall "eigen-proxy.h eigen_sum"         c_sum         :: Ptr CDouble -> CInt -> CInt -> IO CDouble
foreign import ccall "eigen-proxy.h eigen_prod"        c_prod        :: Ptr CDouble -> CInt -> CInt -> IO CDouble
foreign import ccall "eigen-proxy.h eigen_mean"        c_mean        :: Ptr CDouble -> CInt -> CInt -> IO CDouble
foreign import ccall "eigen-proxy.h eigen_norm"        c_norm        :: Ptr CDouble -> CInt -> CInt -> IO CDouble
foreign import ccall "eigen-proxy.h eigen_trace"       c_trace       :: Ptr CDouble -> CInt -> CInt -> IO CDouble
foreign import ccall "eigen-proxy.h eigen_squaredNorm" c_squaredNorm :: Ptr CDouble -> CInt -> CInt -> IO CDouble
foreign import ccall "eigen-proxy.h eigen_blueNorm"    c_blueNorm    :: Ptr CDouble -> CInt -> CInt -> IO CDouble
foreign import ccall "eigen-proxy.h eigen_hypotNorm"   c_hypotNorm   :: Ptr CDouble -> CInt -> CInt -> IO CDouble
foreign import ccall "eigen-proxy.h eigen_determinant" c_determinant :: Ptr CDouble -> CInt -> CInt -> IO CDouble

foreign import ccall "eigen-proxy.h eigen_rank"         c_rank       :: CInt -> Ptr CInt -> Ptr CDouble -> CInt -> CInt -> IO CString
foreign import ccall "eigen-proxy.h eigen_image"        c_image      :: CInt -> Ptr (Ptr CDouble) -> Ptr CInt -> Ptr CInt -> Ptr CDouble -> CInt -> CInt -> IO CString
foreign import ccall "eigen-proxy.h eigen_kernel"       c_kernel     :: CInt -> Ptr (Ptr CDouble) -> Ptr CInt -> Ptr CInt -> Ptr CDouble -> CInt -> CInt -> IO CString
foreign import ccall "eigen-proxy.h free"               c_free       :: Ptr a -> IO ()
foreign import ccall "eigen-proxy.h eigen_solve"        c_solve       :: CInt -> Ptr CDouble -> CInt -> CInt -> Ptr CDouble -> CInt -> CInt -> Ptr CDouble -> CInt -> CInt -> IO CString
foreign import ccall "eigen-proxy.h eigen_relativeError" c_relativeError :: Ptr CDouble -> Ptr CDouble -> CInt -> CInt -> Ptr CDouble -> CInt -> CInt -> Ptr CDouble -> CInt -> CInt -> IO CString
