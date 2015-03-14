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

foreign import ccall "eigen-proxy.h eigen_initParallel" c_initParallel :: IO ()
foreign import ccall "eigen-proxy.h eigen_setNbThreads" c_setNbThreads :: CInt -> IO ()

foreign import ccall "eigen-proxy.h eigen_add"         c_add         :: Ptr CDouble -> CInt -> CInt -> Ptr CDouble -> CInt -> CInt -> IO CString
foreign import ccall "eigen-proxy.h eigen_sub"         c_sub         :: Ptr CDouble -> CInt -> CInt -> Ptr CDouble -> CInt -> CInt -> IO CString
foreign import ccall "eigen-proxy.h eigen_mul"         c_mul         :: Ptr CDouble -> CInt -> CInt -> Ptr CDouble -> CInt -> CInt -> IO CString
foreign import ccall "eigen-proxy.h eigen_transpose"   c_transpose   :: Ptr CDouble -> CInt -> CInt -> Ptr CDouble -> CInt -> CInt -> IO CString
foreign import ccall "eigen-proxy.h eigen_inverse"     c_inverse     :: Ptr CDouble -> CInt -> CInt -> Ptr CDouble -> CInt -> CInt -> IO CString
foreign import ccall "eigen-proxy.h eigen_adjoint"     c_adjoint     :: Ptr CDouble -> CInt -> CInt -> Ptr CDouble -> CInt -> CInt -> IO CString
foreign import ccall "eigen-proxy.h eigen_conjugate"   c_conjugate   :: Ptr CDouble -> CInt -> CInt -> Ptr CDouble -> CInt -> CInt -> IO CString
foreign import ccall "eigen-proxy.h eigen_normalize"   c_normalize   :: Ptr CDouble -> CInt -> CInt -> IO CString
foreign import ccall "eigen-proxy.h eigen_norm"        c_norm        :: Ptr CDouble -> CInt -> CInt -> IO CDouble
foreign import ccall "eigen-proxy.h eigen_squaredNorm" c_squaredNorm :: Ptr CDouble -> CInt -> CInt -> IO CDouble
foreign import ccall "eigen-proxy.h eigen_blueNorm"    c_blueNorm    :: Ptr CDouble -> CInt -> CInt -> IO CDouble
foreign import ccall "eigen-proxy.h eigen_hypotNorm"   c_hypotNorm   :: Ptr CDouble -> CInt -> CInt -> IO CDouble
foreign import ccall "eigen-proxy.h eigen_determinant" c_determinant :: Ptr CDouble -> CInt -> CInt -> IO CDouble
