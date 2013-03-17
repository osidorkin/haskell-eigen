{-# LANGUAGE MultiParamTypeClasses, ForeignFunctionInterface, EmptyDataDecls #-}
module Data.Eigen.Internal where

import Foreign.Ptr
import Foreign.C.Types
import Foreign.C.String
import Control.Monad

data C_MatrixXd

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

foreign import ccall "eigen-proxy.h free" c_freeString :: CString -> IO ()

call :: IO CString -> IO ()
call func = func >>= \c_str -> when (c_str /= nullPtr) $
    peekCString c_str >>= \str -> c_freeString c_str >> error str

foreign import ccall "eigen-proxy.h eigen_initParallel" c_initParallel :: IO ()
foreign import ccall "eigen-proxy.h eigen_setNbThreads" c_setNbThreads :: CInt -> IO ()

foreign import ccall "eigen-proxy.h eigen_create"      c_create :: CInt -> CInt -> IO (Ptr C_MatrixXd)
foreign import ccall "eigen-proxy.h &eigen_destroy"    c_destroy :: FunPtr (Ptr C_MatrixXd -> IO ())
foreign import ccall "eigen-proxy.h eigen_clone"       c_clone :: Ptr C_MatrixXd -> IO (Ptr C_MatrixXd)
foreign import ccall "eigen-proxy.h eigen_get"         c_get :: Ptr CDouble -> Ptr C_MatrixXd -> CInt -> CInt -> IO CString
foreign import ccall "eigen-proxy.h eigen_set"         c_set :: Ptr C_MatrixXd -> CInt -> CInt -> CDouble -> IO CString
foreign import ccall "eigen-proxy.h eigen_data"        c_data :: Ptr C_MatrixXd -> IO (Ptr CDouble)
foreign import ccall "eigen-proxy.h eigen_rows"        c_rows :: Ptr C_MatrixXd -> IO CInt
foreign import ccall "eigen-proxy.h eigen_cols"        c_cols :: Ptr C_MatrixXd -> IO CInt
foreign import ccall "eigen-proxy.h eigen_copy"        c_copy :: Ptr C_MatrixXd -> Ptr C_MatrixXd -> IO ()
foreign import ccall "eigen-proxy.h eigen_resize"      c_resize :: Ptr C_MatrixXd -> CInt -> CInt -> IO ()
foreign import ccall "eigen-proxy.h eigen_add"         c_add :: Ptr C_MatrixXd -> Ptr C_MatrixXd -> Ptr C_MatrixXd -> IO CString
foreign import ccall "eigen-proxy.h eigen_sub"         c_sub :: Ptr C_MatrixXd -> Ptr C_MatrixXd -> Ptr C_MatrixXd -> IO CString
foreign import ccall "eigen-proxy.h eigen_mul"         c_mul :: Ptr C_MatrixXd -> Ptr C_MatrixXd -> Ptr C_MatrixXd -> IO CString
foreign import ccall "eigen-proxy.h eigen_transpose"   c_transpose :: Ptr C_MatrixXd -> IO CString
foreign import ccall "eigen-proxy.h eigen_inverse"     c_inverse :: Ptr C_MatrixXd -> IO CString
foreign import ccall "eigen-proxy.h eigen_adjoint"     c_adjoint :: Ptr C_MatrixXd -> IO CString
foreign import ccall "eigen-proxy.h eigen_normalize"   c_normalize :: Ptr C_MatrixXd -> IO CString
foreign import ccall "eigen-proxy.h eigen_norm"        c_norm :: Ptr C_MatrixXd -> IO CDouble
foreign import ccall "eigen-proxy.h eigen_squaredNorm" c_squaredNorm :: Ptr C_MatrixXd -> IO CDouble
foreign import ccall "eigen-proxy.h eigen_blueNorm"    c_blueNorm :: Ptr C_MatrixXd -> IO CDouble
foreign import ccall "eigen-proxy.h eigen_hypotNorm"   c_hypotNorm :: Ptr C_MatrixXd -> IO CDouble
foreign import ccall "eigen-proxy.h eigen_determinant" c_determinant :: Ptr C_MatrixXd -> IO CDouble
