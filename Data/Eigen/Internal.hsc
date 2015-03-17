{-# LANGUAGE MultiParamTypeClasses, ForeignFunctionInterface, ScopedTypeVariables, FunctionalDependencies, FlexibleInstances, EmptyDataDecls, CPP #-}
module Data.Eigen.Internal where

import Foreign.Ptr
import Foreign.ForeignPtr
import Foreign.Storable
import Foreign.C.Types
import Foreign.C.String
import Control.Monad
import Control.Applicative
import System.IO.Unsafe
import Data.Complex
import qualified Data.Vector.Storable as VS
import qualified Data.ByteString as BS
import qualified Data.ByteString.Internal as BSI

class (Num a, Cast a b, Cast b a, Storable b, Code b) => Elem a b | a -> b where

instance Elem Float CFloat where
instance Elem Double CDouble where
instance Elem (Complex Float) (CComplex CFloat) where
instance Elem (Complex Double) (CComplex CDouble) where

class Cast a b where
	cast :: a -> b

-- | Complex number for FFI with the same memory layout as std::complex\<T\>
data CComplex a = CComplex !a !a

instance Storable a => Storable (CComplex a) where
	sizeOf _ = sizeOf (undefined :: a) * 2
	alignment _ = alignment (undefined :: a)
	poke p (CComplex x y) = do
		pokeElemOff (castPtr p) 0 x
		pokeElemOff (castPtr p) 1 y
	peek p = CComplex
		<$> peekElemOff (castPtr p) 0
		<*> peekElemOff (castPtr p) 1

data CTriplet a = CTriplet !CInt !CInt !a

instance Storable a => Storable (CTriplet a) where
	sizeOf _ = sizeOf (undefined :: a) + sizeOf (undefined :: CInt) * 2
	alignment _ = alignment (undefined :: CInt)
	poke p (CTriplet row col val) = do
		pokeElemOff (castPtr p) 0 row
		pokeElemOff (castPtr p) 1 col
		pokeByteOff p (sizeOf (undefined :: CInt) * 2) val
	peek p = CTriplet
		<$> peekElemOff (castPtr p) 0
		<*> peekElemOff (castPtr p) 1
		<*> peekByteOff p (sizeOf (undefined :: CInt) * 2)

instance Cast CInt Int where; cast = fromIntegral
instance Cast Int CInt where; cast = fromIntegral
instance Cast CFloat Float where; cast (CFloat x) = x
instance Cast Float CFloat where; cast = CFloat
instance Cast CDouble Double where; cast (CDouble x) = x
instance Cast Double CDouble where; cast = CDouble
instance Cast (CComplex CFloat) (Complex Float) where; cast (CComplex x y) = cast x :+ cast y
instance Cast (Complex Float) (CComplex CFloat) where; cast (x :+ y) = CComplex (cast x) (cast y)
instance Cast (CComplex CDouble) (Complex Double) where; cast (CComplex x y) = cast x :+ cast y
instance Cast (Complex Double) (CComplex CDouble) where; cast (x :+ y) = CComplex (cast x) (cast y)

intSize = sizeOf (undefined :: CInt)

encodeInt :: CInt -> BS.ByteString
encodeInt x = BSI.unsafeCreate (sizeOf x) $ (`poke` x) . castPtr

decodeInt :: BS.ByteString -> CInt
decodeInt (BSI.PS fp fo fs)
	| fs == sizeOf x = x
	| otherwise = error "decodeInt: wrong buffer size"
	where x = performIO $ withForeignPtr fp $ peek . (`plusPtr` fo)

data CSparseMatrix a b
type CSparseMatrixPtr a b = Ptr (CSparseMatrix a b)

performIO :: IO a -> a
performIO = unsafeDupablePerformIO

plusForeignPtr :: ForeignPtr a -> Int -> ForeignPtr b
plusForeignPtr fp fo = castForeignPtr fp' where
	vs :: VS.Vector CChar
	vs = VS.unsafeFromForeignPtr (castForeignPtr fp) fo 0
	(fp', _) = VS.unsafeToForeignPtr0 vs

foreign import ccall "eigen-proxy.h free" c_freeString :: CString -> IO ()

call :: IO CString -> IO ()
call func = func >>= \c_str -> when (c_str /= nullPtr) $
	peekCString c_str >>= \str -> c_freeString c_str >> fail str

foreign import ccall "eigen-proxy.h free" free :: Ptr a -> IO ()

foreign import ccall "eigen-proxy.h eigen_setNbThreads" c_setNbThreads :: CInt -> IO ()
foreign import ccall "eigen-proxy.h eigen_getNbThreads" c_getNbThreads :: IO CInt

class Code a where; code :: a -> CInt
instance Code CFloat where; code _ = 0
instance Code CDouble where; code _ = 1
instance Code (CComplex CFloat) where; code _ = 2
instance Code (CComplex CDouble) where; code _ = 3

#let api name, args = "foreign import ccall \"eigen_%s\" c_%s :: CInt -> %s\n%s :: forall b . Code b => %s\n%s = c_%s (code (undefined :: b))", #name, #name, args, #name, args, #name, #name

#api random, 		"Ptr b -> CInt -> CInt -> IO CString"
#api identity, 		"Ptr b -> CInt -> CInt -> IO CString"
#api add, 			"Ptr b -> CInt -> CInt -> Ptr b -> CInt -> CInt -> Ptr b -> CInt -> CInt -> IO CString"
#api sub, 			"Ptr b -> CInt -> CInt -> Ptr b -> CInt -> CInt -> Ptr b -> CInt -> CInt -> IO CString"
#api mul, 			"Ptr b -> CInt -> CInt -> Ptr b -> CInt -> CInt -> Ptr b -> CInt -> CInt -> IO CString"
#api diagonal, 		"Ptr b -> CInt -> CInt -> Ptr b -> CInt -> CInt -> IO CString"
#api transpose, 	"Ptr b -> CInt -> CInt -> Ptr b -> CInt -> CInt -> IO CString"
#api inverse, 		"Ptr b -> CInt -> CInt -> Ptr b -> CInt -> CInt -> IO CString"
#api adjoint, 		"Ptr b -> CInt -> CInt -> Ptr b -> CInt -> CInt -> IO CString"
#api conjugate, 	"Ptr b -> CInt -> CInt -> Ptr b -> CInt -> CInt -> IO CString"
#api normalize, 	"Ptr b -> CInt -> CInt -> IO CString"
#api sum, 			"Ptr b -> Ptr b -> CInt -> CInt -> IO CString"
#api prod, 			"Ptr b -> Ptr b -> CInt -> CInt -> IO CString"
#api mean, 			"Ptr b -> Ptr b -> CInt -> CInt -> IO CString"
#api norm, 			"Ptr b -> Ptr b -> CInt -> CInt -> IO CString"
#api trace, 		"Ptr b -> Ptr b -> CInt -> CInt -> IO CString"
#api squaredNorm, 	"Ptr b -> Ptr b -> CInt -> CInt -> IO CString"
#api blueNorm, 		"Ptr b -> Ptr b -> CInt -> CInt -> IO CString"
#api hypotNorm, 	"Ptr b -> Ptr b -> CInt -> CInt -> IO CString"
#api determinant, 	"Ptr b -> Ptr b -> CInt -> CInt -> IO CString"
#api rank, 			"CInt -> Ptr CInt -> Ptr b -> CInt -> CInt -> IO CString"
#api image, 		"CInt -> Ptr (Ptr b) -> Ptr CInt -> Ptr CInt -> Ptr b -> CInt -> CInt -> IO CString"
#api kernel, 		"CInt -> Ptr (Ptr b) -> Ptr CInt -> Ptr CInt -> Ptr b -> CInt -> CInt -> IO CString"
#api solve, 		"CInt -> Ptr b -> CInt -> CInt -> Ptr b -> CInt -> CInt -> Ptr b -> CInt -> CInt -> IO CString"
#api relativeError, "Ptr b -> Ptr b -> CInt -> CInt -> Ptr b -> CInt -> CInt -> Ptr b -> CInt -> CInt -> IO CString"

#let api2 name, args = "foreign import ccall \"eigen_%s\" c_%s :: CInt -> %s\n%s :: forall a b . Code b => %s\n%s = c_%s (code (undefined :: b))", #name, #name, args, #name, args, #name, #name

#api2 sparse_fromList, 		"CInt -> CInt -> Ptr (CTriplet b) -> CInt -> Ptr (CSparseMatrixPtr a b) -> IO CString"
#api2 sparse_toList, 		"CSparseMatrixPtr a b -> Ptr (CTriplet b) -> CInt -> IO CString"
#api2 sparse_free, 			"CSparseMatrixPtr a b -> IO CString"
#api2 sparse_compress, 		"CSparseMatrixPtr a b -> Ptr (CSparseMatrixPtr a b) -> IO CString"
#api2 sparse_uncompress, 	"CSparseMatrixPtr a b -> Ptr (CSparseMatrixPtr a b) -> IO CString"
#api2 sparse_isCompressed, 	"CSparseMatrixPtr a b -> Ptr CInt -> IO CString"
#api2 sparse_transpose, 	"CSparseMatrixPtr a b -> Ptr (CSparseMatrixPtr a b) -> IO CString"
#api2 sparse_adjoint, 		"CSparseMatrixPtr a b -> Ptr (CSparseMatrixPtr a b) -> IO CString"
#api2 sparse_pruned, 		"CSparseMatrixPtr a b -> Ptr (CSparseMatrixPtr a b) -> IO CString"
#api2 sparse_prunedRef, 	"CSparseMatrixPtr a b -> Ptr b -> Ptr (CSparseMatrixPtr a b) -> IO CString"
#api2 sparse_scale, 		"CSparseMatrixPtr a b -> Ptr b -> Ptr (CSparseMatrixPtr a b) -> IO CString"
#api2 sparse_lowerTriangle, "CSparseMatrixPtr a b -> Ptr (CSparseMatrixPtr a b) -> IO CString"
#api2 sparse_upperTriangle, "CSparseMatrixPtr a b -> Ptr (CSparseMatrixPtr a b) -> IO CString"
#api2 sparse_nonZeros, 		"CSparseMatrixPtr a b -> Ptr CInt -> IO CString"
#api2 sparse_innerSize, 	"CSparseMatrixPtr a b -> Ptr CInt -> IO CString"
#api2 sparse_outerSize, 	"CSparseMatrixPtr a b -> Ptr CInt -> IO CString"
#api2 sparse_coeff, 		"CSparseMatrixPtr a b -> CInt -> CInt -> Ptr b -> IO CString"
#api2 sparse_cols, 			"CSparseMatrixPtr a b -> Ptr CInt -> IO CString"
#api2 sparse_rows, 			"CSparseMatrixPtr a b -> Ptr CInt -> IO CString"
#api2 sparse_norm, 			"CSparseMatrixPtr a b -> Ptr b -> IO CString"
#api2 sparse_squaredNorm, 	"CSparseMatrixPtr a b -> Ptr b -> IO CString"
#api2 sparse_blueNorm, 		"CSparseMatrixPtr a b -> Ptr b -> IO CString"
#api2 sparse_add, 			"CSparseMatrixPtr a b -> CSparseMatrixPtr a b -> Ptr (CSparseMatrixPtr a b) -> IO CString"
#api2 sparse_sub, 			"CSparseMatrixPtr a b -> CSparseMatrixPtr a b -> Ptr (CSparseMatrixPtr a b) -> IO CString"
#api2 sparse_mul, 			"CSparseMatrixPtr a b -> CSparseMatrixPtr a b -> Ptr (CSparseMatrixPtr a b) -> IO CString"
#api2 sparse_block,			"CSparseMatrixPtr a b -> CInt -> CInt -> CInt -> CInt -> Ptr (CSparseMatrixPtr a b) -> IO CString"
#api2 sparse_fromMatrix,	"Ptr b -> CInt -> CInt -> Ptr (CSparseMatrixPtr a b) -> IO CString"
#api2 sparse_toMatrix,		"CSparseMatrixPtr a b -> Ptr b -> CInt -> CInt -> IO CString"


