{-# OPTIONS_HADDOCK hide #-}
{-# OPTIONS_GHC -fno-warn-orphans #-}

{-# LANGUAGE CPP #-} 
{-# LANGUAGE EmptyDataDecls  #-}
{-# LANGUAGE FlexibleInstances  #-}
{-# LANGUAGE ForeignFunctionInterface  #-}
{-# LANGUAGE FunctionalDependencies  #-}
{-# LANGUAGE MultiParamTypeClasses #-}
{-# LANGUAGE ScopedTypeVariables  #-}

module Data.Eigen.Internal where

import Foreign.Ptr
import Foreign.ForeignPtr
import Foreign.Storable
import Foreign.C.Types
import Foreign.C.String
import Control.Monad
#if __GLASGOW_HASKELL__ >= 710
#else
import Control.Applicative
#endif
import System.IO.Unsafe
import Data.Complex
import Data.Bits
import Data.Binary
import Data.Binary.Get
import Data.Binary.Put
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

instance Storable a => Binary (VS.Vector a) where
    put vs = put (BS.length bs) >> putByteString bs where
        (fp,fs) = VS.unsafeToForeignPtr0 vs
        es = sizeOf (VS.head vs)
        bs = BSI.fromForeignPtr (castForeignPtr fp) 0 (fs * es)
        
    get = get >>= getByteString >>= \bs -> let
        (fp,fo,fs) = BSI.toForeignPtr bs
        es = sizeOf (VS.head vs)
        vs = VS.unsafeFromForeignPtr0 (Data.Eigen.Internal.plusForeignPtr fp fo) (fs `div` es)
        in return vs

-- | Complex number for FFI with the same memory layout as std::complex\<T\>
data CComplex a = CComplex !a !a deriving Show

instance Storable a => Storable (CComplex a) where
    sizeOf _ = sizeOf (undefined :: a) * 2
    alignment _ = alignment (undefined :: a)
    poke p (CComplex x y) = do
        pokeElemOff (castPtr p) 0 x
        pokeElemOff (castPtr p) 1 y
    peek p = CComplex
        <$> peekElemOff (castPtr p) 0
        <*> peekElemOff (castPtr p) 1

data CTriplet a = CTriplet !CInt !CInt !a deriving Show

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

instance Cast a b => Cast (CTriplet a) (Int, Int, b) where; cast (CTriplet x y z) = (cast x, cast y, cast z)
instance Cast a b => Cast (Int, Int, a) (CTriplet b) where; cast (x,y,z) = CTriplet (cast x) (cast y) (cast z)

intSize :: Int
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

data CSolver a b
type CSolverPtr a b = Ptr (CSolver a b)

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

newtype MagicCode = MagicCode CInt deriving Eq

instance Binary MagicCode where
    put (MagicCode code) = putWord32be $ fromIntegral code
    get = MagicCode . fromIntegral <$> getWord32be

magicCode :: Code a => a -> MagicCode
magicCode x = MagicCode (code x `xor` 0x45696730)

#let api1 name, args = "foreign import ccall \"eigen_%s\" c_%s :: CInt -> %s\n%s :: forall b . Code b => %s\n%s = c_%s (code (undefined :: b))", #name, #name, args, #name, args, #name, #name

#api1 random,        "Ptr b -> CInt -> CInt -> IO CString"
#api1 identity,      "Ptr b -> CInt -> CInt -> IO CString"
#api1 add,           "Ptr b -> CInt -> CInt -> Ptr b -> CInt -> CInt -> Ptr b -> CInt -> CInt -> IO CString"
#api1 sub,           "Ptr b -> CInt -> CInt -> Ptr b -> CInt -> CInt -> Ptr b -> CInt -> CInt -> IO CString"
#api1 mul,           "Ptr b -> CInt -> CInt -> Ptr b -> CInt -> CInt -> Ptr b -> CInt -> CInt -> IO CString"
#api1 diagonal,      "Ptr b -> CInt -> CInt -> Ptr b -> CInt -> CInt -> IO CString"
#api1 transpose,     "Ptr b -> CInt -> CInt -> Ptr b -> CInt -> CInt -> IO CString"
#api1 inverse,       "Ptr b -> CInt -> CInt -> Ptr b -> CInt -> CInt -> IO CString"
#api1 adjoint,       "Ptr b -> CInt -> CInt -> Ptr b -> CInt -> CInt -> IO CString"
#api1 conjugate,     "Ptr b -> CInt -> CInt -> Ptr b -> CInt -> CInt -> IO CString"
#api1 normalize,     "Ptr b -> CInt -> CInt -> IO CString"
#api1 sum,           "Ptr b -> Ptr b -> CInt -> CInt -> IO CString"
#api1 prod,          "Ptr b -> Ptr b -> CInt -> CInt -> IO CString"
#api1 mean,          "Ptr b -> Ptr b -> CInt -> CInt -> IO CString"
#api1 norm,          "Ptr b -> Ptr b -> CInt -> CInt -> IO CString"
#api1 trace,         "Ptr b -> Ptr b -> CInt -> CInt -> IO CString"
#api1 squaredNorm,   "Ptr b -> Ptr b -> CInt -> CInt -> IO CString"
#api1 blueNorm,      "Ptr b -> Ptr b -> CInt -> CInt -> IO CString"
#api1 hypotNorm,     "Ptr b -> Ptr b -> CInt -> CInt -> IO CString"
#api1 determinant,   "Ptr b -> Ptr b -> CInt -> CInt -> IO CString"
#api1 rank,          "CInt -> Ptr CInt -> Ptr b -> CInt -> CInt -> IO CString"
#api1 image,         "CInt -> Ptr (Ptr b) -> Ptr CInt -> Ptr CInt -> Ptr b -> CInt -> CInt -> IO CString"
#api1 kernel,        "CInt -> Ptr (Ptr b) -> Ptr CInt -> Ptr CInt -> Ptr b -> CInt -> CInt -> IO CString"
#api1 solve,         "CInt -> Ptr b -> CInt -> CInt -> Ptr b -> CInt -> CInt -> Ptr b -> CInt -> CInt -> IO CString"
#api1 relativeError, "Ptr b -> Ptr b -> CInt -> CInt -> Ptr b -> CInt -> CInt -> Ptr b -> CInt -> CInt -> IO CString"

#let api2 name, args = "foreign import ccall \"eigen_%s\" c_%s :: CInt -> %s\n%s :: forall a b . Code b => %s\n%s = c_%s (code (undefined :: b))", #name, #name, args, #name, args, #name, #name

#api2 sparse_new,           "CInt -> CInt -> Ptr (CSparseMatrixPtr a b) -> IO CString"
#api2 sparse_clone,         "CSparseMatrixPtr a b -> Ptr (CSparseMatrixPtr a b) -> IO CString"
#api2 sparse_fromList,      "CInt -> CInt -> Ptr (CTriplet b) -> CInt -> Ptr (CSparseMatrixPtr a b) -> IO CString"
#api2 sparse_toList,        "CSparseMatrixPtr a b -> Ptr (CTriplet b) -> CInt -> IO CString"
#api2 sparse_free,          "CSparseMatrixPtr a b -> IO CString"
#api2 sparse_makeCompressed,"CSparseMatrixPtr a b -> Ptr (CSparseMatrixPtr a b) -> IO CString"
#api2 sparse_uncompress,    "CSparseMatrixPtr a b -> Ptr (CSparseMatrixPtr a b) -> IO CString"
#api2 sparse_isCompressed,  "CSparseMatrixPtr a b -> Ptr CInt -> IO CString"
#api2 sparse_transpose,     "CSparseMatrixPtr a b -> Ptr (CSparseMatrixPtr a b) -> IO CString"
#api2 sparse_adjoint,       "CSparseMatrixPtr a b -> Ptr (CSparseMatrixPtr a b) -> IO CString"
#api2 sparse_pruned,        "CSparseMatrixPtr a b -> Ptr (CSparseMatrixPtr a b) -> IO CString"
#api2 sparse_prunedRef,     "CSparseMatrixPtr a b -> Ptr b -> Ptr (CSparseMatrixPtr a b) -> IO CString"
#api2 sparse_scale,         "CSparseMatrixPtr a b -> Ptr b -> Ptr (CSparseMatrixPtr a b) -> IO CString"
#api2 sparse_nonZeros,      "CSparseMatrixPtr a b -> Ptr CInt -> IO CString"
#api2 sparse_innerSize,     "CSparseMatrixPtr a b -> Ptr CInt -> IO CString"
#api2 sparse_outerSize,     "CSparseMatrixPtr a b -> Ptr CInt -> IO CString"
#api2 sparse_coeff,         "CSparseMatrixPtr a b -> CInt -> CInt -> Ptr b -> IO CString"
#api2 sparse_coeffRef,      "CSparseMatrixPtr a b -> CInt -> CInt -> Ptr (Ptr b) -> IO CString"
#api2 sparse_cols,          "CSparseMatrixPtr a b -> Ptr CInt -> IO CString"
#api2 sparse_rows,          "CSparseMatrixPtr a b -> Ptr CInt -> IO CString"
#api2 sparse_norm,          "CSparseMatrixPtr a b -> Ptr b -> IO CString"
#api2 sparse_squaredNorm,   "CSparseMatrixPtr a b -> Ptr b -> IO CString"
#api2 sparse_blueNorm,      "CSparseMatrixPtr a b -> Ptr b -> IO CString"
#api2 sparse_add,           "CSparseMatrixPtr a b -> CSparseMatrixPtr a b -> Ptr (CSparseMatrixPtr a b) -> IO CString"
#api2 sparse_sub,           "CSparseMatrixPtr a b -> CSparseMatrixPtr a b -> Ptr (CSparseMatrixPtr a b) -> IO CString"
#api2 sparse_mul,           "CSparseMatrixPtr a b -> CSparseMatrixPtr a b -> Ptr (CSparseMatrixPtr a b) -> IO CString"
#api2 sparse_block,         "CSparseMatrixPtr a b -> CInt -> CInt -> CInt -> CInt -> Ptr (CSparseMatrixPtr a b) -> IO CString"
#api2 sparse_fromMatrix,    "Ptr b -> CInt -> CInt -> Ptr (CSparseMatrixPtr a b) -> IO CString"
#api2 sparse_toMatrix,      "CSparseMatrixPtr a b -> Ptr b -> CInt -> CInt -> IO CString"
#api2 sparse_values,        "CSparseMatrixPtr a b -> Ptr CInt -> Ptr (Ptr b) -> IO CString"
#api2 sparse_outerStarts,   "CSparseMatrixPtr a b -> Ptr CInt -> Ptr (Ptr CInt) -> IO CString"
#api2 sparse_innerIndices,  "CSparseMatrixPtr a b -> Ptr CInt -> Ptr (Ptr CInt) -> IO CString"
#api2 sparse_innerNNZs,     "CSparseMatrixPtr a b -> Ptr CInt -> Ptr (Ptr CInt) -> IO CString"
#api2 sparse_setZero,       "CSparseMatrixPtr a b -> IO CString"
#api2 sparse_setIdentity,   "CSparseMatrixPtr a b -> IO CString"
#api2 sparse_reserve,       "CSparseMatrixPtr a b -> CInt -> IO CString"
#api2 sparse_resize,        "CSparseMatrixPtr a b -> CInt -> CInt -> IO CString"

#api2 sparse_conservativeResize,    "CSparseMatrixPtr a b -> CInt -> CInt -> IO CString"
#api2 sparse_compressInplace,       "CSparseMatrixPtr a b -> IO CString"
#api2 sparse_uncompressInplace,     "CSparseMatrixPtr a b -> IO CString"


#let api3 name, args = "foreign import ccall \"eigen_%s\" c_%s :: CInt -> CInt -> %s\n%s :: forall s a b . (Code s, Code b) => s -> %s\n%s s = c_%s (code (undefined :: b)) (code s)", #name, #name, args, #name, args, #name, #name

#api3 sparse_la_newSolver,          "Ptr (CSolverPtr a b) -> IO CString"
#api3 sparse_la_freeSolver,         "CSolverPtr a b -> IO CString"
#api3 sparse_la_factorize,          "CSolverPtr a b -> CSparseMatrixPtr a b -> IO CString"
#api3 sparse_la_analyzePattern,     "CSolverPtr a b -> CSparseMatrixPtr a b -> IO CString"
#api3 sparse_la_compute,            "CSolverPtr a b -> CSparseMatrixPtr a b -> IO CString"
#api3 sparse_la_tolerance,          "CSolverPtr a b -> Ptr CDouble -> IO CString"
#api3 sparse_la_setTolerance,       "CSolverPtr a b -> CDouble -> IO CString"
#api3 sparse_la_maxIterations,      "CSolverPtr a b -> Ptr CInt -> IO CString"
#api3 sparse_la_setMaxIterations,   "CSolverPtr a b -> CInt -> IO CString"
#api3 sparse_la_info,               "CSolverPtr a b -> Ptr CInt -> IO CString"
#api3 sparse_la_error,              "CSolverPtr a b -> Ptr CDouble -> IO CString"
#api3 sparse_la_iterations,         "CSolverPtr a b -> Ptr CInt -> IO CString"
#api3 sparse_la_solve,              "CSolverPtr a b -> CSparseMatrixPtr a b -> Ptr (CSparseMatrixPtr a b) -> IO CString"
-- #api3 sparse_la_solveWithGuess,     "CSolverPtr a b -> CSparseMatrixPtr a b -> CSparseMatrixPtr a b -> Ptr (CSparseMatrixPtr a b) -> IO CString"
#api3 sparse_la_matrixQ,            "CSolverPtr a b -> Ptr (CSparseMatrixPtr a b) -> IO CString"
#api3 sparse_la_matrixR,            "CSolverPtr a b -> Ptr (CSparseMatrixPtr a b) -> IO CString"
#api3 sparse_la_setPivotThreshold,  "CSolverPtr a b -> CDouble -> IO CString"
#api3 sparse_la_rank,               "CSolverPtr a b -> Ptr CInt -> IO CString"
#api3 sparse_la_matrixL,            "CSolverPtr a b -> Ptr (CSparseMatrixPtr a b) -> IO CString"
#api3 sparse_la_matrixU,            "CSolverPtr a b -> Ptr (CSparseMatrixPtr a b) -> IO CString"
#api3 sparse_la_setSymmetric,       "CSolverPtr a b -> CInt -> IO CString"
#api3 sparse_la_determinant,        "CSolverPtr a b -> Ptr b -> IO CString"
#api3 sparse_la_logAbsDeterminant,  "CSolverPtr a b -> Ptr b -> IO CString"
#api3 sparse_la_absDeterminant,     "CSolverPtr a b -> Ptr b -> IO CString"
#api3 sparse_la_signDeterminant,    "CSolverPtr a b -> Ptr b -> IO CString"
