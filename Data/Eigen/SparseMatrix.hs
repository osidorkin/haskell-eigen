{-# LANGUAGE GADTs #-}
{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE FlexibleContexts #-}

module Data.Eigen.SparseMatrix (
    -- * SparseMatrix type
    -- | SparseMatrix aliases follows Eigen naming convention
    SparseMatrix,
    SparseMatrixXf,
    SparseMatrixXd,
    SparseMatrixXcf,
    SparseMatrixXcd,
    -- * Accessing matrix data
    cols,
    rows,
    coeff,
    (!),
    -- * Matrix conversions
    fromList,
    toList,
    fromDenseList,
    toDenseList,
    fromMatrix,
    toMatrix,
    -- * Matrix properties
    norm,
    squaredNorm,
    blueNorm,
    block,
    --compress,
    --uncompress,
    --compressed,
    nonZeros,
    --innerSize,
    --outerSize,
    -- * Basic matrix algebra
    add,
    sub,
    mul,
    -- * Matrix transformations
    pruned,
    scale,
    transpose,
    adjoint,
    --lowerTriangle,
    --upperTriangle,
    -- * Matrix serialization
    encode,
    decode,
) where

import qualified Prelude as P
import Prelude hiding (map)
import qualified Data.List as L
import Data.Complex
import Data.IORef
import Foreign.C.Types
import Foreign.C.String
import Foreign.Storable
import Foreign.Ptr
import Foreign.ForeignPtr
import Foreign.Marshal.Alloc
import Control.Monad
import Control.Applicative
import qualified Data.Eigen.Matrix as M
import qualified Data.Eigen.Matrix.Mutable as MM
import qualified Foreign.Concurrent as FC
import qualified Data.Eigen.Internal as I
import qualified Data.Vector.Storable as VS
import qualified Data.Vector.Storable.Mutable as VSM
import qualified Data.ByteString as BS
import qualified Data.ByteString.Lazy as BSL
import qualified Data.ByteString.Internal as BSI

{-| A versatible sparse matrix representation.

This class implements a more versatile variants of the common compressed row/column storage format.
Each colmun's (resp. row) non zeros are stored as a pair of value with associated row (resp. colmiun) index.
All the non zeros are stored in a single large buffer.
Unlike the compressed format, there might be extra space inbetween the nonzeros of two successive colmuns
(resp. rows) such that insertion of new non-zero can be done with limited memory reallocation and copies.

A call to the function 'makeCompressed' turns the matrix into the standard compressed format compatible with many library.

Implementation deails of SparseMatrix are intentionally hidden behind ForeignPtr bacause Eigen doesn't provide mapping over plain data for sparse matricies.

For more infomration please see Eigen documentation page: <http://eigen.tuxfamily.org/dox/classEigen_1_1SparseMatrix.html>
-}

data SparseMatrix a b where
    SparseMatrix :: I.Elem a b => !(ForeignPtr (I.CSparseMatrix a b)) -> SparseMatrix a b

-- | Alias for single precision sparse matrix
type SparseMatrixXf = SparseMatrix Float CFloat
-- | Alias for double precision sparse matrix
type SparseMatrixXd = SparseMatrix Double CDouble
-- | Alias for single previsiom sparse matrix of complex numbers
type SparseMatrixXcf = SparseMatrix (Complex Float) (I.CComplex CFloat)
-- | Alias for double prevision sparse matrix of complex numbers
type SparseMatrixXcd = SparseMatrix (Complex Double) (I.CComplex CDouble)

-- | Pretty prints the sparse matrix
instance (I.Elem a b, Show a) => Show (SparseMatrix a b) where
    show m = concat [
        "SparseMatrix ", show (rows m), "x", show (cols m),
        "\n", L.intercalate "\n" $ P.map (L.intercalate "\t" . P.map show) $ toDenseList m, "\n"]

-- | Shortcuts for basic matrix math
instance I.Elem a b => Num (SparseMatrix a b) where
    (*) = mul
    (+) = add
    (-) = sub
    fromInteger x = fromList 1 1 [(0,0,fromInteger x)]
    signum = map signum
    abs = map abs
    negate = map negate

-- | Not exposed, For internal use donly
map :: I.Elem a b => (a -> a) -> SparseMatrix a b -> SparseMatrix a b
map f m = fromList (rows m) (cols m) . P.map (\(r,c,v) -> (r,c,f v)) . toList $ m

mk :: I.Elem a b => Ptr (I.CSparseMatrix a b) -> IO (SparseMatrix a b)
mk p = SparseMatrix <$> FC.newForeignPtr p (I.call $ I.sparse_free p)

-- | Number of rows for the sparse matrix
rows :: I.Elem a b => SparseMatrix a b -> Int
rows = _unop I.sparse_rows (return . I.cast)

-- | Number of columns for the sparse matrix
cols :: I.Elem a b => SparseMatrix a b -> Int
cols = _unop I.sparse_cols (return . I.cast)

-- | Matrix coefficient at given row and col
coeff :: I.Elem a b => Int -> Int -> SparseMatrix a b -> a
coeff row col (SparseMatrix fp) = I.performIO $ withForeignPtr fp $ \p -> alloca $ \pq -> do
    I.call $ I.sparse_coeff p (I.cast row) (I.cast col) pq
    I.cast <$> peek pq

-- | Matrix coefficient at given row and col
(!) :: I.Elem a b => SparseMatrix a b -> (Int, Int) -> a
(!) m (row, col) = coeff row col m

{-| For vectors, the l2 norm, and for matrices the Frobenius norm.
    In both cases, it consists in the square root of the sum of the square of all the matrix entries.
    For vectors, this is also equals to the square root of the dot product of this with itself.
-}
norm :: I.Elem a b => SparseMatrix a b -> a
norm = _unop I.sparse_norm (return . I.cast)

-- | For vectors, the squared l2 norm, and for matrices the Frobenius norm. In both cases, it consists in the sum of the square of all the matrix entries. For vectors, this is also equals to the dot product of this with itself.
squaredNorm :: I.Elem a b => SparseMatrix a b -> a
squaredNorm = _unop I.sparse_squaredNorm (return . I.cast)

-- | The l2 norm of the matrix using the Blue's algorithm. A Portable Fortran Program to Find the Euclidean Norm of a Vector, ACM TOMS, Vol 4, Issue 1, 1978.
blueNorm :: I.Elem a b => SparseMatrix a b -> a
blueNorm = _unop I.sparse_blueNorm (return . I.cast)

-- | Extract rectangular block from sparse matrix defined by startRow startCol blockRows blockCols
block :: I.Elem a b => Int -> Int -> Int -> Int -> SparseMatrix a b -> SparseMatrix a b
block row col rows cols = _unop (\p pq -> I.sparse_block p (I.cast row) (I.cast col) (I.cast rows) (I.cast cols) pq) mk

-- | not exposed currently
compress :: I.Elem a b => SparseMatrix a b -> SparseMatrix a b
compress = _unop I.sparse_compress mk

-- | not exposed currently
uncompress :: I.Elem a b => SparseMatrix a b -> SparseMatrix a b
uncompress = _unop I.sparse_uncompress mk

-- | not exposed currently
compressed :: I.Elem a b => SparseMatrix a b -> Bool
compressed = _unop I.sparse_isCompressed (return . (/=0))

-- | Number of non-zeros elements in the sparse matrix
nonZeros :: I.Elem a b => SparseMatrix a b -> Int
nonZeros = _unop I.sparse_nonZeros (return . I.cast)

-- | not exposed currently
innerSize :: I.Elem a b => SparseMatrix a b -> Int
innerSize = _unop I.sparse_innerSize (return . I.cast)

-- | not exposed currently
outerSize :: I.Elem a b => SparseMatrix a b -> Int
outerSize = _unop I.sparse_outerSize (return . I.cast)

-- | Suppresses all nonzeros which are much smaller than reference under the tolerence epsilon
pruned :: I.Elem a b => a -> SparseMatrix a b -> SparseMatrix a b
pruned r = _unop (\p pq -> alloca $ \pr -> poke pr (I.cast r) >> I.sparse_prunedRef p pr pq) mk

-- | Multiply matrix on a given scalar
scale :: I.Elem a b => a -> SparseMatrix a b -> SparseMatrix a b
scale x = _unop (\p pq -> alloca $ \px -> poke px (I.cast x) >> I.sparse_scale p px pq) mk

-- | Transpose of the sparse matrix
transpose :: I.Elem a b => SparseMatrix a b -> SparseMatrix a b
transpose = _unop I.sparse_transpose mk

-- | Adjoint of the sparse matrix
adjoint :: I.Elem a b => SparseMatrix a b -> SparseMatrix a b
adjoint = _unop I.sparse_adjoint mk

-- | not exposed currently
lowerTriangle :: I.Elem a b => SparseMatrix a b -> SparseMatrix a b
lowerTriangle = _unop I.sparse_lowerTriangle mk

-- | not exposed currently
upperTriangle :: I.Elem a b => SparseMatrix a b -> SparseMatrix a b
upperTriangle = _unop I.sparse_upperTriangle mk

-- | Adding two sparse matrices by adding the corresponding entries together. You can use @(+)@ function as well.
add :: I.Elem a b => SparseMatrix a b -> SparseMatrix a b -> SparseMatrix a b
add = _binop I.sparse_add mk

-- | Subtracting two sparse matrices by subtracting the corresponding entries together. You can use @(-)@ function as well.
sub :: I.Elem a b => SparseMatrix a b -> SparseMatrix a b -> SparseMatrix a b
sub = _binop I.sparse_sub mk

-- | Matrix multiplication. You can use @(*)@ function as well.
mul :: I.Elem a b => SparseMatrix a b -> SparseMatrix a b -> SparseMatrix a b
mul = _binop I.sparse_mul mk

-- | Construct sparse matrix of given size from the list of triplets (row, col, val)
fromList :: I.Elem a b => Int -> Int -> [(Int, Int, a)] -> SparseMatrix a b
fromList rows cols xs = I.performIO $ VS.unsafeWith vs $ \p -> alloca $ \pq -> do
    I.call $ I.sparse_fromList (I.cast rows) (I.cast cols) p (I.cast $ VS.length vs) pq
    peek pq >>= mk
    where vs = VS.fromList $ P.map (\(row,col,val) -> I.CTriplet (I.cast row) (I.cast col) (I.cast val)) xs

-- | Convert sparse matrix to the list of triplets (row, col, val). Compressed elements will not be included
toList :: I.Elem a b => SparseMatrix a b -> [(Int, Int, a)]
toList m@(SparseMatrix fp) = I.performIO $ do
    let s = nonZeros m
    xs <- VSM.new s
    withForeignPtr fp $ \p ->
        VSM.unsafeWith xs $ \q ->
            I.call $ I.sparse_toList p q (I.cast s)
    let f (I.CTriplet row col val) = (I.cast row, I.cast col, I.cast val)
    P.map f . VS.toList <$> VS.unsafeFreeze xs

-- | Construct sparse matrix of two-dimensional list of values. Matrix dimensions will be detected automatically. Zero values will be compressed.
fromDenseList :: (I.Elem a b, Eq a) => [[a]] -> SparseMatrix a b
fromDenseList list = fromList rows cols $ do
    (row, vals) <- zip [0..] list
    (col, val) <- zip [0..] vals
    guard $ val /= 0
    return (row, col, val)
    where
        rows = length list
        cols = L.foldl' max 0 $ P.map length list

-- | Convert sparse matrix to (rows X cols) dense list of values
toDenseList :: I.Elem a b => SparseMatrix a b -> [[a]]
toDenseList m = [[coeff row col m | col <- [0 .. cols m - 1]] | row <- [0 .. rows m - 1]]

-- | Construct sparse matrix from dense matrix. Zero elements will be compressed
fromMatrix :: I.Elem a b => M.Matrix a b -> SparseMatrix a b
fromMatrix m1 = I.performIO $ alloca $ \pm0 ->
    M.unsafeWith m1 $ \vals rows cols -> do
        I.call $ I.sparse_fromMatrix vals rows cols pm0
        peek pm0 >>= mk

-- | Construct dense matrix from sparse matrix
toMatrix :: I.Elem a b => SparseMatrix a b -> M.Matrix a b
toMatrix m1@(SparseMatrix fp) = I.performIO $ do
    m0 <- MM.new (rows m1) (cols m1)
    MM.unsafeWith m0 $ \vals rows cols ->
        withForeignPtr fp $ \pm1 ->
            I.call $ I.sparse_toMatrix pm1 vals rows cols
    M.unsafeFreeze m0

-- | Encode the sparse matrix as a lazy byte string
encode :: I.Elem a b => SparseMatrix a b -> BSL.ByteString
encode m@(SparseMatrix fp) = I.performIO $ do
    let size = nonZeros m
    tris <- VSM.new size
    withForeignPtr fp $ \p ->
        VSM.unsafeWith tris $ \q ->
            I.call $ I.sparse_toList p q (I.cast size)
    tris <- VS.unsafeFreeze tris
    let 
        tri@(I.CTriplet _ _ val) = VS.head tris

    return $ BSL.fromChunks [
        encodeInt (I.code val),
        encodeInt (I.cast $ rows m),
        encodeInt (I.cast $ cols m),
        encodeInt (I.cast $ size),
        let (fp, fs) = VS.unsafeToForeignPtr0 tris in BSI.PS (castForeignPtr fp) 0 (fs * sizeOf tri)]

    where
        encodeInt :: CInt -> BS.ByteString
        encodeInt x = BSI.unsafeCreate (sizeOf x) $ (`poke` x) . castPtr
        

-- | Decode sparse matrix from the lazy byte string
decode :: forall a b . I.Elem a b => BSL.ByteString -> SparseMatrix a b
decode st = I.performIO $ do
    ref <- newIORef st
    let next size = readIORef ref >>= \a ->
            let (b,c) = BSL.splitAt (fromIntegral size) a
            in if BSL.length b /= fromIntegral size
                then error "SparseMatrix.decode: stream exhausted"
                else do
                    writeIORef ref c
                    return . BS.concat . BSL.toChunks $ b

        val = undefined :: b
        tri = undefined :: I.CTriplet b
        triSize = sizeOf tri

    code <- I.decodeInt <$> next I.intSize
    when (code /= I.code val) $
        fail "SparseMatrix.decode: wrong matrix type"
    
    rows <- I.decodeInt <$> next I.intSize
    cols <- I.decodeInt <$> next I.intSize
    size <- I.decodeInt <$> next I.intSize
    
    BSI.PS fp fo _ <- next (I.cast size * triSize)
    BSL.null <$> readIORef ref >>= (`unless` fail "SparseMatrix.decode: stream underrun")

    let tris = VS.unsafeFromForeignPtr0 (I.plusForeignPtr fp fo) (I.cast size)
    
    VS.unsafeWith tris $ \p -> alloca $ \pq -> do
        I.call $ I.sparse_fromList rows cols p size pq
        peek pq >>= mk

_unop :: Storable c => (I.CSparseMatrixPtr a b -> Ptr c -> IO CString) -> (c -> IO d) -> SparseMatrix a b -> d
_unop f g (SparseMatrix fp) = I.performIO $
    withForeignPtr fp $ \p ->
        alloca $ \pq -> do
            I.call (f p pq)
            peek pq >>= g

_binop :: Storable c => (I.CSparseMatrixPtr a b -> I.CSparseMatrixPtr a b -> Ptr c -> IO CString) -> (c -> IO d) -> SparseMatrix a b -> SparseMatrix a b -> d
_binop f g (SparseMatrix fp1) (SparseMatrix fp2) = I.performIO $
    withForeignPtr fp1 $ \p1 ->
        withForeignPtr fp2 $ \p2 ->
            alloca $ \pq -> do
                I.call (f p1 p2 pq)
                peek pq >>= g
