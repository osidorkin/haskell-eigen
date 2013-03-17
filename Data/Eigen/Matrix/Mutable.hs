{-# LANGUAGE RecordWildCards, MultiParamTypeClasses #-}
module Data.Eigen.Matrix.Mutable (
    MMatrix(..),
    new,
    rows,
    cols,
    copy,
    clone,
    resize,
    get,
    set,
    setRow,
    setCol,
    add,
    sub,
    mul,
    inverse,
    adjoint,
    transpose,
    normalize,
    norm,
    blueNorm,
    hypotNorm,
    squaredNorm,
    determinant,
    with
) where

import Foreign.Ptr
import Foreign.ForeignPtr
import Foreign.C.String
import Foreign.Storable
import Foreign.Marshal.Alloc
import Control.Monad
import Data.Eigen.Internal

data MMatrix = MMatrix { mm_fp :: ForeignPtr C_MatrixXd };

new :: Int -> Int -> IO MMatrix
new rows cols = do
    pm <- c_create (cast rows) (cast cols)
    fp <- newForeignPtr c_destroy pm
    return $ MMatrix fp

rows :: MMatrix -> IO Int
rows mm = fmap cast $ with mm c_rows

cols :: MMatrix -> IO Int
cols mm = fmap cast $ with mm c_cols

copy :: MMatrix -> MMatrix -> IO ()
copy dst src =
    with dst $ \pdst ->
    with src $ \psrc ->
        c_copy pdst psrc

clone :: MMatrix -> IO MMatrix
clone mm = with mm $ \src -> do
    dst <- c_clone src
    fp <- newForeignPtr c_destroy dst
    return $ MMatrix fp

resize :: MMatrix -> Int -> Int -> IO ()
resize mm rows cols = with mm $ \pm ->
    c_resize pm (cast rows) (cast cols)

get :: MMatrix -> Int -> Int -> IO Double
get mm row col =
    with mm $ \pm ->
    alloca $ \pv -> do
        call $ c_get pv pm (cast row) (cast col)
        fmap cast $ peek pv

set :: MMatrix -> Int -> Int -> Double -> IO ()
set mm row col val = with mm $ \pm -> call $ c_set pm (cast row) (cast col) (cast val)

setRow :: MMatrix -> Int -> [Double] -> IO ()
setRow mm row vals = zipWithM_ (\col val -> set mm row col val) [0..] vals

setCol :: MMatrix -> Int -> [Double] -> IO ()
setCol mm col vals = zipWithM_ (\row val -> set mm row col val) [0..] vals

binop :: (Ptr C_MatrixXd -> Ptr C_MatrixXd -> Ptr C_MatrixXd -> IO CString) -> MMatrix -> MMatrix -> MMatrix -> IO ()
binop f ret lhs rhs  =
    with ret $ \pret ->
    with lhs $ \plhs ->
    with rhs $ \prhs ->
        call $ f pret plhs prhs

add :: MMatrix -> MMatrix -> MMatrix -> IO ()
add = binop c_add

sub :: MMatrix -> MMatrix -> MMatrix -> IO ()
sub = binop c_sub

mul :: MMatrix -> MMatrix -> MMatrix -> IO ()
mul = binop c_mul

inpop :: (Ptr C_MatrixXd -> IO CString) -> MMatrix -> IO ()
inpop f mm = with mm $ \pm -> call $ f pm

inverse :: MMatrix -> IO ()
inverse = inpop c_inverse

adjoint :: MMatrix -> IO ()
adjoint = inpop c_adjoint

transpose :: MMatrix -> IO ()
transpose = inpop c_transpose

normalize :: MMatrix -> IO ()
normalize = inpop c_normalize

norm :: MMatrix -> IO Double
norm mm = fmap cast $ with mm c_norm

squaredNorm :: MMatrix -> IO Double
squaredNorm mm = fmap cast $ with mm c_squaredNorm

--stableNorm :: MMatrix -> IO Double
--stableNorm mm = fmap cast $ with mm c_stableNorm

blueNorm :: MMatrix -> IO Double
blueNorm mm = fmap cast $ with mm c_blueNorm

hypotNorm :: MMatrix -> IO Double
hypotNorm mm = fmap cast $ with mm c_hypotNorm

determinant :: MMatrix -> IO Double
determinant mm = fmap cast $ with mm c_determinant

with :: MMatrix -> (Ptr C_MatrixXd -> IO a) -> IO a
with (MMatrix fp) f = withForeignPtr fp f
