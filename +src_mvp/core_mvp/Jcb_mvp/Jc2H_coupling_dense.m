function [Hout] = Jc2H_coupling_dense(Jin, Zbc_Kop)
% function Hout = Jc2H_coupling_dense(Jin, Zbc_Kop)

Hout = Zbc_Kop * Jin;