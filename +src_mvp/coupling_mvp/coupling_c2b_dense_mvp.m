function Vout = coupling_c2b_dense_mvp(Jin, Zbc)
% function coupling_c2b_dense_mvp(Jin, Zbc)
% provides an interface for computing fields, incident on the tissue 
% INPUT:
%       - Jin - surface coil currents
%       - Zbc - coupling operator
% OUTPUT:
%       - Vout - fields, incident on the tissue
% Author: Georgy Guryev, Cambridge, MA, 2019    


Vout = Zbc * Jin;
