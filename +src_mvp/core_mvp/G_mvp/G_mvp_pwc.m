function Jout = G_mvp_pwc(Jin, dx, dom_dims)
% function G_mvp_pwc(Jin, r, dims) implements Gram matrix vector product 
% Inputs:
%        - Jin: a vector of volumetrix polarizarion currents
%        - dx:  resolution (i.e. size of voxel)

% reshape input vector to tensor
Jin_tensor = reshape(Jin, dom_dims);

% Gram scaling factor for PWC is voxel volume
G = dx^3;

% scale currents by Gram scaling
Jout = G * Jin_tensor;
     
