function Jout = G_mvp_precond(Jin, dx, dims)
% function set_G_mvp_precond(Jin, dx, dims)
% used to assemble preconditioner for the coupled solver
% unlike G_mvp_pwx() it accepts a vector of body tissue mat. properties and
% scales it by a grammian; 

% ----------------------------------------------------------------------- %
% Input parameters:                                                       %
% ----------------------------------------------------------------------- %
%   - Jin  - a vector of material properties of scatterer [Nb_vox, 1]
%   - dx   - voxel resolution 
%   - dims - object, that keeps track of dimensions in MARIE 2.0
% ----------------------------------------------------------------------- %
% Output parameters:                                                      %
% ----------------------------------------------------------------------- %
%   - Jout  - a vector of material scaled by gram matrix [ql * Nb_vox, 1]
% ----------------------------------------------------------------------- %

% set output dimensions
Jtemp = repmat(Jin, dims.l, 1);

% voxel volume
V = dx.^3;

if 1 == dims.l
    Jout = repmat(Jtemp .* V, dims.q, 1);
else
    Gram = V .* [ones(size(Jin)); 
                 ones(size(Jin)) ./ 12; 
                 ones(size(Jin)) ./ 12;
                 ones(size(Jin)) ./ 12];
                           
    Jout = repmat(Gram .* Jtemp, dims.q,1);
end