function [Xext_bb, Xc_bb, Xb_bb] = get_BB_limits(coil, scatterer, projector, dims)

% get body voxel domain
idxS = scatterer.index_vie.S_1d;

% copy RHBM coordinated to vectors
x  = scatterer.dom_vie.x_tensor;
y  = scatterer.dom_vie.y_tensor;
z  = scatterer.dom_vie.z_tensor;

N_exp_1D = dims.N_exp_1D;
res     = scatterer.dom_vie.res;

% copy SCOIL coordinates to vectors
xc = coil.node(1,:);
yc = coil.node(2,:);
zc = coil.node(3,:);

% get near zone frame width 
% nbw = projector.near_boundary_width;

% define min distance from the coil to the boundary of the ext. domain
N_coil2bound = (N_exp_1D / 2) - 1;
N_exp_half   = (N_exp_1D - 1) / 2;

d_coil2bound = res * N_coil2bound;

%% get body domain + 2 voxels

Xb_bb = zeros(3,2);

Xb_bb(1,1) = min(x(idxS))-res; Xb_bb(1,2) = max(x(idxS))+res;
Xb_bb(2,1) = min(y(idxS))-res; Xb_bb(2,2) = max(y(idxS))+res;
Xb_bb(3,1) = min(z(idxS))-res; Xb_bb(3,2) = max(z(idxS))+res;


%% get dimensions of the SCOIL

Xc_bb = zeros(3,2);

Xc_bb(1,1) = min(xc);  Xc_bb(1,2) = max(xc);
Xc_bb(2,1) = min(yc);  Xc_bb(2,2) = max(yc);
Xc_bb(3,1) = min(zc);  Xc_bb(3,2) = max(zc);

%% get bounding box for SCOIL + RHBM

Xext_bb = zeros(3,2);

%% define a bounding box for the SCOIL + RHBM region

for i = 1:size(Xext_bb,1)

    % case 1: xb_min >= xc_min or xb_max <= xc_max
    % case 2: yb_min >= yc_min or yb_max <= yc_max
    
    if(Xb_bb(i,1) < Xc_bb(i,1) - d_coil2bound)
        Xext_bb(i,1) = Xb_bb(i,1)- eps;
    else
        n_shift = round(abs(Xb_bb(i,1) - Xc_bb(i,1)) / res);
        Xext_bb(i,1) = Xb_bb(i,1) - res * (n_shift + N_exp_half) - eps;
    end
    
    
     if(Xb_bb(i,2) > Xc_bb(i,2) + d_coil2bound)
        Xext_bb(i,2) = Xb_bb(i,2) + eps;
     else
        n_shift = round(abs(Xb_bb(i,2) - Xc_bb(i,2)) / res);
        Xext_bb(i,2) = Xb_bb(i,2) + res * (n_shift + N_exp_half) + eps;
     end


end



