function [dom_ext, prop_ext] = extend_domain(coil, scatterer, projector, dims, freq)


% copy RHBM coordinated to vectors, get resolution
x  = scatterer.dom_vie.x;
y  = scatterer.dom_vie.y;
z  = scatterer.dom_vie.z;

res = scatterer.dom_vie.res;

%% generate extended domain
%  produce new extended grid for the bounding box
%  new grid is genereted with respect to the VIE voxelization

emu = src_utils.EM_utils(freq);

ge = @(a,b) a >= (b - eps); 
le = @(a,b) a <= (b + eps);

% get bounding box for SCOIL + RHBM domain
[Xext_bb, ~, Xb_bb] = src_scatterer.get_BB_limits(coil, scatterer, projector, dims);

% get indices withing the bounding box
idx_body_x = find(ge(x,Xb_bb(1,1)) & le(x, Xb_bb(1,2)));
idx_body_y = find(ge(y,Xb_bb(2,1)) & le(y, Xb_bb(2,2)));
idx_body_z = find(ge(z,Xb_bb(3,1)) & le(z, Xb_bb(3,2)));


% generate coordinates of extended domain
dom_ext = src_scatterer.Domain(generate_extended_domain(Xext_bb, res));

% update dims sizes with extended domain dimensions
dims.get_extended_body_dims(dom_ext.r);


% allocate indices for material properties
prop_ext.epsilon_r = ones(dims.ext(1:3));
prop_ext.sigma_e   = zeros(dims.ext(1:3));
prop_ext.rho       = zeros(dims.ext(1:3));

% get bounding box for body in extended domain
xb_min_new = max(Xb_bb(1,1),min(x)) - res / 3;
xb_max_new = min(Xb_bb(1,2),max(x)) + res / 3;
yb_min_new = max(Xb_bb(2,1),min(y)) - res / 3;
yb_max_new = min(Xb_bb(2,2),max(y)) + res / 3;
zb_min_new = max(Xb_bb(3,1),min(z)) - res / 3;
zb_max_new = min(Xb_bb(3,2),max(z)) + res / 3;


% mapping indexes from extended to VIE domain
idx_vie_x = find((dom_ext.x >= xb_min_new) & ...
                 (dom_ext.x <= xb_max_new));
idx_vie_y = find((dom_ext.y >= yb_min_new) & ...
                 (dom_ext.y <= yb_max_new));
idx_vie_z = find((dom_ext.z >= zb_min_new) & ...
                 (dom_ext.z <= zb_max_new));

% map material properties to the extended grid
prop_ext.rho(idx_vie_x,idx_vie_y,idx_vie_z)        = scatterer.prop_vie.rho(idx_body_x,idx_body_y,idx_body_z);
prop_ext.sigma_e (idx_vie_x,idx_vie_y,idx_vie_z)   = scatterer.prop_vie.sigma_e(idx_body_x,idx_body_y,idx_body_z);
prop_ext.epsilon_r (idx_vie_x,idx_vie_y,idx_vie_z) = scatterer.prop_vie.epsilon_r(idx_body_x,idx_body_y,idx_body_z);
 
