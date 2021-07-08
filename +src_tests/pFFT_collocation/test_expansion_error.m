close all;

%% constants and parameters 
vertex_labels = {'Vp', 'Vn', 'V3', 'V4'};
patches = {[1 3 4], [2 3 4]};

tr_scale = 0.0025 * [1, -1, 0.5, 0.5];

% generate random pair of triangles
r_tr = tr_scale .* rand(3,4);
r_o = (r_tr(:,3) + r_tr(:,4)) / 2;

[Xs,Ys,Zs] = sphere; 

L = [norm(r_tr(:,1) - r_tr(:,3)), norm(r_tr(:,1) - r_tr(:,4)), ...
     norm(r_tr(:,3) - r_tr(:,4)), norm(r_tr(:,2) - r_tr(:,3)), ...
     norm(r_tr(:,2) - r_tr(:,4))];



freq = 298e6;

Quad_order_C2Col = 1;
Quad_order_sie   = 10;
Quad_order_vie   = 3;

% get quadrature points and weights
[wp_vie, z_vie]                          = gauss_1d(Quad_order_C2Col);
[Np_sie, z1_sie, z2_sie, z3_sie, wp_sie] = dunavant_rule(Quad_order_sie);

% form a vector of quadratures
VIE_quads = [Quad_order_C2Col; wp_vie; z_vie];
SIE_quads = [Np_sie; wp_sie; z1_sie; z2_sie; z3_sie];

%% 

res = 0.003;


dims.ql = 3;
dims.q  = 3;
dims.l  = 1;

% expansion voxels along 1D axis
N_exp = 5;
N_R_near = 3 * max(ceil(N_exp/2), ceil(max(L/res)));
N_dist = 10 * N_R_near;

% get electromagnetic constants
emu = src_utils.EM_utils(freq);



%% 


switch N_exp
    % select number of quad points on a sphere
    case 3
        N_colloc_pts = 26;
    case 5
        N_colloc_pts = 86;
    case 7
        N_colloc_pts = 230;
    otherwise
        N_colloc_pts = 230;
end

N_colloc_pts = 86;

% generate collocation points 
col_points = src_assembly.generate_collocation_points(res, N_R_near, N_colloc_pts);

% generete cube centers for elementary Cell
[vox_centers] = src_assembly.generate_cube_centers(r_o, res, N_exp);

% compute Lop/Nop matrix that maps volumetrix currents to E in colloc. points      
[A] = E_DGF_PWX_Col(vox_centers.', res, col_points, Quad_order_vie, freq, dims);

% form current collocation points
cur_col_points = col_points + r_o;  % + repmat(r_c, 1, size(col_points,2));

% get fields produced by unitary excitation of current rwg triangle
rhs = Assemble_rwg_coupling(cur_col_points, r_tr, SIE_quads, VIE_quads, res, dims.l, emu.k0);

% find weigths/ scale factors of volumetric currents
Ic = A \ rhs;

norm(rhs - A * Ic) / norm(rhs)

 %% test accuracy as a function of distance 

N_colloc_pts = 86;
rel_error    = cell(N_colloc_pts,1);
abs_error    = cell(N_colloc_pts,1);
for i = 1:N_dist 
    
    % generate collocation points
    col_points = src_assembly.generate_collocation_points(res, i, N_colloc_pts);
    
    % generete cube centers for elementary Cell
    [vox_centers] = src_assembly.generate_cube_centers(r_o, res, N_exp);
    
    % compute Lop/Nop matrix that maps volumetrix currents to E in colloc. points
    [A_test] = E_DGF_PWX_Col(vox_centers.', res, col_points, Quad_order_vie, freq, dims);
    
    % form current collocation points
    cur_col_points = col_points + r_o;  % + repmat(r_c, 1, size(col_points,2));
    
    % get fields produced by unitary excitation of current rwg triangle
    rhs_test = Assemble_rwg_coupling(cur_col_points, r_tr, SIE_quads, VIE_quads, res, dims.l, emu.k0);
    
    rel_error{i} = norm(rhs_test - A_test * Ic) / norm(rhs_test);
    abs_error{i} = norm(rhs_test - A_test * Ic) / norm(rhs);
end


%% visualization


figure(); 

Xs = res * N_R_near * Xs;
Ys = res * N_R_near * Ys;
Zs = res * N_R_near * Zs;

% surf(Xs,Ys,Zs);
% alpha 0.01
% hold on;

scatter3(vox_centers(:,1), vox_centers(:,2), vox_centers(:,3));
hold on;
scatter3(r_tr(1,:), r_tr(2,:), r_tr(3,:));

for i=1:4
    text(r_tr(1,i), r_tr(2,i), r_tr(3,i), vertex_labels{i});
end

for i = 1:2
    patch(r_tr(1,patches{i}), r_tr(2,patches{i}), r_tr(3,patches{i}), 'g');
end

figure(); semilogy([rel_error{:}]);



