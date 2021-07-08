function generate_basis(task_settings, dims, scatterer, coil, freq)

tic; 

% get basis settings 
basis = task_settings.basis;
basis_type = basis.Type;

tol = basis.Tolerance;

dom_basis = scatterer.dom_vie;
index_basis = scatterer.index_vie;

% x,y,z coordinates of the domain
xd = dom_basis.x_tensor;
yd = dom_basis.y_tensor;
zd = dom_basis.z_tensor;
res = dom_basis.res;

emu = src_utils.EM_utils(freq);
idxS = index_basis.S_1d;

xs = xd(idxS);
ys = yd(idxS);
zs = zd(idxS);

Scoord = [xs,ys,zs].';

switch basis_type
    
    case 'Dipole'

        % define number of sampling points
        N_samp  = basis.Number_of_points;
        
        % order of basis
        N_l = dims.l;
        
        % order of quadrature
        Quad_order_vie = task_settings.vsie.Np_quad_coup_vie;

        %% generate equivalent surface
        
        r_pts = generate_CSEP_samples(basis, N_samp);
        
        % plot generated points
        figure(); scatter3(r_pts(:,1), r_pts(:,2), r_pts(:,3));
        
        %% get quadrature points
        
        [wp_vie, z_vie] = gauss_1d(Quad_order_vie);
        
        VIE_quads = [Quad_order_vie; wp_vie; z_vie];        
       
        %% compute fields produced by dipoles
        
        Zbc_approx = src_coupling.Assemble_p2v_coupling(Scoord, r_pts, VIE_quads, res, N_l, emu.k0);
        
    case 'RWG'
        profile on;
        
        tic;
%         Zbc_approx = src_coupling.assemble_coupling_matrix(task_settings, dims, scatterer, coil, freq);
        [U_N,U_K] = src_assembly.generate_basis_CS(task_settings, dims, scatterer, coil, freq);
        toc;
        profile off;
        profile viewer;
end

%% apply DEIM to generated basis

tic; 

% check if DEIM++ is available
if (3 == exist('mexDeim','file'))
    [~, Pin] = mexDeim(U_N);
else
    [~, Pin] = src_numeric.deim(U_N);
end

toc;
%% postprocessing 

% get S domain coordinates
nS = size(xd(idxS),1);

% find all the voxels with at least one component selected by deim
Pin = reshape(Pin, nS, dims.ql * idx);
idxD = find(sum(Pin,2));

% get deim points coordinates
xds = xd(idxS(idxD));
yds = yd(idxS(idxD));
zds = zd(idxS(idxD));

Dcoord = [xds(:), yds(:), zds(:)];

%  build the Pin for the DEIM approach
idxD= repmat(idxD, dims.ql,1) + kron(nS*[0:dims.ql-1].', ones(size(idxD)));

e   = speye(dims.ql * nS);
Pin = sparse(e(:, idxD));
ndeim = size(Pin,2);


% Obtain the DEIM extended inverse matrix
Xin_N = (Pin.'*U_N)\speye(ndeim,ndeim);
Xin_K = (Pin.'*U_K)\speye(ndeim,ndeim); 

%%

save(fullfile(basis.Path,basis.Filename),...
    'Pin', 'Xin_N', 'Xin_K', 'U_N', 'U_K', 'Dcoord','dom_basis', 'index_basis','-v7.3');

Time_Basis_Assembly = toc;

fprintf(fid, '\n Time spent on generating basis:            %.2f sec', Time_Basis_Assembly);



