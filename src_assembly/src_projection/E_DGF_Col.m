function    [E] =  E_DGF_Col(SCOIL, i_sie, r_col, GL_order, freq, res) 
% function  [E] =  E_DGF_Col(SCOIL, r_col, GL_order, f)  
% this function computes three components of the electric fields 
% created by SCOIL coil model at collocation locations r_col
%
% -----------------------------------------------------------------------
%                      CAUTION! 
% -----------------------------------------------------------------------
% The resulting vector E is the actual electric field 
% The field components in vector E are stored in format: 
%
%   E = [E_1x, E_1y, E_1z, ...,E_nx, E_ny, E_nz].'         
%
% -----------------------------------------------------------------------

node  = SCOIL.node;
elem  = SCOIL.elem;
edge  = SCOIL.edge;
etod  = SCOIL.etod;
index = SCOIL.index;


%% get size of output vector
[rows,cols] = size(r_col);

% M = rows * cols;

E = zeros(cols,rows);

%% EM constants
EM_params;
% mu = 4*pi*1e-7;
% co = 299792458;
% eo = 1/co^2/mu;
% 
% % frequency-dependend parameters
% lambda = co / f;
% omega  = 2 * pi * f;
% ko     = 2 * pi / lambda; 

%% get SIE points 

[Np,wt,Z,z1,z2,z3] = gauss_2d_sie(GL_order);

%% preprocess COIL data

% find physical index of common edge
phys_edge_idx = find(index == i_sie);

% find parent triangles
[idx_edge_1,parent_el_1] = find(etod == phys_edge_idx); 
[idx_edge_2,parent_el_2] = find(etod == -phys_edge_idx); 

% put supporting triangles into the list
parent_els = [parent_el_1, parent_el_2];
idx_edges  = [idx_edge_1, idx_edge_2];

%% this part is developed specificaly for a single RWG sopport

% length of internal edge: l = norm(r_3 - r_2);
l = norm(node(:,edge(1,phys_edge_idx)) - node(:,edge(2,phys_edge_idx)));

scaling_factor = 1 / ce * l;

for point = 1:cols
    
    % get coordinates of current collocation point 
    r_cur = r_col(:,point);
    
    % create current electric field components, init with zeros
    E_cur = zeros(3,1);
    
    % Loop over RWG supporting triangles 
    for sup = 1:2
        
        parent_el = parent_els(sup);
        idx_edge  = idx_edges(sup);
        
        % get vertex of current supporting triangle
        r_1 = node(:,elem(1,parent_el));     
        r_2 = node(:,elem(2,parent_el));
        r_3 = node(:,elem(3,parent_el));
        r_v = node(:,elem(idx_edge,parent_el));
        
        % sign of common edge vector
        sgn = (-1)^(sup + 1);
        
        % Integrate RWG support
        for i_sie = 1:Np
            
            % get coordinates of current source 
            r_src = r_1 .* z1(i_sie) + r_2 .* z2(i_sie) + r_3 .* z3(i_sie);
            
            % find vector and distance from current source point to
            % collocation point
            R_vec =  r_src - r_cur;
            R     =  norm(R_vec);
            
            % define multiplication coeeficients
            kappa = exp(-1i * ko * R) / (4 * pi * R^3);
            P     = (3 + 3 * 1i * ko * R - ko^2 * R^2) / R^2;  
            Q     = (1 + 1i * ko * R - ko^2 * R^2) / P;
            mult  = kappa * P;
            
            % compute DGF entries
            Gxx   = R_vec(1) * R_vec(1) - Q;
            Gxy   = R_vec(1) * R_vec(2);
            Gxz   = R_vec(1) * R_vec(3);
            Gyy   = R_vec(2) * R_vec(2) - Q;
            Gyz   = R_vec(2) * R_vec(3);
            Gzz   = R_vec(3) * R_vec(3) - Q;
            
            % form dyadic Green's function
            DGF = [Gxx, Gxy, Gxz;
                   Gxy, Gyy, Gyz;
                   Gxz, Gyz, Gzz];
           
           
            E_cur = E_cur + sgn .* wt(i_sie) .* mult .* DGF * (r_src - r_v);
        end
    end
    
    % multiply by scaling factor
%     E_cur =  E_cur; 
    
    E(point,:) = E_cur;
%     E(point,2) = E_cur(2);
%     E(point,3) = E_cur(3);

end

E = reshape(E,[rows * cols,1]);

E = scaling_factor * E;

% E = E * res^3;

