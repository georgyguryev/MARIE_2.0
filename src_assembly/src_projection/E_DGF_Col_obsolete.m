function    [E] =  E_DGF_Col_obsolete(SCOIL, r_col, GL_order, f) 
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

node = SCOIL.node;
elem = SCOIL.elem;
edge = SCOIL.edge;


E = [];

%% get size of output vector
[rows,cols] = size(r_col);

M = rows * cols;

%% EM constants
mu = 4*pi*1e-7;
co = 299792458;
eo = 1/co^2/mu;

% frequency-dependend parameters
lambda = co / f;
omega  = 2 * pi * f;
ko     = 2 * pi / lambda; 

%% get SIE points 

[Np,wt,Z,z1,z2,z3] = gauss_2d_sie(GL_order);

%% this part is developed specificaly for a single RWG sopport

% length of internal edge: l = norm(r_3 - r_2);
l = norm(node(:,edge(1,1)) - node(:,edge(1,2)));

for point = 1:cols
    
    % get coordinates of current collocation point 
    r_cur = r_col(:,point);
    
    % create current electric field components, init with zeros
    E_cur = zeros(3,1);
    
    % Loop over RWG supporting triangles 
    for sup = 1:2
        
        % get vertex of current supporting triangle
        r_1 = node(:,elem(1,sup));     
        r_2 = node(:,elem(2,sup));
        r_3 = node(:,elem(3,sup));
        
        % sign of common edge vector
        sgn = (-1)^(sup + 1);
        
        % Integrate RWG support
        for i_sie = 1:Np
            
            % get coordinates of current source 
            r_src = r_1 .* Z(i_sie, 1) + r_2 .* Z(i_sie, 2) + r_3 .* Z(i_sie, 3);
            
            % find vector and distance from current source point to
            % collocation point
            R_vec =  r_src - r_cur;
            R     =  norm(R_vec);
            
            % define multiplication coeeficients
            kappa = exp(-1i * ko * R) / (4 * pi * R^3);
            P     = (3 + 3 * 1i * ko * R - ko^2 * R^2) / R^2;  
            Q     = (ko^2 * R^2  - 1i * ko * R - 1) / P;
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
            
            
           
            E_cur = E_cur + sgn .* wt(i_sie) .* mult .* DGF * (r_src - r_1);
            
        end;
    end;
    
    % multiply by scaling factor
    E_cur = 1i * omega * 1/(ko^2)* mu * l * E_cur; 
    
    E = [E; E_cur];

end;


