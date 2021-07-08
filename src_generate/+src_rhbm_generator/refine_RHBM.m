function RHBM_refined = refine_RHBM(RHBM, factor)

% get underlying tissue parameters
r         = RHBM.r;
epsilon_r = RHBM.epsilon_r;
sigma_e   = RHBM.sigma_e;
rho       = RHBM.rho;

% get problem dimensions
[L,M,N] = size(epsilon_r);

old_res = abs(r(2,1,1,1) - r(1,1,1,1));
new_res = old_res / factor; 

epsilon_r_new = zeros(L * factor, M * factor, N * factor);
sigma_e_new   = zeros(L * factor, M * factor, N * factor);
rho_new       = zeros(L * factor, M * factor, N * factor);

% generate new cartesian grid
x_new = zeros(factor*L,1);
y_new = zeros(factor*L,1);
z_new = zeros(factor*L,1);

L_f = zeros(L,factor);
M_f = zeros(M,factor);
N_f = zeros(N,factor);

for i = 1:factor
    x_new(i:factor:factor*L) = squeeze(r(:,1,1,1)) + (new_res - old_res)/2  + new_res * (i-1);
    y_new(i:factor:factor*M) = squeeze(r(1,:,1,2)) + (new_res - old_res)/2  + new_res * (i-1);
    z_new(i:factor:factor*N) = squeeze(r(1,1,:,3)) + (new_res - old_res)/2  + new_res * (i-1);
end

r_new = src_geo.grid3d(x_new, y_new, z_new);


% refine material properties
for i = 1:factor
    L_f(:,i) = i:factor:factor*L-(factor-i);
    M_f(:,i) = i:factor:factor*M-(factor-i);
    N_f(:,i) = i:factor:factor*N-(factor-i);
end

for l = 1:factor
    for m = 1:factor
        for n = 1:factor
            
%             l_step = -(stepi-(2*l-1))/(2*stepi);
%             m_step = -(stepi-(2*m-1))/(2*stepi);
%             n_step = -(stepi-(2*n-1))/(2*stepi);
            

                
            epsilon_r_new(L_f(:,l),M_f(:,m),N_f(:,n)) = epsilon_r;
            sigma_e_new(L_f(:,l),M_f(:,m),N_f(:,n))   = sigma_e;
            rho_new(L_f(:,l),M_f(:,m),N_f(:,n))       = rho;

            
        end
    end
end

% form a new structure
RHBM_refined.r = r_new;
RHBM_refined.epsilon_r = epsilon_r_new;
RHBM_refined.sigma_e = sigma_e_new;
RHBM_refined.rho = rho_new;
