function E_f = interpolate_fields(Ein, scatterer, dims, ratio)
% function interpolate_fields(E, dims, ratio) 
% interpolates fields on a finer grid, converting fields tested with linear
% basis to constant components on a finer grid

E = zeros(dims.vie);

idxS = scatterer.index_vie.index_ql(dims.ql, dims.Nvox_vie);

E(idxS) = Ein;


[L,M,N,~] = size(E);

L_f = zeros(L,ratio);
M_f = zeros(M,ratio);
N_f = zeros(N,ratio);

for i = 1:ratio
    L_f(:,i) = i:ratio:ratio*L-(ratio-i);
    M_f(:,i) = i:ratio:ratio*M-(ratio-i);
    N_f(:,i) = i:ratio:ratio*N-(ratio-i);
end

E_f = zeros(ratio*L,ratio*M,ratio*N,3);

for l = 1:ratio
    for m = 1:ratio
        for n = 1:ratio
            
            l_step = -(ratio-(2*l-1))/(ratio);
            m_step = -(ratio-(2*m-1))/(ratio);
            n_step = -(ratio-(2*n-1))/(ratio); 
            
            for i = 1:3
                
                E_f (L_f(:,l),M_f(:,m),N_f(:,n),i) = E(:,:,:,(i-1)*4+1) / ratio^3 + ...
                                                     l_step * E(:,:,:,(i-1)*4+2) + ...
                                                     m_step * E(:,:,:,(i-1)*4+3) + ...
                                                     n_step * E(:,:,:,(i-1)*4+4);
                                                 
            end
            
        end
    end
end