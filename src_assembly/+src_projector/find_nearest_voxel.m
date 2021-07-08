function [idx_nv] = find_nearest_voxel(r_sie, scatterer)
% function [i_vie, j_vie, k_vie] = find_nearest_voxel(r_sie, RHBM)
% function finds the nearest voxel to given position r_sie
% ------------------------------------------------------------------------
% INPUT:
%           r_sie - position of current RWG edge center
%           RHBM  - Extended voxelized domain VIE + SIE 
% ------------------------------------------------------------------------
% OUTPUT:
%           idx_nv = [i_vie; j_vie; k_vie] - tuple of indices of the nearest voxel     
% ------------------------------------------------------------------------


% get 1D coordinates of the domain
xd = scatterer.dom_ext.x;
yd = scatterer.dom_ext.y;
zd = scatterer.dom_ext.z;


% get size of cube side
res = scatterer.dom_ext.res; 

% sanity check: r_sie has to be inside extended region!!!
if (( min(xd) > r_sie(1)) || (max(xd) < r_sie(1)) ||...
    ( min(yd) > r_sie(2)) || (max(yd) < r_sie(2)) || ...
    ( min(zd) > r_sie(3)) || (max(zd) < r_sie(3)))

    error('find_nearest_voxel(): the r_sie value does not belong to the VSIE region');
end


%%

% form a list of potential candidates
[list_idx_x] = find((xd >= (r_sie(1) - res)) & (xd <= (r_sie(1) + res)));
[list_idx_y] = find((yd >= (r_sie(2) - res)) & (yd <= (r_sie(2) + res)));
[list_idx_z] = find((zd >= (r_sie(3) - res)) & (zd <= (r_sie(3) + res)));


for i_x = 1:length(list_idx_x)
    for j_y = 1:length(list_idx_y)
        for k_z = 1:length(list_idx_z)
            
            dist = norm(r_sie - [xd(list_idx_x(i_x)); yd(list_idx_y(j_y)); zd(list_idx_z(k_z))],Inf);
            if dist - res/2 < eps
                % return computed indices
                i_vie = list_idx_x(i_x);
                j_vie = list_idx_y(j_y);
                k_vie = list_idx_z(k_z);
                
                idx_nv = [i_vie; j_vie; k_vie];
                return;
            end
        end
    end
end

