function [Near_RWG]  = form_RWG_near_list(SCOIL, RHBM, Ncell, i_sie, idx_x_ctr, idx_y_ctr, idx_z_ctr) 
% function []  = form_RWG_near_list(SCOIL, RHBM, idx_x_ctr, idx_y_ctr, idx_z_ctr) 
% function creates a list of all triangles that are enclosed within near
% zone


% copy SIE and VIE data to local variables
node  = SCOIL.node;
elem  = SCOIL.elem;
etod  = SCOIL.etod;
index = SCOIL.index;
edge  = SCOIL.edge;


% get extended VIE grid coordinates
r = RHBM.r;

[Lx,My,Nz,~] = size(r);

x_vie = squeeze(r(:,1,1,1));
y_vie = squeeze(r(1,:,1,2));
z_vie = squeeze(r(1,1,:,3));


% get total number of SIE dofs
N_sie = max(index);


% preallocate sparse vector of near triangles
Near_RWG_buf = zeros(1, N_sie);


% get "near"-zone limits


% we assume that Ncell is n_1d^3 , where n_1d is number of cubes per Cell
% for now the total number of near expansion voxels 1D is given by 
n_1d     =  ceil(Ncell^(2/3));

% indices of exterior points
idx_min  = -(n_1d - 1) / 2;
idx_max  =  (n_1d - 1) / 2;

idx_x_min = max(idx_min + idx_x_ctr,1);
idx_y_min = max(idx_min + idx_y_ctr,1);
idx_z_min = max(idx_min + idx_z_ctr,1);
idx_x_max = min(idx_max + idx_x_ctr,Lx);
idx_y_max = min(idx_max + idx_y_ctr,My);
idx_z_max = min(idx_max + idx_z_ctr,Nz);


% define "near" zone limits 
x_min    = x_vie(idx_x_min);
y_min    = y_vie(idx_y_min);
z_min    = z_vie(idx_z_min);
x_max    = x_vie(idx_x_max);
y_max    = y_vie(idx_y_max);
z_max    = z_vie(idx_z_max);


x = x_vie(idx_x_min:idx_x_max);
y = y_vie(idx_y_min:idx_y_max);
z = z_vie(idx_z_min:idx_z_max);


r_near = grid3d(x,y,z);

x_near = r_near(:,:,:,1); 
y_near = r_near(:,:,:,2);
z_near = r_near(:,:,:,3);

% plot3(x_near(:),y_near(:),z_near(:),'.r');

% if any of 4 points is within the limits and RWG to the "Near" list

% loop over RWG basises
for i = 1:N_sie
    
    % by default current RWG is far
    flag_near = 0;
    
    % get idx of current physical edge
    phys_edge_idx = find(index == i);
    
    % get parent triangles 
%     first_parent_tr  = ceil(find(etod == phys_edge_idx) / 3 );
%     secnd_parent_tr  = ceil(find(etod == -phys_edge_idx) / 3 );
%     
%     % parent triangles have coordinates of 4 unique vertices
%     
%     first_parent_nodes = elem(1:3,first_parent_tr);
%     secnd_parent_nodes = elem(1:3,secnd_parent_tr);
%     
%     cur_nodes_list = unique([first_parent_nodes; secnd_parent_nodes]);
%     
%     flag_near = 0;
%     
%     % check if any vecrtice belongs to 
%     for cur_node = 1:length(cur_nodes_list) 
%         cur_r = node(:,cur_nodes_list(cur_node));
%         
%         if (  (cur_r(1) >= x_min && cur_r(1) <= x_max)...
%            && (cur_r(2) >= y_min && cur_r(2) <= y_max)...
%            && (cur_r(3) >= z_min && cur_r(3) <= z_max))
%        
%             flag_near = 1;
%             break;
%             
%         end;
%         
%     end;


    % get points that form current basis centre
    point_1 = edge(1,phys_edge_idx);
    point_2 = edge(2,phys_edge_idx);
    
    % get coordinates of central vertices
    r_1 = node(:,point_1);
    r_2 = node(:,point_2);

    % get coordinates of center for given basis function
    r_c = (r_2 + r_1) ./ 2;
    
    if ((r_c(1) >= x_min && r_c(1) <= x_max)...
     && (r_c(2) >= y_min && r_c(2) <= y_max)...
     && (r_c(3) >= z_min && r_c(3) <= z_max))
       
         flag_near = 1;
         
    end
    
    % 
    if (1 == flag_near)
        Near_RWG_buf(i) = 1.0;
    end
    
end

% get indicies of near RWG basis functions
Near_RWG_buf = (find(Near_RWG_buf)).';

% form pairs of near interations
Near_RWG = [repmat(i_sie,size(Near_RWG_buf)), Near_RWG_buf];



