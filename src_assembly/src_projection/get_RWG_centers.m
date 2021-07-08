function [RWG_cntr] = get_RWG_centers(SCOIL)
% function [RWG_cntr] = get_RWG_centers(SCOIL)
% forms a matrix of coordinates of centers of 
% teminal and inner edges (dofs)

% copy vars 
index = SCOIL.index;
edge  = SCOIL.edge;
node  = SCOIL.node;

% get number of non-boundary edges
N_dofs = max(index);

% allocate memory for matrix of center coordinates
RWG_cntr = zeros(N_dofs, 3); 


for i_sie = 1:N_dofs
    
    % get index of physical edge associated with current dof edge
    
    phys_idx = find(index == i_sie);
    
    % get indexes of boundary vertices
    node1 = edge(1, phys_idx);
    node2 = edge(2, phys_idx);
    
    r1   = node(:,node1);
    r2   = node(:,node2);
    
    % get coordinates of edge center
    r_c   = (r1(:) + r2(:)) ./ 2;
    
    % fill in the matrix
    RWG_cntr(i_sie, :) = r_c.';
end

