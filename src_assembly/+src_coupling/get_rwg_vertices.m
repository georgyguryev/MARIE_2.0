function r_rwg = get_rwg_vertices (coil, i_sie)
% function get_rwg_vertices(coil, i_sie)
% extracts vertices that form the support of the i_sie-th RWG basis function
% INPUT:
%        - coil object
%        - i_sie - index of the RWG basis function
% OUTPUT:
%        - coordinates of vertices of the i_sie-th support
% Author: Georgy Guryev, Cambridge, MA, 2019    


% find physical edge for given dof(i_sie)
phys_edge = find(coil.index == i_sie);

% find tirangles that share this physical edge
[idx_tr1,tr_1] = find(coil.etod == phys_edge);
[idx_tr2,tr_2] = find(coil.etod == -phys_edge);

% get coordinates of oposit vertices
rv_p = coil.node(:,coil.elem(idx_tr1,tr_1)); 
rv_n = coil.node(:,coil.elem(idx_tr2,tr_2));

% get local index of vertices, that form common edge
edge_node = setdiff([1,2,3], idx_tr1);

% get coordinates of vertices, that form common edge
r23 = coil.node(:,coil.elem(edge_node,tr_1));

% disp(i_sie);

% stack coordinates in a single vector
r_rwg = [rv_p; rv_n; r23(:)];

