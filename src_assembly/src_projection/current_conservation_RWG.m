function [I] = current_conservation_RWG(SCOIL, GL_order, i_sie)
% function [I] = current_conservation_RWG(SCOIL, GL_order)
% computes three components of net current produced by single RWG support

% allocate vector forn net current components
I = zeros(3,1);

%% extract RWG parameters

node  = SCOIL.node;
elem  = SCOIL.elem;
edge  = SCOIL.edge;
index = SCOIL.index;
etod  = SCOIL.etod;

phys_edge = find(index == i_sie);

[idx_edge_1, fisrt_parent_tr]  = find(etod == phys_edge);
[idx_edge_2, second_parent_tr] = find(etod == -phys_edge);

% store parent triangles in a vector
parent_tr = [fisrt_parent_tr, second_parent_tr];
idx_edges = [idx_edge_1, idx_edge_2];

% length of internal edge: l = norm(r_3 - r_2);
l = norm(node(:,edge(1,phys_edge)) - node(:,edge(2,phys_edge)));


%% quadrature points and weights for SIE
[Np,wt,Z,z1,z2,z3] = gauss_2d_sie(GL_order);


%% Compute Current components by integrating over RWG support
for sup = 1:2
    
    % 
    current_tr      = parent_tr(sup);
    current_vertice = idx_edges(sup);
    
    % get vertex of current supporting triangle
    r_1 = node(:,elem(1,current_tr));
    r_2 = node(:,elem(2,current_tr));
    r_3 = node(:,elem(3,current_tr));
    r_v = node(:,elem(current_vertice, current_tr));
    
    % sign of common edge vector
    sgn = (-1)^(sup + 1);
    
    for j_quad = 1:Np

        r_src = r_1 .* Z(j_quad, 1) + r_2 .* Z(j_quad, 2) + r_3 .* Z(j_quad, 3);
        I = I + sgn .* wt(j_quad) .* (r_src - r_v);
    end
end

I = I .* l; 
