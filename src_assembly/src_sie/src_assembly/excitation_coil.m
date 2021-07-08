function [V] = excitation_coil(B_o,index,node,edge,port)
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                        Excitation Vector                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% first_node  = [3 1 2];
% second_node = [2 3 1];
%

% get total number of unknowns for SIE
Nsie = max(index);

% allocate memory for port excitation
V  = zeros(Nsie,1);

% get the number of terminal edges that form a port
Nfeed_terms = length(port.t);


% loop over terminals
for i = 1:Nfeed_terms
    
     % get current terminal edge
     index_edge = port.t(i);
     real_edge  = find(index == index_edge);
     
     % get points that form a terminal
     point1 = edge(1,real_edge);
     point2 = edge(2,real_edge);
     
     % get coordinates of terminal nodes
     rpoint1 = node(:,point1);
     rpoint2 = node(:,point2);
     
     % compute length of the inner edge
     L_edge = norm(rpoint2-rpoint1);
     
     %
     V(index_edge,1) = L_edge;    
end

V = B_o*V;