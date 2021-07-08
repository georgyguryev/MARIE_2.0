function [Q] = first_moment_RWG(RWG, GL_order)
% function [Q] = first_moment_RWG(RWG, GL_order)
% function computes first order RWG moments for 


Q = zeros(9,1);


%%
node = RWG.node;
elem = RWG.elem;
edge = RWG.edge;
ct   = RWG.Ct;

% length of internal edge: l = norm(r_3 - r_2);
l = norm(node(:,edge(1,1)) - node(:,edge(1,2)));


%% quadrature points and weights for SIE
[Np,wt,Z,~,~,~] = gauss_2d_sie(GL_order);


%% Compute Current components by integrating over RWG support
for sup = 1:2
    
    % get vertex of current supporting triangle
    r_1 = node(:,elem(1,sup));
    r_2 = node(:,elem(2,sup));
    r_3 = node(:,elem(3,sup));
    
    % sign of common edge vector
    sgn = (-1)^(sup + 1);
    
    for i_sie = 1:Np

        r_src = r_1 .* Z(i_sie, 1) + r_2 .* Z(i_sie, 2) + r_3 .* Z(i_sie, 3);
        
        Rho   = r_src - r_1; 
        dr    = r_src - ct;
            
        Q = Q + sgn .* wt(i_sie) .* kron(Rho, dr);
    end;
end;

Q = Q .* l; 