function [Z] = lumped_elem_assembly(Z,index,node,edge,port,f)

% setup EM constants
omega = 2 * pi * f;

% get number of terminals corresponding to given port
Nterms = length(port.t);

%% select load type

switch port.tuning_param.LE
    
    % Resistor
    case 'resistor'
        Zx_val = port.tuning_param.val_1;
   
    % Inductor
    case 'inductor'
        Zx_val = 1i * omega * port.tuning_param.val_1;
        
    % Capacitor
    case 'capacitor'
        Zx_val = 1 / (1i * omega * port.tuning_param.val_1);
        
    case 'mutual_inductor'
        Zx_val   = 1i * omega * port.tuning_param.val_1;
        Zx_M_val = 1i * omega * port.tuning_param.val_2;
        
    % Unknown load type
    otherwise
        error('The provided load type is unknown');
end

%% update impedance matrix

L_edges = zeros(Nterms,1);

% compute legths of all terminals, that form given port
for i = 1:Nterms
    
    % get current terminal
    index_edge = port.t(i);

    % find corresponding physical edge
    phys_edge = find(index == index_edge);
    
    % get corresponding terminal nodes 
    point1 = edge(1,phys_edge);
    point2 = edge(2,phys_edge);
    
    % get nodes' coordinates
    rpoint1 = node(:, point1);
    rpoint2 = node(:, point2);
   
    % find current terminal edge length
    L_edges(i) = norm(rpoint1-rpoint2);
end

% compute impedance update values 
for i = 1:Nterms
    
    % get row index that corresponds to i-th terminal
    index_i = port.t(i);
    
    for j = 1:Nterms
        
        % get row index that corresponds to i-th terminal
        index_j = port.t(j);
        
        % compute update impedance value
        Z_update = Zx_val * L_edges(i) * L_edges(j);
        
        % update resulting impedance matrix
        Z(index_i, index_j) = Z(index_i, index_j) - Z_update;
    end
end

%% add contribution of the mutual inductors (if any)
if strcmp(port.tuning_param.LE,'mutual_inductor')
    
    Nterm_mutual = length(port.t_mutual);
    Lm_edges     = zeros(Nterm_mutual,1);
    
    % 
    for i = 1:Nterm_mutual
        
        index_edge = port.t_mutual(i);
        phys_edge = find(index == index_edge);
        
        point1 = edge(1,phys_edge);
        point2 = edge(2,phys_edge);
        
        rpoint1 = node(:, point1);
        rpoint2 = node(:, point2);
        
        Lm_edges(i) = norm(rpoint1-rpoint2);
    end
    
    
    % add contribution of mutual inductancies
    for i = 1:Nterms
        index_i = port.t(i);
        for j = 1:Nterm_mutual
            index_j = port.t_mutual(j);
            Z_update = Zx_M_val * L_edges(i) * Lm_edges(j);
            Z(index_i, index_j) = Z(index_i, index_j) - Z_update;
        end
    end
end