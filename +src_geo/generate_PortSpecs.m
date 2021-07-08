function generate_PortSpecs(fid, N_ports, port_specs)
% function generate_PortSpecs(fid, N_ports, port_specs)
% The function generates file with port loads based on port_specs 
% ------------------------------------------------------------------------
% Input:
%        - fid     - file id with port-specific info
%        - N_ports - number of ports 
%        - port_specs - structure with port loads to be written to file
% Output:
%        - .txt file with port specifications
% ------------------------------------------------------------------------


% check if the file is opened successfully  
assert( -1 ~= fid, 'Invalid fid!');

% the port_specs structure is not specified
if nargin < 3
    
    tuning_param = struct('init_val', [], 'min_val', [], 'max_val', []);
    
    load   = struct('type', [], 'value', [], 'mutual_val', [],...
        'mutual_port', [],'matching',[], 'tuning', [],...
        'tuning_param', tuning_param);
    
    port   = struct('id', [], 'type', [],'load', load);
    
    % final structure - array of ports
    port_specs = repmat(port, N_ports, 1);
    
    
    
    for i_port = 1:N_ports
        
        port_specs(i_port).id = i_port;
        port_specs(i_port).type = 'feed';
        port_specs(i_port).load.type = 'capacitor';
        port_specs(i_port).load.value = 10e-12;
        port_specs(i_port).load.mutual_val = '[]';
        port_specs(i_port).load.mutual_port = '[]';
        port_specs(i_port).load.matching = 1;
        port_specs(i_port).load.tuning   = 0;
        port_specs(i_port).load.tuning_param.init_val = '[]';
        port_specs(i_port).load.tuning_param.min_val  = '[]';
        port_specs(i_port).load.tuning_param.max_val  = '[]';
        
    end
end


%% init structures


% Write number of ports
fprintf(fid, '%d \n', N_ports);

% 
for i_port = 1:N_ports
    
    % write port id, port type, load type
    fprintf(fid, '%d %s %s \t', port_specs(i_port).id,...
                         port_specs(i_port).type,...
                         port_specs(i_port).load.type);
                     
                     
    % write load specs
    fprintf(fid, '%s %s %s \t', port_specs(i_port).load.value,...
                             port_specs(i_port).load.mutual_val,...
                             port_specs(i_port).load.mutual_port);
                         
    % tuning-matching flags
    fprintf(fid, '%d %d \t', port_specs(i_port).load.matching,...
                          port_specs(i_port).load.tuning);
                      
    % range of optimized values (tuning)                   
    fprintf(fid, '%s %s %s \t', port_specs(i_port).load.tuning_param.init_val,...
                                port_specs(i_port).load.tuning_param.min_val,...
                                port_specs(i_port).load.tuning_param.max_val);

    fprintf(fid, '\n');

end




