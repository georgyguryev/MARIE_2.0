function [ports] = parse_PortSpecs(specs_filename)
% function [ports] = parse_PortSpecs(specs_filename)
% The function parses "specs_filename" text file with ports and loads 
% speccifications and outputs the ports structure
% ------------------------------------------------------------------------
% Input:
%        - specs_filename - name of file that contains ports specification
% Output:
%        - ports structure 
% ------------------------------------------------------------------------
% Author: Georgy Guryev, Cambridge, MA, 2019    

% check if the file exists 
if ~exist(specs_filename, 'file')
    error ('\n No match found for specified file %s \n', specs_filename);
    
% check if the file has a correct extension
elseif ~strcmp(specs_filename(end-3:end),'.txt')
    error ('\n The provided file %s has a wrong format \n', specs_filename);
end

%% open Port specs file 

fid  = fopen(specs_filename);

% get number of ports 
port_line = strsplit(fgetl(fid),':');
N_ports = str2double(port_line{2});

%% init structures

tuning_param  = struct('LE', [], 'val_1', [], 'val_2', [], 'mutual_port', [], ...
                       'tuning_flag',0, 'ini_val', [], 'min_val', [], 'max_val', []);
matching_param = struct('LE_1',[],  'LE_2', [], 'LE_3', [], 'val_1', [], 'val_2', [], 'val_3', []);

port   = struct('id', [], 'type', [],'load_type', [], 'ref_impedance', [], ...
                'tuning_param', tuning_param, 'matching_param',  matching_param);

% final structure - array of ports
ports = repmat(port, N_ports, 1);

counter = 0; 

%% read port section 

PS_token  = fgetl(fid);
PS_format = fgetl(fid);

for i = 1:N_ports
    str = strsplit(fgetl(fid));
    
    ports(i).id   = str2double(str{1});
    ports(i).type = str{2};
    ports(i).load_type = str{3};
    ports(i).ref_impedance = str2double(str{4});
end

% end PS_token
End_PS_token = fgetl(fid);
Space_token  = fgetl(fid);

%% read tuning section

TS_token  = fgetl(fid);
TS_format = fgetl(fid);

End_TS_token = '$End_TuningSection';

str = fgetl(fid);

while ~strcmp(End_TS_token, str)
    
    % str split 
    str = strsplit(str);
        
    % port id
    id = str2double(str{1});
    
    ports(id).tuning_param.LE    = str{2};
    ports(id).tuning_param.val_1 = str2double(str{3});
    ports(id).tuning_param.val_2 = str2double(str{4});
    ports(id).tuning_param.mutual_port = str2double(str{5});
    ports(id).tuning_param.tuning_flag = str2double(str{6});
    ports(id).tuning_param.ini_val = str2double(str{7});
    ports(id).tuning_param.min_val = str2double(str{8});
    ports(id).tuning_param.max_val = str2double(str{9});
    
    counter = counter + 1;
    
    str = fgetl(fid);
end

Space_token  = fgetl(fid);

%% read matching section

MS_token  = fgetl(fid);
MS_format = fgetl(fid);

End_MS_token = '$End_MatchingSection';

str = fgetl(fid);

while ~strcmp(End_MS_token, str)
    
    % str split 
    str = strsplit(str);

    % port id
    id = str2double(str{1});
    
    ports(id).matching_param.LE_1 = str{2};
    ports(id).matching_param.LE_2 = str{3};
    ports(id).matching_param.LE_3 = str{4};
    ports(id).matching_param.val_1 = str2double(str{5});
    ports(id).matching_param.val_2 = str2double(str{6});
    ports(id).matching_param.val_3 = str2double(str{7});

    counter = counter + 1;
    
    str = fgetl(fid);
end

%% sanity check for port number

assert((N_ports == counter ), 'PortSpecs format error: Port count mismatch!'); 
