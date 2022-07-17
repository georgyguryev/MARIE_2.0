function [F,feed_tune_ports] = Assembly_SIE_excitation(SCOIL) 
%%    Assembly the SIE system (parallel version)
% _________________________________________________________________________
%
%   Fucntion to generate the free-space Z matrix for the SIE
%   RWG functions are applied
%
% _________________________________________________________________________
%
%
%% INPUT
%       SCOIL structure with
%           index - mapping of the internal edge number to dof number
%           etod - etod numbering of elements and edges
%           node - coordinates of the nodes 
%           edge - numbering of edges
%           elem - 3 indexes of the nodes defining an element
%           Index_elem - mapping of index to elements
%           port - port definition
%       freq - frequency
%
%
%% OUTPUT
%       Z is the impedance matrix
%       F is the rhs for the port definition
%
%
% -------------------------------------------------------------------------
%
%   J. Fernandez Villena -- jvillena@mit.edu
%   A.G. Polimeridis -- thanos_p@mit.edu
%   Computational Prototyping Group, RLE at MIT
%
% _________________________________________________________________________


% extract struct data
index = SCOIL.index;
node = SCOIL.node;
edge = SCOIL.edge;
port = SCOIL.port;


% -------------------------------------------------------------------------
%            Info 
% -------------------------------------------------------------------------

Nvars  = max(index);
Nfeeds = nnz(strcmp({port(:).type}, 'feed'));

% Data
Eo     = 1;
B_o    = -Eo ;   % - sign because e^sca = -e^inc

F = zeros(Nvars,Nfeeds);

%% form list of feed and load port numbers

port_ids  = [port.id];
tuning_param = [port.tuning_param];

% create feed/load mask
feed_mask = strcmp({port(:).type}, 'feed');
tune_mask = [tuning_param.tuning_flag];

% get feed/load port indices
feed_tune_ports  = port_ids(feed_mask | tune_mask);


% -------------------------------------------------------------------------
%             V assembly 
% -------------------------------------------------------------------------

tic

for pnum = 1:length(feed_tune_ports)
    
    % get current excitation port
    port_num = feed_tune_ports(pnum);
    
    % assembly excitation 
    [Vp] = excitation_coil(B_o,index,node,edge,port(port_num));
    F(:,pnum) = Vp;
end