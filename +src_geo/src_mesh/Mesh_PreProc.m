function [edge,etod,index,port,index_elem] = Mesh_PreProc(e,elem,port)
% _________________________________________________________________________
%
%
%   Pre-process a mesh file from gmsh into SIE-friendly data
%
% _________________________________________________________________________
%
%
%
% -------------------------------------------------------------------------
%
%   J. Fernandez Villena -- jvillena@mit.edu
%   A.G. Polimeridis -- thanos_p@mit.edu
%   Computational Prototyping Group, RLE at MIT
%
%   Modified in 2019 by Georgy Guryev -- georgy@mit.edu 
%   Computational Prototyping Group, RLE at MIT
% _________________________________________________________________________
%

% Nn  = size(node,2); % Number of nodes
Nbd     = size(e,2);    % Number of exterior (boundary) edges
Ne      = size(elem,2); % Number of elements
NS_mat  = zeros(Ne, Ne);

v_counter_NS = zeros(Ne, 1);
% -------------------------------------------------------------------------
%   Pre-processing: generate edge and etod: edge[2xNd], etod[3xNe]
%   Finds adjancencies
%   Assigns boundaries
% -------------------------------------------------------------------------

% t1 = clock;


first_node  = [3 1 2]; % define edge vectors clockwise
second_node = [2 3 1];

% for edge and etod generation
edge = zeros(2,3*Ne); % allocate maximun possible size
etod = zeros(3,Ne); % allocate size

% define array for boundary
% % kn=-1*ones(3*Ne,1); % allocate maximun possible size
kn=zeros(3*Ne,1); % allocate maximun possible size

% define arrays for adjacency
eparent = zeros(3*Ne,1); % allocate maximum possible size for adjacency check

NS = zeros(floor(Ne*Ne/2),2); nscount = 0;
VA = zeros(100*Ne,2); vacount = 0;
EA = zeros(3*Ne,2); eacount = 0;

% counter of terminal edges
term_count = 0;

% loop on elements
Nd = 0;
for ii=1:Ne

    count = Nd;
    flagpad = zeros(ii-1,1);
    flagead = zeros(ii-1,1);


    for in=1:3

        node1 = elem(first_node(in),ii);
        node2 = elem(second_node(in),ii);
        flage  = 0;

        % check vertex adjacency
        R = elem(1,1:ii-1) - node1; flagpad(R == 0) = 1;
        R = elem(2,1:ii-1) - node1; flagpad(R == 0) = 1;
        R = elem(3,1:ii-1) - node1; flagpad(R == 0) = 1;

        for jj=1:Nd

            if (node1 == edge(1,jj)) && (node2 == edge(2,jj))
                etod(in,ii) = jj; % positive match found
                flage = 1; % edge adjacency
                flagead(eparent(jj)) = 2;

                % update kn if kn(jj) was exterior edge
                if 0 == kn(jj)
                    kn(jj) = -1; % set boundary flag as internal node
                end
            end

            if (node2 == edge(1,jj)) &&  (node1 == edge(2,jj))
                etod(in,ii) = -jj; % negative match found
                flage = 1; % edge adjacency
                flagead(eparent(jj)) = 2;

                % update kn if kn(jj) was exterior edge
                if 0 == kn(jj)
                    kn(jj) = -1; % set boundary flag as internal node
                end
            end

        end


        if (flage == 0) % new edge

            count = count+1;
            edge(1,count) = node1;
            edge(2,count) = node2;
            etod(in,ii) = count;

            eparent(count) = ii; % store the element to which the edge belongs


            for kk=1:Nbd % assign boundary
                bn1=e(1,kk);
                bn2=e(2,kk);
                if( (bn1 == node1) && (bn2 == node2) )
                    kn(count)=e(3,kk); % zero if external, positive integer if port
                end
                if( (bn2 == node1) && (bn1 == node2))
                    kn(count)=e(3,kk); % zero if external, positive integer if port
                end
            end


        end

    end


    % %     % check adjacency

    for kk = 1:length(flagead)

        adval = max(flagead(kk),flagpad(kk));
        switch adval
            case 1
                vacount = vacount+1;
                VA(vacount,1) = ii;
                VA(vacount,2) = kk;
            case 2
                eacount = eacount+1;
                EA(eacount,1) = ii;
                EA(eacount,2) = kk;
            case 0
                nscount = nscount+1;
                
                idx_1 = min(ii, kk);
                idx_2 = max(ii, kk);
                NS(nscount,1) = idx_1;
                NS(nscount,2) = idx_2;
                v_counter_NS(idx_1) =  v_counter_NS(idx_1) + 1;
                NS_mat(idx_1, v_counter_NS(idx_1)) = idx_2;
        end

    end

    Nd = count;

end

edge = edge(:,1:Nd);
kn = kn(1:Nd);

% remove replicas in adjacency
ST = [1:Ne; 1:Ne].';
EA = EA(1:eacount,:);
VA = VA(1:vacount,:);
NS = NS(1:nscount,:);

% -------------------------------------------------------------------------
%   Indexing & Number of unknown dofs (Nff)
% -------------------------------------------------------------------------

% allocate space to index array
index = zeros(Nd,1);

% different types of edges, and sort them in ascending order
etype = sort(unique(kn));
idx   = find(etype > 0);
etype = etype(idx);

% number of port edges
Nports = length(etype);

% check if ports
if Nports ~= size(port,1)
    error(' Dimension mismatch: number of ports in .txt does not match Nports in .msh');
end

% form feed-load mask
tE_list_mask = strcmp({port.type}, 'feed');
tL_list_mask = strcmp({port.type}, 'load');


% perform split
tE_list = etype(tE_list_mask);
tL_list = etype(tL_list_mask);

%% loop over excitation ports first

for i = 1:length(tE_list)
    
    % find positions of corresponding terminals
    idx = find(kn == tE_list(i));
    
    p_num = find(etype == tE_list(i));
    
    dofnum = term_count+1:term_count+length(idx); % create the dof number for each element
    
    index(idx) = dofnum; % assign the dof number to the index position
    
    port(p_num).t = dofnum; % store the dof number for the positive edge of port
    
    term_count = term_count+length(idx); % increase counter
    
end

%% loop over load ports

for i = 1:length(tL_list)
    
    % find positions of corresponding terminals
    idx = find(kn == tL_list(i));
    
    p_num = find(etype == tL_list(i));
    
    dofnum = term_count+1:term_count+length(idx); % create the dof number for each element
    
    index(idx) = dofnum; % assign the dof number to the index position
    
    port(p_num).t = dofnum; % store the dof number for the positive edge of port
    
    term_count = term_count+length(idx); % increase counter
    
end


%% Post-Process for mutual inductance

for port_num = 1:Nports
    if strcmp(port(port_num).tuning_param.LE,'mutual_inductor')
        
        % assign mutual terminals
        port(port_num).t_mutual = port(port(port_num).tuning_param.mutual_port).t;
    end
end


% now for the internal edges!
idx = find(kn == -1);
dofnum = term_count+1:term_count+length(idx); % create the dof number for each element
index(idx) = dofnum;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%               Find adjacency                                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

index_elem.ST = ST;
index_elem.VA = VA;
index_elem.EA = EA;
index_elem.NS = NS;
index_elem.NS_cell = mat2cell(NS_mat, ones(1,Ne), Ne);
