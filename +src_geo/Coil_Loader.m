classdef Coil_Loader < handle
% Class Coil Loader provides an interface to parse geo files and  
% Author: Georgy Guryev, Cambridge, MA, 2019    
    
    properties (Constant)
        
        % invalid extension(s) error
        ext_err_msg = 'Coil file(s) has/have wrong extension';
        
        % invalid port_specs file format
        specs_err_msg = 'Wrong port specs file format';
        
        % invalid coil msh file format
        msh_err_msg = 'Wrong coil mesh file format';
        
    end
% ====================================================================== %    
    
    properties (Access = private)
        
        % file with port specifications
        port_specs_file_ = '';
        
        % file with coil mesh
        coil_msh_file_   = '';
        
        % coil path
        coil_path_       = '';
        
        % loader status
        loader_status_ = false;
                
    end
    
% ====================================================================== %    
    
    methods
        
        
        function obj = Coil_Loader(sie_settings)
            % Coil Loader constructor

            % setting properties
            obj.port_specs_file_ = sie_settings.CoilPortsFile;
            obj.coil_msh_file_   = sie_settings.CoilMeshFile;
            obj.coil_path_       = sie_settings.Coil_Path;
            
            % validate files' extensions
            obj.validate_file_extensions_()
            assert(obj.loader_status_, obj.ext_err_msg);
        end
                
        % ------------------------------------------------------------- %

        function COIL = load_coil(obj)
            % method loads ports_specs and mesh file
            
            % read port specs
            ports = load_port_specs_(obj);
            
            % read coil mesh geometry; arange ports according to port_specs
            COIL_struct = load_coil_mesh_(obj, ports);
            
            % instantiate Surface_Coil with a COIL structure
            COIL = src_geo.Surface_Coil(COIL_struct);
        end
        
        % ------------------------------------------------------------- %

        
    end
    
% ====================================================================== %    
    
    methods (Access = private)
  
        % ------------------------------------------------------------- %

        function ports = load_port_specs_(obj)
            % private method reads port specs from file
              ports = src_geo.parse_PortSpecs(obj.port_specs_file_);
        end
        
        % ------------------------------------------------------------- %
        
        
        function COIL_struct = load_coil_mesh_(obj, ports)
            % private method reads mesh from .msh file
                        
            % parse .msh file; extract elementary entities 
            [node,~,e,elem,~] = Mesh_Parse(fullfile(obj.coil_path_,...
                                           obj.coil_msh_file_));
                         

            % extract higher-level geometric entities
            [edge,etod,index,ports,index_elem] = Mesh_PreProc(e,elem,ports);
            
            
            % validate coil mesh orientation         
            obj.validate_coil_mesh_(index_elem, etod);
            
            % "validate" and "fix" orientation of triangles (XOR with validate_coil_mesh_()!!!)
%             [elem_new] = obj.validate_coil_mesh_old_(node, elem);
            
            % 
            [Ct,Ln,Pn] = Mesh_CLP(node,elem);
            
            
            % save COIL data to structure
            COIL_struct = struct('index',  index,...
                                 'etod', etod, ...
                                 'node', node, ...
                                 'edge', edge, ...
                                 'elem', elem, ...
                                 'index_elem', index_elem, ...
                                 'Ct', Ct, ...
                                 'Ln', Ln, ...
                                 'Pn', Pn, ...
                                 'port', ports);
                             
            % form lists for vertex locations
            COIL_struct = obj.get_vertex_locations_(COIL_struct);
        end
        
        % ------------------------------------------------------------- %
        
        function [] = validate_file_extensions_(obj)
            % private function checks if files have valid extension
            
            % get files' extensions
            specs_file_ext = obj.port_specs_file_(end-2:end);
            msh_file_ext   = obj.coil_msh_file_(end-2:end); 
 
            % validate extensions
            if(strcmp(specs_file_ext,'txt') && strcmp(msh_file_ext,'msh'))
                obj.loader_status_ = true;
            else
                obj.loader_status_ = false;
            end
     
        end
        
        % ------------------------------------------------------------- %
        
        function validate_coil_mesh_(obj, index_elem, etod)
        % verify consistancy of orientation of mesh triangles
            
            ie_EA = index_elem.EA(:,1);
            je_EA = index_elem.EA(:,2);
            n_EA_elem = length(ie_EA);

            for index_EA = 1:n_EA_elem
                ie = ie_EA(index_EA);
                je = je_EA(index_EA);

                if sum(ismember(etod(:,ie),etod(:,je)))
                    error('There are both clock-wise and counter-clock-wise triangles in the mesh. Shift the incosistnet line loops in the .geo file, by visualizing the normals in GMSH.');
                end

            end
        end
        
        % ------------------------------------------------------------- % 
        function [elem] = validate_coil_mesh_old_(obj, node, elem)
        % validate_coil_mesh_old_ "fixes" inconsistances associated with
        % wrong mesh orientation
        
            mesh_center = sum(node,2);
            new_elem    = elem;
            number_flipped = 0;
            
            for k=1:size(elem,2)
                
                % get radius vectors for vertices
                r_k = reshape(node(:,elem(1:3,k)),3,3);
                
                % compute a vector, orthogonal to current triangular patch
                normal = cross(r_k(:,2) - r_k(:,1), r_k(:,3) - r_k(:,2));
                
                % get radius vector pointing from face inwards (ie towards mesh center)
                to_center  = mesh_center - sum(r_k,2);
                
                % project normal on radius-vector facing inwards
                projection = dot(normal, to_center);
                
                % if orientation of the normal is different - FIX it
                if projection < 0
                    new_elem(1:3,k) = elem(3:-1:1,k);
                    number_flipped  = number_flipped + 1;
                end
            end
            
            elem = new_elem;
        
        end
        
       % ---------------------------------------------------------------- % 
        
        function [coil] = get_vertex_locations_(obj, coil)
            
            N_sie = max(coil.index);
            r_v = zeros(12, N_sie);
            
            for i = 1:N_sie
                % get current rwg vertices
                r_v(:,i) = src_coupling.get_rwg_vertices(coil, i);
            end         
            r_v = r_v.';
            coil.rp = r_v(:,1:3).';
            coil.rn = r_v(:,4:6).';
            coil.r2 = r_v(:,7:9).';
            coil.r3 = r_v(:,10:12).';
        end
                
    end
    
  % ====================================================================== %    
  
end