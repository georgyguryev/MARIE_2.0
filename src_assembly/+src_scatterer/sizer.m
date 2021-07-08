classdef sizer < handle
    
    properties 
        
        % dimensions of VIE
        L_vie  
        M_vie
        N_vie
        vie 
        Nvox_vie
        
        % dimensions of exteded domain
        L_ext
        M_ext
        N_ext
        ext
        Nvox_ext
        
        % expansion basis dimensions 
        N_exp_1D
        N_exp_3D
        exp
        
        % dimensions of the near zone
        N_near_3D
        N_near_1D
        near
        
        N_scat
        
        % dimensions SIE 
        N_t
        N_i
        N_sie
        N_ports
        N_feeds
        N_loads
        N_elems                 % number of triangles 
        
        % operator dimensions
        L_Nop_vie
        M_Nop_vie
        N_Nop_vie
        op_vie
        
        % operator dimensions
        L_Nop_ext
        M_Nop_ext
        N_Nop_ext
        op_ext
        
        % operator dimensions
        L_Nop_near
        M_Nop_near
        N_Nop_near
        op_near
        
        % dimensions 
        N_freqs
        
        % basis type
        q = 3
        l 
        ql
        
        % 
        N_scat_vox  
        
        
        
    end
    
    methods
        
        function obj = sizer()
            % empty default constructor of class sizer
                        
        end
        
        function get_coil_dims(obj, coil)
            % function gets coil dimensions 
            
            obj.N_sie   = max(coil.index);
            obj.N_elems = size(coil.elem,2);
            obj.N_ports = size(coil.port,1);
            obj.N_feeds = nnz(strcmp({coil.port(:).type}, 'feed'));
            obj.N_loads = obj.N_ports - obj.N_feeds;
        end
        
        function get_body_dims(obj, r)
            
            % call private method to get domain dimensions
            dom_vie = obj.update_domain_dims_(r);
            
            % update domain dimensions
            obj.vie   = dom_vie.dims;
            obj.L_vie = dom_vie.L;
            obj.M_vie = dom_vie.M;
            obj.N_vie = dom_vie.N;
            obj.Nvox_vie = dom_vie.Nvox;
            
        end
        
        function get_extended_body_dims(obj, r)
            % function updates dimensions for extended domain
            
            % call private method to get domain dimensions
            dom_ext = obj.update_domain_dims_(r);
            
            % update domain dimensions
            obj.ext   = dom_ext.dims;
            obj.L_ext = dom_ext.L;
            obj.M_ext = dom_ext.M;
            obj.N_ext = dom_ext.N;
            obj.Nvox_ext = dom_ext.Nvox;
        end
        
        function set_nonair_voxel_dims(obj, index_1d)
            % function sets up the size of non-air voxels
            
            obj.N_scat_vox = size(index_1d,1);
        end
        
    end
    
    
    methods (Access = private)
        
        function dom = update_domain_dims_(obj, r)
            % read r dimensions
            dom.dims = size(r);
            
            % update obj.dims components
            dom.L = dom.dims(1);
            dom.M = dom.dims(2);
            dom.N = dom.dims(3);
            dom.Nvox = dom.L * dom.M * dom.N;
            
            % update obj.dims.vie to match currents' dimensions
            dom.dims(4) = obj.ql;

        end
        
    end
    
end