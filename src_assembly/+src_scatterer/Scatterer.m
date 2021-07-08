classdef Scatterer < handle 
    
    properties (SetAccess = protected)
       
        % VIE/VSIE domains 
        dom_vie
        dom_ext
        dom_near
        
        % sizer object with all dimensions
        dims
        
        % material properties 
        prop_vie
        prop_ext
        
        % indices of of non-air voxels
        index_vie
        index_ext
        
    
        % coordinates of non-air voxels (VIE)
        Scoord
        
        
        % store frequency to infer scatterer voxels
        freq
    end
    
    methods 
        
        function obj = Scatterer(RHBM, dims, freq)
            % contructor class Scatterer()
            
            
            crop_RHBM    = obj.crop_air_voxels_(RHBM, freq);
            obj.dom_vie  = src_scatterer.Domain(crop_RHBM.r);
            obj.prop_vie = src_scatterer.Scatterer_properties(crop_RHBM);
            obj.dims     = dims;
            obj.freq     = freq;
            
            % get body dimensions           
            obj.get_body_dimensions_()
            
            % get nonair voxel indices
            index_vie_1d  = obj.get_nonair_voxels_idx_(obj.prop_vie);
            obj.index_vie = src_scatterer.Scatterer_index(index_vie_1d, obj.dims.Nvox_vie);
            
            % update scatterer non-air voxel dimensions
            obj.dims.set_nonair_voxel_dims(index_vie_1d);
            
            obj.dims.N_scat = size(index_vie_1d,1);
            
            % get coordinates of non-air voxels
            obj.set_scatterrer_coords_();
        end
        
        % ------------------------------------------------------------- %
        
        function set_extended_domain(obj, dom_ext, prop_ext)
            % function set_extended domain
            
            % set extended domain and material properties
            obj.dom_ext   = dom_ext;
            obj.prop_ext  = src_scatterer.Scatterer_properties(prop_ext);
            
            % get indices of non-air voxels
            index_ext_1d    = obj.get_nonair_voxels_idx_(prop_ext);
            
            % construct index structure with idx Mapping from ext2vie
            obj.index_ext = src_scatterer.Scatterer_index(index_ext_1d,...
                                                          obj.dims.Nvox_ext,...
                                                          obj.index_vie.S_1d);
            
            % update dims with a number of scatterer voxels
            obj.dims.N_scat = size(index_ext_1d,1);

        end
        
        % ------------------------------------------------------------- %
       
        function generate_near_domain(obj)
            
            % get dimensions of near domain
            N_near_1D = obj.dims.N_near_1D;
            
            % get coordinates of a 'near zone' domain
            r_near = obj.dom_ext.r(1:N_near_1D, 1:N_near_1D, 1:N_near_1D,:);
            
            % method set_near_domain(dom_near) sets the property dom_near
            obj.dom_near = src_scatterer.Domain(r_near);
        end
        
        
        
    end
    
    methods(Access = private)
        
        function index_1d = get_nonair_voxels_idx_(obj, properties)
            % method get_idx_4_nonair_voxels()
            
            % get em constants
            emu = src_utils.EM_utils(obj.freq);
            
            % get complex permittivities
            e_r = properties.epsilon_r + properties.sigma_e / (emu.ce);
            
            % find non-air voxels (check if e_r is below the threshold
            index_1d = find(abs(e_r(:) - 1) > 1e-12); 
            
        end
        
        % ------------------------------------------------------------- %
        
        function set_scatterrer_coords_(obj)
            obj.Scoord = [obj.dom_vie.x_tensor(obj.index_vie.S_1d), ...
                          obj.dom_vie.y_tensor(obj.index_vie.S_1d), ...
                          obj.dom_vie.z_tensor(obj.index_vie.S_1d)];
        end
        
        % ------------------------------------------------------------- %
        
        
        function get_body_dimensions_(obj)
        % private function updates obj.dims fields that correspond to body
        % dimensions
            
            % read r dimensions
            obj.dims.vie = size(obj.dom_vie.r);
            
            % update obj.dims components
            obj.dims.L_vie = obj.dims.vie(1);
            obj.dims.M_vie = obj.dims.vie(2);
            obj.dims.N_vie = obj.dims.vie(3);
            obj.dims.Nvox_vie = obj.dims.L_vie * obj.dims.M_vie * obj.dims.N_vie;

            
        end
        
        % ------------------------------------------------------------- %
        
        function crop_RHBM = crop_air_voxels_(~, RHBM, freq)
            % function crops RHBM to reduce dimensions of the FFT domain (N/K operators)
            
            emu = src_utils.EM_utils(freq);
            
            % get angular frequency 
            omega = 2 * pi * freq;
            
            % compute relative permittivity
            e_r = RHBM.epsilon_r - 1j * RHBM.sigma_e / (emu.e0 * omega);
            
            % find non-air voxels
            idx = find(abs(e_r - 1) > 1e-12);
            
            xd = RHBM.r(:,:,:,1);
            yd = RHBM.r(:,:,:,2);
            zd = RHBM.r(:,:,:,3);
                        
            % find non-air voxels' limits
            xmin = min(xd(idx)); xmax = max(xd(idx));
            ymin = min(yd(idx)); ymax = max(yd(idx));
            zmin = min(zd(idx)); zmax = max(zd(idx));
            
            % resize RHBM
            idxX_new = find(RHBM.r(:,1,1,1) <= (xmax + 1) & RHBM.r(:,1,1,1) >= (xmin - 1));
            idxY_new = find(RHBM.r(1,:,1,2) <= (ymax + 1) & RHBM.r(1,:,1,2) >= (ymin - 1));
            idxZ_new = find(RHBM.r(1,1,:,3) <= (zmax + 1) & RHBM.r(1,1,:,3) >= (zmin - 1));
            
            % return cropped data
            crop_RHBM.r         = RHBM.r(idxX_new,idxY_new,idxZ_new,:);
            crop_RHBM.rho       = RHBM.rho(idxX_new,idxY_new,idxZ_new);
            crop_RHBM.sigma_e   = RHBM.sigma_e(idxX_new,idxY_new,idxZ_new);
            crop_RHBM.epsilon_r = RHBM.epsilon_r(idxX_new,idxY_new,idxZ_new);
            
        end
        
    end
    
end