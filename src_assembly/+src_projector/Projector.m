classdef Projector < handle
   
    properties
        
        
        central_cell_id
        
        cell_id
        
        near_vox_list
                
        C2B_near_list       % list of near scatterer voxels from near list
        
        C2C_near_list
                
        index_exp2near
        
        collocation_mode
        
        near_boundary_width
        
        PS
        
    end
    
    methods
        
        function form_near_lists(obj, coil, scatterer, projector, dims)
            
            % call method form_near_list(); update projector properties
            [obj.central_cell_id, obj.cell_id, obj.near_vox_list,...
             obj.C2C_near_list, obj.C2B_near_list, obj.index_exp2near] = src_projector.form_near_lists(coil, scatterer, dims);
        
        end
        
        % --------------------------------------------------------------- %
        
        function assemble_projection_matrix(obj, setup, freq)
            
            % get collocation mode
            projector_assembly_mode = setup.task_settings_.vsie.Proj_assembly_mode;
            
            
            % select appropriate method for assembly
            switch projector_assembly_mode
                
                % point collocation
                case 'Sphere'
                    obj.PS = src_projector.projection_assembly_sphere(setup.coil, setup.scatterer, setup.projector,...
                                                                 setup.dims, setup.task_settings_, freq);
                % voxel quadrature collocation
                case 'Cube'
                    
                    % assemble projection_mvp
                    fP_mvp_near = setup.assemble_near_mvp();
                    obj.PS = src_projector.projection_assembly_vox(fP_mvp_near, setup.coil, setup.scatterer,setup.projector,...
                                                              setup.dims, setup.task_settings_, freq);
            end
           
        % End of method assemble_projection_matrix()    
        end
        
    end
    
end