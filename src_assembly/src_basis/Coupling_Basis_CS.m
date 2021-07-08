classdef Coupling_Basis_CS < Coupling_Basis_Base  
    
    properties
        
        column_sampling

    end
    
    
    methods
        
        function obj = Coupling_Basis_CS(task_settings, scatterer, coil, operator, dims, Zc_inv, freq)
            
            obj = obj@Coupling_Basis_Base(task_settings, scatterer, coil, operator, dims, freq);
            
        end
        
        % --------------------------------------------------------------- %
        
        function  construct_basis(obj)
            
            % instantiate cross approximation 
            obj.column_sampling = Column_sampling_approx(obj.task_settings, obj.scatterer, obj.coil, obj.dims, obj.freq);
            
            % compute N and K crosses
            [obj.U_N, obj.U_K] = obj.column_sampling.compute_approximation('NK');
            
            obj.operator.X_N = obj.operator.alpha_N * obj.Zc_inv;
            obj.operator.W_N = obj.operator.X_N * obj.operator.alpha_N.';
            
            % save basis
            obj.save_basis();
        end
                
    end
   
end