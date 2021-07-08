classdef Coupling_Basis_Cross < Coupling_Basis_Base  
    
    properties
        
        cross
        
        % store left singular vectors for N and K operators
        V_N
        V_K
    end
    
    
    methods
        
        function obj = Coupling_Basis_Cross(task_settings, scatterer, coil, operator, dims, Zc_inv, freq)
            
            obj = obj@Coupling_Basis_Base(task_settings, scatterer, coil, operator, dims, freq);
            
            obj.Zc_inv = Zc_inv;
            
        end
        
        % --------------------------------------------------------------- %
        
        function  construct_basis(obj)
            
            % instantiate cross approximation 
            obj.cross = Cross_approx(obj.task_settings, obj.scatterer, obj.coil, obj.dims, obj.freq);
            
            % compute N and K crosses
            [obj.operator.U_N, obj.operator.alpha_N] = obj.cross.compute_approximation('N');
            [obj.operator.U_K, obj.operator.alpha_K] = obj.cross.compute_approximation('K');
            
            obj.operator.X_N = obj.operator.alpha_N * obj.Zc_inv;
            obj.operator.W_N = obj.operator.X_N * obj.operator.alpha_N.';
            
            % save basis
            obj.save_basis();
        end
                
    end
   
end