classdef  Scatterer_properties < handle
   
    properties (SetAccess = immutable)
        epsilon_r
        sigma_e
        rho
    end
    
    % ================================================================== %
    
    properties %(SetAccess = private)
        
        % resulting material properties matrices
        Mr
        Mc
        Mr_inv
        Mc_inv
        Mcr_inv 
        Mcr
        Mrc
        
    end
    
    % ================================================================== %

    methods
        
        function obj = Scatterer_properties(RHBM)
            % constructor of material properties
            
            % set initial material properties
            obj.epsilon_r = RHBM.epsilon_r;
            obj.sigma_e   = RHBM.sigma_e;
            obj.rho       = RHBM.rho;
        end
        
        % --------------------------------------------------------------- %

        function set_material_matrices(obj, Mr, Mc, Mr_inv, Mc_inv, Mcr_inv)
            % function sets Mr, Mc and Mc_inv, used in final mvp and
            % preconditioners
            
            % copy matrices with material properties
            obj.Mr = Mr;
            obj.Mc = Mc;
            obj.Mr_inv  = Mr_inv;
            obj.Mc_inv  = Mc_inv;
            obj.Mcr     = Mc .* Mr_inv;
            obj.Mrc     = Mr .* Mc_inv;
            obj.Mcr_inv = Mcr_inv;
        end
        
    end
    
end