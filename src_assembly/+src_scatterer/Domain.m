classdef Domain < handle 
    
    properties (SetAccess = immutable)
        
        r
        x
        y
        z
        x_tensor
        y_tensor
        z_tensor
        res
    end
    
    
    methods
          
        function obj = Domain(r)
            % constructor Domain(r) 
            
            obj.r = r;
            
            % get 3D domain coordinates
            obj.x_tensor = squeeze(r(:,:,:,1));
            obj.y_tensor = squeeze(r(:,:,:,2));
            obj.z_tensor = squeeze(r(:,:,:,3));
            
            % get 1D grid vector
            obj.x = squeeze(r(:,1,1,1));
            obj.y = squeeze(r(1,:,1,2));
            obj.z = squeeze(r(1,1,:,3));
            
            % get resolution
            obj.res = abs(obj.x(2) - obj.x(1));
            
        end
    end
    
    
end