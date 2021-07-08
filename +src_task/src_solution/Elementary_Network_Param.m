classdef Elementary_Network_Param
    
    properties
        
        Y
        Z
        S
        
    end
    
    methods
        
        function obj = Elementary_Network_Param(N_feeds, N_freqs)
            
            obj.Y = zeros(N_feeds, N_feeds, N_freqs);
            obj.Z = zeros(N_feeds, N_feeds, N_freqs);
            obj.S = zeros(N_feeds, N_feeds, N_freqs);
            
        end
    
    end
    
end