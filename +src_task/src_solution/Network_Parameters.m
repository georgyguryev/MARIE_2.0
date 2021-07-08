classdef Network_Parameters
    
    properties
      
        coil 
      
        src 
        
        tuning
        
        dims
    end
    
    
    methods
        
        function obj = Network_Parameters(dims)
            
            obj.dims = dims;
            
            N_feeds = obj.dims.N_feeds;
            N_freqs = obj.dims.N_freqs;
            
            obj.coil   = Elementary_Network_Param(N_feeds, N_freqs);
            obj.src    = Elementary_Network_Param(N_feeds, N_freqs); 
            obj.tuning = Elementary_Network_Param(N_feeds, N_freqs);
        end
        
    end
    
end