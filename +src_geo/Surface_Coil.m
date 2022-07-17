classdef Surface_Coil < handle
% Class for surface RF coils
% Author: Georgy Guryev, Cambridge, MA, 2019    

    properties (SetAccess = immutable)
        
        % geometric parameters
        index = []
        etod  = []
        edge  = []                         
        elem  = []
        index_elem = []

        
        % Number of unknowns/ports/feed elements
        Nc      = 0;
        N_feed  = 0;
        N_tuning = 0;
 
    end
    
    properties 
        
       % port structure and nodes might be updated (change of
       % loads/matching circuits for port; rotation or displacement for node)
        port  = [];
        name  = [];
        node  = [];                         
        Ct = []
        Ln = []
        Pn = []
        rp = []
        rn = []
        r2 = []
        r3 = []
    end
    
    methods
        
        function obj = Surface_Coil(COIL_struct)
            % constructor for class Surface coil
            
            obj.index = COIL_struct.index;
            obj.etod  = COIL_struct.etod;
            obj.node  = COIL_struct.node;
            obj.edge  = COIL_struct.edge;
            obj.elem  = COIL_struct.elem;
            obj.port  = COIL_struct.port;
            obj.index_elem = COIL_struct.index_elem;
            obj.Ct = COIL_struct.Ct;
            obj.Ln = COIL_struct.Ln;
            obj.Pn = COIL_struct.Pn;
            obj.rp = COIL_struct.rp;
            obj.rn = COIL_struct.rn;
            obj.r2 = COIL_struct.r2;
            obj.r3 = COIL_struct.r3;
            
            % update relevant problem size information
            obj.Nc     = max(obj.index);
            obj.N_feed = nnz(strcmp({obj.port.type},'feed'));
        end
        
%-------------------------------------------------------------------------%        
        function obj = update_coil_position(coil)
            % update_coil_position()
            
            % update coil coordinates
            obj.node = coil.node;
            
            % update port parameters
            obj.port = coil.port;
            
            % update Ct, Pn
            obj.Ct = coil.Ct;
            obj.Ln = coil.Ln;
            obj.Pn = coil.Pn;
            
        end
               
    end
end