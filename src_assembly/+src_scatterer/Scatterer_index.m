classdef Scatterer_index < handle
    
    properties (SetAccess = immutable) 
        S_1d
        S_3d
        S_1d_vie2ext
    end
    
    methods 
        function obj = Scatterer_index(index_1d, Nvox, S_1d_vie)
            % constructor of class Scatterer_index keeps track of indices
            
            % store 1D and 3D scatterer indicies
            obj.S_1d = index_1d;
            obj.S_3d = [index_1d; index_1d + Nvox; index_1d + 2 * Nvox];
            
            % if the scatterer index is generated form a dictionary to map
            % 1D indices from extended to vie domain
            if nargin > 2
                obj.S_1d_vie2ext = containers.Map(obj.S_1d, S_1d_vie);
            end
            
        end
        
        
        function idx_ql = index_ql(obj, lq, dom_dims)
            
            idx_ql = repmat(obj.S_1d, lq,1) + kron(dom_dims * [0:lq-1].', ones(size(obj.S_1d)));
        end
    end
end