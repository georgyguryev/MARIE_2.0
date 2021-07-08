classdef EM_operators < handle
    % class EM_operators is a container of EM operators, that had to be
    % assenbled
    
    properties
        % N operators for vie, extended (pFFT) and near (pFFT) domains
        N_vie_
        N_ext_
        N_near_
        N_vie_cell_
        N_ext_cell_
        N_near_cell_
        
        N_cell_idx
        
        
        % GPU versions of N operators 
        N_vie_gpu
        N_ext_gpu
        N_near_gpu
        
        
        % K operators for vie, extended (pFFT) and near (pFFT) domains
        K_vie_
        K_ext_
        K_near_
        
        K_vie_gpu
        K_ext_gpu
        K_near_gpu
        
        P_near_
        
        Jc_ini
        
        U_N
        U_K
        
        UU_N
        V_N
        X_N
        X_K
        W_N
        W_K
        alpha_N
        alpha_K
        alpha_Zc_inv
        M_q2ql
        domain
        index
        
        Zbc
        Zbc_Nop
        Zbc_Kop
        
        dims;
    end
    
    methods 
 %----------------------------------------------------------------------%
        function obj = EM_operators(dims)
            obj.dims = dims;
        end
        
 %----------------------------------------------------------------------%
       
        
        function N_to_cell(obj, dom_str, gpu_flag)
            
            obj.N_cell_idx   = obj.expansion_indices_();

            switch dom_str
                case 'vie'
                    obj.N_vie_cell_ = obj.tensor_to_cell_tensor_(obj.N_vie_, gpu_flag);
                case 'ext'
                    obj.N_ext_cell_ = obj.tensor_to_cell_tensor_(obj.N_ext_, gpu_flag);
                case 'near'
                    obj.N_near_cell_ = obj.tensor_to_cell_tensor_(obj.N_near_, gpu_flag);
            end
        end      
        
        
    end
    
 %----------------------------------------------------------------------%
    
    methods(Access=private)
        
        function cell_op = tensor_to_cell_tensor_(obj,operator, gpu_flag)
            
            
            idx_ql = obj.N_cell_idx.idx_lq_vec;

            ql = size(idx_ql,1); 
            
            cell_op = cell(1,ql);
           
       
            if(gpu_flag)
                for i = 1:ql
                    cell_op{i} = gpuArray(operator(:,:,:,idx_ql(:,i)));
                end
            else
                for i = 1:ql
                    cell_op{i} = operator(:,:,:,idx_ql(:,i));
                end
            end
        end
        
        %-----------------------------------------------------%
        
        function [idx] = expansion_indices_(obj)
            
            
            if 12 == obj.dims.ql
                
                map_pq   = [1 2 3; 2 4 5; 3 5 6];
                map_lplq = [1 5 6 7; -5 2 8 9; -6 8 3 10; -7 9 10 4];
                
                % generate indices to unfold loops
                [p_msh, lp_msh] = meshgrid(1:obj.dims.q, 1:obj.dims.l);
                
                p_vec  = kron(p_msh(:), ones(obj.dims.ql,1));
                lp_vec = kron(lp_msh(:), ones(obj.dims.ql,1));
                q_vec  = repmat(p_msh(:), obj.dims.ql,1);
                lq_vec = repmat(lp_msh(:), obj.dims.ql,1);
                
                
                pq_vec     = sub2ind([obj.dims.q, obj.dims.q], p_vec, q_vec);
                lplq_vec   = sub2ind([obj.dims.l, obj.dims.l], lp_vec, lq_vec);
                idx.left_vec   = sub2ind([obj.dims.l, obj.dims.q], lp_msh, p_msh);
                
                % define final indexing for N_mvp_pwl
                idx_pq_vec       = reshape(map_pq(pq_vec), [obj.dims.ql, obj.dims.ql]);
                idx_lplq_vec     = reshape(map_lplq(lplq_vec), [obj.dims.ql, obj.dims.ql]);
                
                idx.idx_lq_vec   = reshape(sub2ind([10, 6], abs(idx_lplq_vec), idx_pq_vec),...
                                           obj.dims.ql, obj.dims.ql);
                idx.sgn_vec      = sign(idx_lplq_vec);
            else
                
                idx.idx_lq_vec = [1 2 3; 2 4 5; 3 5 6];
                idx.left_vec   = [1 2 3];
            end
        end
        
    end
    
end