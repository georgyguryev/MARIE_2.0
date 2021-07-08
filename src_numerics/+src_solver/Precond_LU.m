classdef Precond_LU < handle

% ====================================================================== %

    properties (Access = private)
        
        L_sie
        U_sie
        P_sie
        
        L_vie
        U_vie
        P_vie
        
        L_vsie
        U_vsie
        P_vsie
        
        basis_vsie
        
        dims
        
    end
   
% ====================================================================== %
    
    methods 
        
        % --------------------------------------------------------------- %
        
        function obj = Precond_LU(dims)
            % constructor of class Precond_LU
            obj.dims = dims;
        end
        
        % --------------------------------------------------------------- %
        
        function assemble_precond_sie(obj, Zc)
            % method assemble_precond_sie()
            
            % compute lu decomposition with permutation vector pc
            [obj.L_sie, obj.U_sie, obj.P_sie] = lu(Zc, 'vector');
            
            % get scaling factor from Uc
            D_sqrt = diag(sqrt(diag(obj.U_sie)));
            
            % symmetrize Lc and Uc for better convergence of GMRES
            obj.U_sie = D_sqrt \ obj.U_sie;
            obj.L_sie = obj.L_sie * D_sqrt;
            obj.P_sie = obj.P_sie;
            
        end 
        
        % --------------------------------------------------------------- %
        
        function assemble_precond_vie(obj, idx_1d, Mcr_inv, mvp, freq)
            
            emu = src_utils.EM_utils(freq);
            
            % define block-diagonal preconditioner of size [N_scat * l]
            M_vie = speye(length(idx_1d) * obj.dims.l);
            
            % 
            M_vie = spdiags(sqrt((1 / emu.ce) .* mvp.G_prec(Mcr_inv(idx_1d))), 0 , M_vie);
            
            % store results
            obj.P_vie = 1:obj.dims.q * size(M_vie);
            obj.L_vie = sparse(blkdiag(M_vie, M_vie, M_vie));
            obj.U_vie = obj.L_vie;
            
%             % store results
%             obj.P_vie = []; 
%             obj.L_vie = full(M_vie(1,1)); 
%             obj.U_vie = full(M_vie(1,1)); 
            
            
        end
        
        % --------------------------------------------------------------- %
        
        function assemble_precond_vsie_coupled(obj, idx_1d, Mcr_inv, mvp, freq)
            
            emu = src_utils.EM_utils(freq);

            % replicate material properties l times
            Mb = 1 / emu.ce .* mvp.G_prec(Mcr_inv(idx_1d));
            
            % form diagonal block with material properties
            Mb_prec = spdiags(sqrt(Mb), 0, speye(length(Mb)));
            
            % store results
            obj.P_vsie = [obj.P_sie, max(obj.P_sie) + (1:size(Mb_prec))];
            obj.L_vsie = sparse(blkdiag(obj.L_sie, Mb_prec));
            obj.U_vsie = sparse(blkdiag(obj.U_sie, Mb_prec));
            
        end
        
        % --------------------------------------------------------------- %
       
        function assemble_precond_vsie_decoupled(obj, idx_1d, Mcr_inv, mvp, freq)
       
            emu = src_utils.EM_utils(freq);
                        
            % replicate material properties l times
            Mb = (emu.ce ./ mvp.G_prec(Mcr_inv(idx_1d)));
            
            % form diagonal block with material properties
            Mb_prec = spdiags(Mb(:), 0, speye(length(Mb(:))));
            
            obj.basis_vsie = sparse(blkdiag(eye(obj.dims.N_sie), Mb_prec));
            
            
        end
        
        
        % --------------------------------------------------------------- %
        
        function precond = sie(obj)
            % function returns sie preconditioner int precond struct
            
            precond.L = obj.L_sie;
            precond.U = obj.U_sie;
            precond.P = obj.P_sie;
        end
        
        % ---------------------------------------------------------------
        
        function precond = vie(obj)
            % function returns sie preconditioner int precond struct
            
            precond.L = obj.L_vie;
            precond.U = obj.U_vie;
            precond.P = obj.P_vie;
        end
        
        % --------------------------------------------------------------- %
        
        function precond = vsie(obj)
            % function returns sie preconditioner int precond struct
            
            precond.L = obj.L_vsie;
            precond.U = obj.U_vsie;
            precond.P = obj.P_vsie;
        end
        % --------------------------------------------------------------- %
        
        function precond = vsie_basis(obj)
            
            Nb = obj.dims.N_scat * obj.dims.ql + obj.dims.N_sie;
            
            % function returns sie preconditioner int precond struct
            precond.L = obj.basis_vsie;
            precond.U = speye(Nb);
            precond.P = 1:Nb;
        end

    end
% ====================================================================== %
   
end