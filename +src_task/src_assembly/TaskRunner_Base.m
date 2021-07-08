classdef (Abstract) TaskRunner_Base < handle
% declare an abstract Task class 
    
% ====================================================================== %
    % private properties
    properties (SetAccess = protected)
        
        
        % coil loader
        coil_loader
        
        task_settings_
        
        task_solution_ 
        
        % temporary var. for resulting electric fields
        Etot
        Etot_spec
        
        % temporary var. for resulting magnetic fields
        Htot
        Htot_spec
        B1_plus
        B1_minus
        
    end
    
% ====================================================================== %
    
    properties (Constant)
        % sizer object keeps track of system sizes
        dims = src_scatterer.sizer();
    end
    
% ====================================================================== %    

    methods
        % define constructor
        function obj = TaskRunner_Base(task_settings)
            obj.task_settings_ = Task_settings(task_settings);
            obj.task_solution_ = Task_Solution();
            
            % update basis dimensions
            obj.update_dims_();
        end
        
        % ------------------------------------------------------------- %
        
        function sol = get_solution(obj)
            % method returns solution
            sol = obj.task_solution_;
        end
                
        
    end
    
% ====================================================================== %
    
    methods (Abstract)
        
        % run task
        run(obj)
        
    end
    
% ====================================================================== %
    
    methods (Access = protected)
        
        function [freqs,N_freq] = get_working_frequencies_(obj)
            % get frequencies from general settings
            [freqs, N_freq]  = obj.task_settings_.get_working_frequencies();
            
            % update dims
            obj.dims.N_freqs = N_freq;
        end
        
        
        function update_dims_(obj)
            % function updates dimensions (ql, l)
            pwx = obj.task_settings_.vie.PWX;
            
            % update basis type
            switch pwx
                case 0
                    obj.dims.l  = 1;
                case 1
                    obj.dims.l  = 4;
                otherwise
                    error("Unsupported/unknown basis order (VIE)");
            end
            
            % update resulting size
            obj.dims.ql = obj.dims.q * obj.dims.l;
            
        end
        
        % ------------------------------------------------------------------------
        
        function [] = compute_E_total_base_(obj, Jb_sol, mvp, i_freq)
            % computes the electric field
 
            % iterate over ports
            for feed_port = 1:obj.dims.N_feeds
                
                % get polarisation currents
                Jb = Jb_sol(1:end,feed_port, i_freq);
                
                t_Etotal_spec = mvp.Jb2Eb(Jb);
                t_Etotal = zeros(size(t_Etotal_spec));
                
                obj.Etot_spec(:,:,:,:,feed_port, i_freq) = t_Etotal_spec;
                
                % PWL scaling
                if obj.dims.ql == 12
                    
                    G = repmat([1; 12; 12; 12], obj.dims.q,1);
                    
                    for i = 1:obj.dims.ql
                        t_Etotal(:,:,:,i) = G(i) * t_Etotal_spec(:,:,:,i);
                    end
                else
                    t_Etotal = t_Etotal_spec;
                end
                
                % store only the total electric field
                obj.Etot(:,:,:,:,feed_port,i_freq) = t_Etotal;
                
            end
            
        end
        %------------------------------------------------------------------------
        function [] = compute_Htot_pfft_(obj, Jc, Jb, mvp,i_freq)
            
            res = obj.scatterer.dom_vie.res;

            % iterate over ports
            for feed_port = 1:obj.dims.N_feeds
                   
                idxS    = obj.assembly_vsie.scatterer.index_vie.index_ql(obj.dims.ql, obj.dims.Nvox_vie);
                H_tot   = zeros(obj.dims.vie);
                H_tot(idxS) = mvp.Jcb2Htot(Jc(:,feed_port), Jb(:,feed_port));
                
                % store only the total electric field
                obj.Htot(:,:,:,:,feed_port,i_freq) = H_tot / res^3;
                
            end
            
        end
       %------------------------------------------------------------------------
       function [] = compute_Htot_basis_(obj, Jc, Jb, mvp, i_freq)
            
             % get dimensions of the problem
            dims = obj.dims;
            n1 = dims.L_vie;
            n2 = dims.M_vie;
            n3 = dims.N_vie;
            n_feeds = dims.N_feeds;
            
            res = obj.scatterer.dom_vie.res;
            
            % compute tested total electric field
            t_Htotal = zeros(n1,n2,n3,dims.ql);
            
            idxS = obj.assembly_vsie.scatterer.index_vie.index_ql(dims.ql, dims.Nvox_vie);
            gpu_flag = obj.task_settings_.general.GPU_flag;
                        
            Jsol = zeros(n1 * n2 * n3 * dims.ql, n_feeds);
            Jsol(idxS,:) = Jb;
            Jsol = reshape(Jsol,[n1, n2, n3, dims.ql, n_feeds]);
                        
            for feed_port = 1:dims.N_feeds
                
                % compute tested incident magnetic field
                t_Hinc = mvp.get_Hinc_coupling(Jc(:,feed_port));
                
                % compute tested scattered magnetic field
                t_Hsca = mvp.K_vie(Jsol(:,:,:,:,feed_port));
                
                % compute total tested magnetic field
                if gpu_flag
                    t_Htotal(idxS) = gather(t_Hsca(idxS) + t_Hinc);
                else
                    t_Htotal(idxS) = t_Hsca(idxS) + t_Hinc;
                end
                
                % PWL scaling
                if dims.ql == 12
                    
                    G = repmat([1; 1/12; 1/12; 1/12;], obj.dims.q,1);

                    for i = dims.ql
                        t_Htotal(:,:,:,i) = G(i).*t_Htotal(:,:,:,i);
                    end
                    
                end
                

                % store only the total magnetic field
                obj.Htot(:,:,:,:,feed_port,i_freq) = t_Htotal / res^3;
            end
       end
             
    end
    
end