classdef TaskRunner_SIE < TaskRunner_Base
    
% ====================================================================== %
    
    properties 
        
        % assembly_sie
        assembly_sie_
        
        % surface coil
        coil
        
        % edge to edge impedance matrix
        Zcoil_
        Zcoil_inv_
        
        % port to terminal projection matrix
        Fcoil
        
        % coil edge currents 
        Jcoil
        
        Z_mt
        Z_L
        Z_L_hat
        Z_ref
        Z_sh
        T1 
        T2
        
        
        % port network parameters
        NP
        
        % feed and port indices
        feed_tune_ports_
          
    end
    
% ====================================================================== %
   
    methods
        
        % define Task_SIE constructor
        function [obj] = TaskRunner_SIE(task_settings)
            
            % call superclass constructor
            obj = obj@TaskRunner_Base(task_settings);
            
            obj.coil_loader = src_geo.Coil_Loader(task_settings.sie);
            
            % load COIL
            obj.load_coil_model_();
            
            % get coil dimensions
            obj.dims.get_coil_dims(obj.coil);
            
            obj.assembly_sie_ = src_assembly.Assembly_SIE(task_settings, obj.coil, obj.dims);
            
            % init solution
            obj.init_solution_();

        end
              
        % ------------------------------------------------------------- %
                
        % solve the problem
        function run(obj)
            
            % solve SIE task
            obj.run_frequency_sweep_();
                        
        end    

    end
    
% ====================================================================== %
   
    methods (Access = private)
        
        function [] = load_coil_model_(obj)
            % private method loads coil to task runner
            obj.coil = obj.coil_loader.load_coil();
        end
        
        % ------------------------------------------------------------- %
        
        % function for initialization of the solution
        function init_solution_(obj)
            
            % init task_solution object
            obj.task_solution_.init_SIE_solution(obj.dims);
            
            obj.Fcoil = [];
                                            
            % reset network params 
            obj.NP = Network_Parameters(obj.dims);
        end

        % ------------------------------------------------------------- %
        
        % assemble SIE system
        function  assemble_SIE_system_(obj, freq)
            
            % call some assebly method
            obj.assembly_sie_.assemble_system(freq);
         
            % copy assembled matrices
            obj.Zcoil_    = obj.assembly_sie_.Zcoil_;
            obj.Zcoil_inv_ = obj.assembly_sie_.Zcoil_inv_;
            obj.Fcoil    = obj.assembly_sie_.Fcoil_;
            obj.Z_mt      = obj.assembly_sie_.Z_mt;
            obj.Z_L_hat = obj.assembly_sie_.Z_L_hat;
            obj.Z_ref     = obj.assembly_sie_.Z_ref;
            obj.Z_sh      = obj.assembly_sie_.Z_sh;

         
            obj.feed_tune_ports_ = obj.assembly_sie_.feed_tune_ports_;
        end 
        
        % ------------------------------------------------------------- %
        
        function compute_surface_currents_(obj, freq_num)
            % function compute_surface_currents_() computes edge current
            % densities on the surface of the RF coil
            
            % if SIE runs at single frequency
            if nargin < 2 
                freq_num = 1;
            end
            
            w_src_I  = obj.assembly_sie_.w_src_I;
            w_coil_V = obj.assembly_sie_.w_coil_V;
            Zw_coil  = obj.assembly_sie_.Zw_coil;
            Zw_src   = obj.assembly_sie_.Zw_src;
            Yw_coil  = obj.assembly_sie_.Yw_coil;
             
            Jct = obj.assembly_sie_.Zcoil_ \ obj.Fcoil;
            
            
            obj.Jcoil = (Jct + Jct / ( inv(obj.Z_L_hat) - obj.Fcoil.' * Jct) * ...
                        obj.Fcoil.' * Jct) * obj.assembly_sie_.rhs_cp; 
                    
%             obj.Jcoil = obj.Zcoil_ \ obj.assembly_sie_.rhs_c;          
                                 
            % Compute coil currents and voltages 
            Icoil = - obj.Fcoil.' * obj.Jcoil;
            Vcoil = obj.assembly_sie_.rhs_cp + obj.Z_L_hat * obj.Fcoil.' * obj.Jcoil;
            
            % Compute currents and voltages seen from feed ports
            Isrc  = w_src_I * Icoil + Yw_coil * Vcoil;
            Vsrc  = w_coil_V * Vcoil + Zw_coil * Icoil + Zw_src * Isrc;
            
            % find NP seeing from the voltage sourse (i.e. coil + external circuit + ref. impedance)
            obj.NP.src.Y(:,:,freq_num) = Isrc / Vsrc;         
            obj.NP.src.Z(:,:,freq_num) = y2z(obj.NP.src.Y(:,:,freq_num));
            obj.NP.src.S(:,:,freq_num) = y2s(obj.NP.src.Y(:,:,freq_num));
            
            % find NP seeing from the coil ports
%             obj.NP.coil.Y(:,:,freq_num) = Icoil / Vcoil;
            obj.NP.coil.Y(:,:,freq_num) = - obj.Fcoil.' * Jct; 
            obj.NP.coil.Z(:,:,freq_num) = y2z(obj.NP.coil.Y(:,:,freq_num));
            obj.NP.coil.S(:,:,freq_num) = y2s(obj.NP.coil.Y(:,:,freq_num));
            
        end
        
        % ------------------------------------------------------------- %
        
        function tuning_port_ids = get_tuning_load_ids_(obj)
            % method get_tuning_load_ids_() returns the tuning load ids
            
            % get port and load info 
            port = obj.coil.port;
            load = [port(:).tuning_param];
            
            % get matching load mask
            tuning_mask = [load.tuning_flag];
            
            % find feed port ids
            tuning_port_ids = find(tuning_mask);
            
        end
        
        % ------------------------------------------------------------- %

        function [] = tune_coil_(obj, working_freq, optimizer)
            % method performes optimization of tuning loads s.t. nominal
            % range constraints; minimizes the Frobenious norm of the 
            % reduced scattering matrix
            
             % assemble SIE
            obj.assemble_SIE_system_(working_freq, 1);
            
            % compute network parameters at the specified frequency
            obj.compute_network_params_(working_freq, 1);
            
            
            % get tuning and matching port ids
            tuning_port   = obj.get_tuning_load_ids_();
                        
            % extract port and load structures
            loads         = [obj.coil.port.load];
            tuning_params = [loads(tuning_port).tuning_param];
            
            % get constraints on nominal values for tuning loads
            av_optim = [tuning_params.init_val];
            lb_optim = [tuning_params.min_val];
            ub_optim = [tuning_params.max_val];
            
            % optimization flags
            options = optimoptions('fmincon','Algorithm','active-set','MaxIter',1000,'TolFun',1e-16,...
                                   'StepTolerance',1e-16,'MaxFunctionEvaluations',inf,...
                                   'ConstraintTolerance',1e-16,'TypicalX',ones(size(tuning_port))*1e-16);
                               
            % run optimization problem
            [lumped_loads, ~] = fmincon(@(av_optim)obj.tuning_objective_func_(working_freq, av_optim),...
                                   av_optim,[],[],[],[],lb_optim,ub_optim,[],options);
            
            % loads                   
            for i = 1:length(tuning_port)
                
                % update port settings
                obj.coil.port(tuning_port(i)).load.value  = lumped_loads(i);
                obj.coil.port(tuning_port(i)).load.tuning = 0;
                
            end
                               
        end
        
        % ------------------------------------------------------------- %

        function [Smp_norm] = tuning_objective_func_(obj, working_freq, av_optim)
            % method tuning_objective_func_() computes the objective
            % function for current tuning load values
            
            % get number of network ports
            N = size(obj.Sport_,1);
            
            % get tuning port ids 
            tuning_ports   = obj.get_tuning_load_ids_();
            matching_ports = obj.get_matching_port_ids_();
            
            % get local indices of tuning portster
            tuning_indices   = find(ismember(obj.feed_tune_ports_, tuning_ports));
            matching_indices = find(ismember(obj.feed_tune_ports_, matching_ports));
            
            % allocate terminating loads
            Zt = {};
            Z0 = 50;

            % init terminating loads with 50 Ohm impedance
            Zt(1:N) = {Z0};

            
            % update 
            for tuning_id = 1:length(tuning_ports)
                
                tuning_port  = tuning_ports(tuning_id);
                tuning_index = tuning_indices(tuning_id);
                
                % update termin
                Zt{tuning_index} = obj.port_impedance_(tuning_port, working_freq, av_optim(tuning_id));
            end
            
            % reduce dimensions of the multiport network to the number of
            % feed ports 
            Smp = snp2smp(obj.Sport_, Z0, matching_indices, Zt);
            
            % compute the resulting frobenius norm 
            Smp_norm = norm(Smp, 'fro');
            
        end
        
        % ------------------------------------------------------------- %

        function update_solution_(obj)
            
            % set surface currents 
            obj.task_solution_.set_surface_currents(obj.Jcoil);
            
            % set network parameters
            obj.task_solution_.set_Network_param(obj.NP);

        end 
        
        % ------------------------------------------------------------- %


        function run_frequency_sweep_(obj)
                                  
             % get SIE properties/parameters
            [freqs, N_freq] = obj.get_working_frequencies_();
            
            % if not tuned run tuning 
            if ~isempty(obj.get_tuning_load_ids_())
                obj.tune_coil_(freqs);
            end
            
            % reset all parameters
            obj.init_solution_();
                        
            % loop over frequencies
            for i_freq = 1:N_freq
                
                % get current frequency
                freq = freqs(i_freq);
                            
                % assemble SIE
                obj.assemble_SIE_system_(freq);
                
                % find surface currents for current 
                obj.compute_surface_currents_(i_freq);
                
            end

            % update solution object
            obj.update_solution_();
            
        end

    end
    
end