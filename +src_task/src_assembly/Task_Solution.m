classdef Task_Solution < handle
    methods
        % define interface methods 
        % --------------------------------------------------------------- %
        
        function [Jb] = get_plrz_currents_vec(obj)
        % get vector of polarisation currents (output from iter. solver)
            Jb = obj.Jb_plrz_vec_; 
        end
        % --------------------------------------------------------------- %
        
        function [Jb_tensor] = get_plrz_currents_tensor(obj)
        % get tensor of polarisation currents (output from iter. solver)
            Jb_tensor = obj.Jb_plrz_tensor_;
        end
        
        % --------------------------------------------------------------- %
        
        function [Jc] = get_surface_currents(obj)
        % get solution for the coil surface currents (vec. or mat.) 
            Jc = obj.Jc_surface_;
        end
        
        % --------------------------------------------------------------- %
        
        function [E_tot] = get_E_total(obj)
        % get solution for total electric fields 
            E_tot = obj.E_tot_;
        end
        
        % --------------------------------------------------------------- %
        
        function [H_tot] = get_H_total(obj)
        % get solution for total magnetic fields 
            H_tot = obj.H_tot_;
        end
        
        % --------------------------------------------------------------- %
        
        function [NP] = get_Network_param(obj)
            NP = obj.Network_param_;
        end
        % --------------------------------------------------------------- %
        
        function [B1_p] = get_B1_p(obj)
            B1_p = obj.B1_plus_;
        end
        
        % --------------------------------------------------------------- %
        
        function [B1_m] = get_B1_m(obj)
            B1_m = obj.B1_minus_;
        end
        
        % --------------------------------------------------------------- %
        
        function [P_abs] = get_P_abs(obj)
            P_abs = obj.P_absorbed_;
        end
        
        % --------------------------------------------------------------- %
        
        function set_plrz_currents_vec(obj, Jb)
            obj.Jb_plrz_vec_ = Jb;
        end
        
        % --------------------------------------------------------------- %
        
        function set_plrz_currents_tensor(obj, Jb)
            obj.Jb_plrz_tensor_ = Jb;
        end
        
        % --------------------------------------------------------------- %

        function set_surface_currents(obj, Jc)
            obj.Jc_surface_ = Jc;
        end
        
        % --------------------------------------------------------------- %

        function set_relative_residual(obj, rel_res)
            obj.rel_res_ = rel_res;
        end
        
        % --------------------------------------------------------------- %
        
        function set_residual_vector(obj, res_vec)
            obj.res_vec_ = res_vec;
        end
        
        % --------------------------------------------------------------- %
        
        function set_Network_param(obj, NP)
            obj.Network_param_ = NP;
        end
        
        % --------------------------------------------------------------- %
        
        function set_E_total(obj, E_tot)
            obj.E_tot_ = E_tot;
        end
        
        % --------------------------------------------------------------- %
        
        function set_H_total(obj, H_tot)
            obj.H_tot_ = H_tot;
        end
        
        % --------------------------------------------------------------- %
        
        function set_lapse_time(obj, lapse_time)
            obj.lapse_time_ = lapse_time; 
        end
        
        % --------------------------------------------------------------- %
        
        function set_B1_plus(obj, B1_plus)
            obj.B1_plus_ = B1_plus;
        end
        
        % --------------------------------------------------------------- %
        
        function set_B1_minus(obj, B1_minus)
            obj.B1_minus_ = B1_minus;
        end
        
        % --------------------------------------------------------------- %
        
        function set_absorbed_power(obj, P_absorbed)
            obj.P_absorbed_ = P_absorbed;
        end
        % --------------------------------------------------------------- %
        
        function set_SAR(obj, SAR)
            obj.SAR_ = SAR;
        end
        % --------------------------------------------------------------- %
        
        function set_epsilon_r(obj, epsilon_r)
            obj.epsilon_r_ = epsilon_r;
        end
        % --------------------------------------------------------------- %
        
        function set_sigma_e(obj, sigma_e)
            obj.sigma_e_ = sigma_e;
        end
        
        % --------------------------------------------------------------- %
        
        function set_freqs(obj, freqs)
            obj.freqs_ = freqs;
        end
        
        % --------------------------------------------------------------- %
        
        function init_SIE_solution(obj, dims)
            obj.Network_param_ = Network_Parameters(dims);
        end
        
        % --------------------------------------------------------------- %
       
        function init_VSIE_solution(obj,dims)
            
            obj.lapse_time_ = zeros(dims.N_feeds, dims.N_freqs);
            obj.res_vec_    = cell(dims.N_feeds, dims.N_freqs);
            obj.rel_res_    = zeros(dims.N_feeds, dims.N_freqs);
            obj.E_tot_      = zeros([dims.vie,dims.N_feeds,dims.N_freqs]);
            obj.H_tot_      = zeros([dims.vie,dims.N_feeds,dims.N_freqs]);
            
        end
        
        % --------------------------------------------------------------- %
        
        function ave_time = average_runtime(obj)
            
            ave_time = round(mean(obj.lapse_time_), 3);
        end
        
        % --------------------------------------------------------------- %
        
        function ave_iter = average_iteration_count(obj)
            
            N_ports = size(obj.res_vec_,1);
            
            iter_count = 0;
            for i = 1:N_ports
                iter_count = iter_count + size(obj.res_vec_{i},1);
            end
            
            ave_iter = round(iter_count / N_ports);
        end
        
    end
    
    
        
    % list general solution types
    properties %(Access = private)
        
        lapse_time_
        
        res_vec_
        rel_res_
        
      % electric field measurements
        E_tot_
        P_absorbed_
        SAR_
        
        % magnetic field measurements
        H_tot_
        B1_plus_
        B1_minus_
        
        % network parameters
        Network_param_
        
        % coil surface currents
        Jc_surface_

        % coil polarization currents
        Jb_plrz_vec_
        Jb_plrz_tensor_
        
        % properties
        epsilon_r_
        sigma_e_
        
        % frequencies
        freqs_
    end
   
end