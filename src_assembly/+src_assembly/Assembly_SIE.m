classdef Assembly_SIE < src_assembly.Assembly_Base    
% Class for coil-specific problem assembly
% Author: Georgy Guryev, Cambridge, MA, 2019
    
    properties
        
         
        Z_L_hat
        Z_ref
        Z_sh
        Z_ser
        Z_mt
        rhs_cp
        w_src_I
        w_coil_V
        Zw_coil
        Zw_src
        Yw_coil
    end
    
% ====================================================================== %

    methods 
        function [obj] = Assembly_SIE(task_settings, coil, dims)
            
            % call superclass constructor
            obj = obj@src_assembly.Assembly_Base(task_settings, dims);

            % store coil in properties
            obj.coil= coil;

            % get prev coil
            prev_coil = obj.setgetPrevCoil();
            
            % if the prev coil does not exist or is 
            if isempty(prev_coil) 
                obj.setgetPrevCoil(coil);
            end
        end
        
        % --------------------------------------------------------------- %
        
        function assemble_system(obj, freq)
            
            % assemble edge2edge SIE matrix 
            obj.assemble_SIE_system_(freq)
            
            % take into account external circuitry
            obj.add_matching_impedance_(freq);
            
            % form rhs taking into account ext. circuitry
            obj.form_rhs_sie_();
            
            % update impedance matrix and inverse impedance matrix;
            obj.Zcoil_     = obj.Zcoil_ + obj.Z_mt;
            obj.Zcoil_inv_ = eye(size(obj.Zcoil_)) / obj.Zcoil_;
            
        end

    end
% ====================================================================== %
    
    methods (Sealed = true, Access = protected)
           
        % function assembles Z_coil_ and F_coil_
        function  assemble_SIE_system_(obj, freq)
            
            % if SIE runs at single frequency
            if nargin < 2
                freq = 3e8;
            end
            
            % get previous settings
            prevFreq  = obj.setgetPrevFreq_SIE();
            prevCoil  = obj.setgetPrevCoil();
            [prevZcoil,~,~] = obj.setgetPrevZcoil();
            
            % if the coil and frequency are the same skip assembly
            if ~isequal(prevCoil.elem, obj.coil.elem) || any(prevFreq ~= freq) || isempty(prevZcoil)
                
 
                % call some assebly method
                [obj.Zcoil_, obj.Fcoil_, obj.feed_tune_ports_] = Assembly_SIE_par(obj.coil, freq);
                
                % update persistent variables
                obj.setgetPrevFreq_SIE(freq);
                obj.setgetPrevCoil(obj.coil);
                obj.setgetPrevZcoil(obj.Zcoil_, obj.Fcoil_, obj.feed_tune_ports_);
            else
                [obj.Zcoil_, obj.Fcoil_, obj.feed_tune_ports_] = obj.setgetPrevZcoil();
            end

        end
        
        % ------------------------------------------------------------- %
        
        function [] = add_matching_impedance_(obj, freq)
            % function add_matching_impedance_() constructs a sparse matrix
            % of impedance contributions due to matching loads
            
            % get ports and find matching
            port_ids = obj.get_matching_port_ids_();
            
            % allocate memory for matching contribution
            obj.Z_L_hat   = zeros(obj.dims.N_feeds);
            obj.Z_ref     = zeros(obj.dims.N_feeds);
            obj.Z_sh      = zeros(obj.dims.N_feeds);
            obj.rhs_cp    = zeros(obj.dims.N_feeds);            
            obj.w_src_I   = zeros(obj.dims.N_feeds);
            obj.w_coil_V  = zeros(obj.dims.N_feeds);
            
            % loop over matching ports
            for pnum = 1:obj.dims.N_feeds
                
                % get current load port
                port_num = port_ids(pnum);
                
                [Z_L_hat_i, rhs_c_i, w_src_i, w_coil_i, Zw_coil_i, Zw_src_i, Yw_coil_i] = obj.get_load_values_(port_num, freq);
                
                obj.Z_L_hat(pnum, pnum)  = Z_L_hat_i;
                obj.rhs_cp(pnum, pnum)   = rhs_c_i;
                obj.w_src_I(pnum, pnum)  = w_src_i;
                obj.w_coil_V(pnum, pnum) = w_coil_i;
                obj.Zw_coil(pnum, pnum)  = Zw_coil_i;
                obj.Zw_src(pnum, pnum)   = Zw_src_i;
                obj.Yw_coil(pnum, pnum)  = Yw_coil_i;
            end
            
            obj.Z_mt  = -obj.Fcoil_* obj.Z_L_hat * obj.Fcoil_.';
            
        end
        
        % ------------------------------------------------------------- %
        function form_rhs_sie_(obj)
            
            % define rhs
            obj.rhs_c = obj.Fcoil_ * obj.rhs_cp;
            
        end
        
        % ------------------------------------------------------------- %
        
        function  [Z_L_hat, rhs_c, w_src_I,w_coil_V, Zw_coil, Zw_src, Yw_coil] = get_load_values_(obj, port_id, freq)
            % method port_impedance_(id) returns the nominal value of the
            % load at id-th port
            
            % get load at the id-th port
            port = obj.coil.port(port_id);
            load_type = port.load_type;
                       
            switch load_type
                
                case 'shunt'
                    
                    Z_ref = port.ref_impedance;
                    
                    Z_sh = obj.LE_impedance_(port.matching_param.LE_1,...
                        port.matching_param.val_1, freq);
                    
                    Z_ser = 0;
                    
                    Z_ref_hat = Z_ref;
                    
                    Z_L = Z_ref * Z_sh / (Z_sh + Z_ref);
                    
                    Z_L_hat = Z_L; 
                    w_src_I = 1;
                    w_coil_V = 1;
                    Zw_coil = 0;
                    Zw_src  = 0;
                    Yw_coil = 1 / Z_sh;
                    rhs_c = Z_L / Z_ref_hat;

                    
                case 'shunt-series'
                    
                    Z_ref = port.ref_impedance;
                    
                    Z_sh = obj.LE_impedance_(port.matching_param.LE_1,...
                        port.matching_param.val_1, freq);
                    
                    Z_ser = obj.LE_impedance_(port.matching_param.LE_2,...
                        port.matching_param.val_2, freq);
                    
                    Z_ref_hat = Z_ref;
                    
                    Z_L = Z_ref * Z_sh / (Z_sh + Z_ref);
                    
                    Z_L_hat = Z_L + Z_ser;
                    w_src_I = (Z_sh + Z_ser) / Z_sh;
                    w_coil_V = 1;
                    Zw_coil = Z_ser;
                    Zw_src  = 0;
                    Yw_coil = 1 / Z_sh;
                    rhs_c = Z_L / Z_ref_hat;

                    
                case 'series-shunt'
                    
                    Z_ref = port.ref_impedance;
                    
                    Z_sh = obj.LE_impedance_(port.matching_param.LE_2,...
                        port.matching_param.val_2, freq);
                    
                    Z_ser = obj.LE_impedance_(port.matching_param.LE_1,...
                        port.matching_param.val_1, freq);
                    
                    Z_ref_hat = Z_ref + Z_ser;
                    
                    
                    Z_L = Z_ref_hat * Z_sh / (Z_sh + Z_ref_hat);
                    
                    Z_L_hat = Z_L;
                    w_src_I = 1;
                    w_coil_V = 1;
                    
                    Zw_coil = 0;
                    Zw_src  = Z_ser;
                    Yw_coil = 1 / Z_sh; 
                    rhs_c = Z_L / Z_ref_hat;

                    
                case 'shunt-series-shunt'
                    Z_ref = port.ref_impedance;
                    
                    Z_sh = obj.LE_impedance_(port.matching_param.LE_1,...
                        port.matching_param.val_1, freq);
                    
                    Z_ser = obj.LE_impedance_(port.matching_param.LE_2,...
                        port.matching_param.val_2, freq);
                    
                    Z_par = obj.LE_impedance_(port.matching_param.LE_3,...
                        port.matching_param.val_3, freq);
                    
                    Z_tilde = Z_ref + Z_ser + Z_ser * Z_ref / Z_sh;
                    
                    w = 1 / (1 + Z_tilde / Z_par + Z_ref / Z_sh);
                    
                    Z_L_hat = Z_tilde * w;
                    Zw_coil = Z_ser;
                    Zw_src  = 0; 
                    Yw_coil = (Z_sh + Z_par + Z_ser) / (Z_par * Z_sh);
                    
                    w_coil_V = 1 + Z_ser / Z_par;
                    w_src_I = (Z_sh + Z_ser) / Z_sh;
                    rhs_c = w;
                    
                    
                otherwise
                    error('Unknown load type!');
            end
            

        end
        
        % ------------------------------------------------------------- %
        
        function Z_LE = LE_impedance_(obj, LE, val, freq)
            
            % pick load type
            switch LE
                case 'resistor'
                    Z_LE = val;
                case 'inductor'
                    Z_LE = 2 * pi * freq * 1i * val;
                case 'capacitor'
                    Z_LE = 1 / (2 * pi * freq * 1i * val);
                otherwise
                    error('Unknown load type!');
            end
        end
        
        % ------------------------------------------------------------- %
        
        function port_ids = get_matching_port_ids_(obj)
            % function get_matching_port_ids_() returns indices of matching
            % ports
            
            % get port and load info 
            port = obj.coil.port;
            load = strcmp({port(:).type}, 'feed');
            
            % find feed port ids
            port_ids = find(load);
            
        end
 
    end
    
% ====================================================================== %
    
end