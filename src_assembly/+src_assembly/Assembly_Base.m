classdef (Abstract) Assembly_Base < handle
% define base class for system assembly

% Abstract class for EM problem assembly
% Author: Georgy Guryev, Cambridge, MA, 2019

% specify general properties for all assembly subclasses
properties 
    
    % store coil and body models
    coil
    
    %  
    task_settings_
    
    % edge2edge impedance matrix
    Zcoil_ 
    
    Zcoil_inv_
    
    % port2edge projection matrix (i.e. V_t = F * V_p  and I_p = F.' * I_t)
    Fcoil_
    
    % ids of feed and tune ports
    feed_tune_ports_
    
    % assemble K, N, G operators
    operators
    
    rhs
    rhs_Hinc
    
    rhs_c
    rhs_b
    
    solver
    
    % coupling matrix for (explicit coupling)
    Z_bc
    
    % class instance keeps track of all sizes
    dims
        
end

% ====================================================================== %

methods 
    
    % define constructor
    function obj = Assembly_Base(task_settings, dims)
        
        % copy task settings
        obj.task_settings_ = task_settings;
        
        % instantiate operator class EM_operators
        obj.operators = src_operator.EM_operators(dims);
        
        % get a handle on dimentions
        obj.dims = dims;               
    end
    
end

% ====================================================================== %


methods(Access = protected)
    
    % function assembles Z_coil_ and F_coil_
    assemble_SIE_system_(obj)
    
    % assemble VIE/VSIE operators
    assemble_operators_(obj)
    
    % assemble mvp's
    assemble_core_mvps_(obj)
    
    % assemble preconditioners for given problem
    assemble_preconditioners_(obj)
end

% ====================================================================== %

% define abstract method 
methods (Abstract)
       
    % interface method
    assemble_system(obj)
       
end

% ====================================================================== %

methods(Static)
    
    function prevCoil = setgetPrevCoil(coil)
        
        % define persistent variable
        persistent PrevCoil;
        
        % update PrevCoil 
        if nargin
            PrevCoil = coil;
        end
        
        % return previous coil
        prevCoil = PrevCoil;
        
    end
    
    % ------------------------------------------------------------- %

    function prevFreq = setgetPrevFreq_SIE(freq)
        
        persistent PrevFreq_SIE;
        
        if nargin
            PrevFreq_SIE = freq;
        end
       
        prevFreq = PrevFreq_SIE;
    end
    % ------------------------------------------------------------- %
    
    function prevFreq = setgetPrevFreq_VIE(freq)
        
        persistent PrevFreq_VIE;
        
        if nargin
            PrevFreq_VIE = freq;
        end
       
        prevFreq = PrevFreq_VIE;
    end
    
    % ------------------------------------------------------------- %
    
    function [prevZcoil, prevFcoil, prevPorts] = setgetPrevZcoil(Zcoil, Fcoil, ports)
        
        persistent PrevZcoil PrevFcoil PrevPorts;
        
        if nargin
            PrevZcoil = Zcoil;
            PrevFcoil = Fcoil;
            PrevPorts = ports;
        end
       
        prevZcoil     = PrevZcoil;
        prevFcoil     = PrevFcoil;
        prevPorts     = PrevPorts;
    end
    
    % ------------------------------------------------------------- %
    
    function [prevZcoil_inv] = setgetPrevZcoil_inv(Zcoil_inv)
        
        persistent PrevZcoil_inv;
        
        if nargin
            PrevZcoil_inv = Zcoil_inv;
        end
       
        prevZcoil_inv     = PrevZcoil_inv;
    end    
    % ------------------------------------------------------------- %
    
    function prev_N_vie = setget_Prev_N_vie(N_vie)
        
        persistent Prev_N_vie;
        
        if nargin
            Prev_N_vie = N_vie;
        end
        
        prev_N_vie = Prev_N_vie;
    end
    
    % ------------------------------------------------------------- %
    
    function prev_N_ext = setget_Prev_N_ext(N_ext)
        
        persistent Prev_N_ext;
        
        if nargin
            Prev_N_ext = N_ext;
        end
        
        prev_N_ext = Prev_N_ext;
    end
    
    % ------------------------------------------------------------- %

    function prev_dom_vie = setget_PrevDomain_vie(dom_vie)
        
        persistent PrevDomain_VIE;
        
        if nargin
            PrevDomain_VIE = dom_vie;
        end
        
        prev_dom_vie = PrevDomain_VIE;
    end
    
    % ------------------------------------------------------------- %

    function prev_dom_ext = setget_PrevDomain_ext(dom_ext)
        
        persistent PrevDomain_EXT;
        
        if nargin
            PrevDomain_EXT = dom_ext;
        end
        
        prev_dom_ext = PrevDomain_EXT;
    end
    
    % ------------------------------------------------------------- %
    
    function prev_rhs = setget_Prev_rhs(rhs)
        
        persistent Prev_rhs;
        
        if nargin
            Prev_rhs = rhs;
        end
        
        prev_rhs = Prev_rhs;
    end
    
    % ------------------------------------------------------------- %
    function [out_basis] = setget_Prev_basis(in_basis)
        
        persistent Basis;
        
        if nargin
            Basis = in_basis;
        end
        
        out_basis = Basis;
    end

    
    
end
    
end