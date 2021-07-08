classdef Einc_Basis < handle
    
    
    properties
        
        basis_file_name
        
        Dcoord
        Xin
        Pin
        U
        
        dom_basis
        index_basis
    end
    
    % --------------------------------------------------------- %
    methods
        
        function obj = Einc_Basis(basis_fname)
            
            % get name of file with basis
            obj.basis_file_name = basis_fname;
            
        end

        function  load_basis(obj)
            
            m = matfile(obj.basis_file_name);
            
            obj.Dcoord      = m.Dcoord;
            obj.Xin         = m.Xin;
            obj.Pin         = m.Pin;
            obj.U           = m.U;
            obj.dom_basis   = m.dom_basis;
            obj.index_basis = m.index_basis;
        end
        
    end
    
end