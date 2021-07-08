classdef Visualizer < handle
    
    % Class for the Visualization of the electromagnetic field related
    % measurements and the geometry of the problem
    % Author: Ilias I. Giannakopoulos, Cambridge, MA, 2019
    
    
    properties 
        
        % variable keeps track of the dimensions 
        dims
        
    end
    
    methods
        
        function obj = Visualizer()
            
        end

        function visualize_normals(~,coil)
            % Visualizes the normals of the triangular mesh of the coil structure
            % INPUT: COIL: a struct that contains the information for the
            % coil structure
            src_visualizer.visualize_normals_(coil);
        end

        function visualize_coil(~,coil)
            % Visualizes the coil structure and the ports
            % INPUT: COIL: a struct that contains the information for the
            % coil structure
            src_visualizer.visualize_coil_(coil);
        end
        
        function visualize_body(~,xd,yd,zd,id)
            % Visualizes the body model structure
            % INPUT: RHBM: a struct that contains the information for the
            % body model structure
            src_visualizer.visualize_body_(xd,yd,zd,id);
        end
        
        function visualize_coil_and_body(~,xd,yd,zd,id,coil)
            % Visualizes the body and coil models structure
            % INPUT: RHBM: a struct that contains the information for the
            % body model structure
            % COIL: a struct that contains the information for the
            % coil structure
            src_visualizer.visualize_coil_and_body_(xd,yd,zd,id,coil)
        end
        
        function visualize_properties(~,EP,sac,cut)
            % Visualizes a slice of one electrical property of the body model
            % INPUT: EP the electrical property (relative permittivity or
            % conductivity)
            % sac: String: Sagittal, Coronal or Axial
            % cut: Cropping point
            src_visualizer.visualize_properties_(EP,sac,cut)
        end
        
        function h = visualize_fields(~,E,sac,cut,port_cut,freq_cut)
            % Visualizes the slice of an electromagnetic field measurement magnitude
            % INPUT: E: Electromagnetic field measurement
            % sac: String: Sagittal, Coronal or Axial
            % cut: Cropping point
            % port_cut: Port to visualize for
            % freq_cut: Cropping frequency point
            h = src_visualizer.visualize_fields_(E,sac,cut,port_cut,freq_cut);
        end
        
        function visualize_Network_param_matrix(~,E,freq_num, title)
            % Visualizes the matrix of the S,Z or Y paramateres for a frequency point
            % INPUT: E: The S,Z or Y parameter matrices for all frequency points
            % freq_cut: Cropping frequency point
            src_visualizer.visualize_Network_param_matrix(E,freq_num, title)
        end
        
        function show_frequency_sweep(~,E,port_1,port_2,freqs)
            % Visualizes one S,Z or Y paramater S(i,j) for all the frequency points
            % INPUT: E: The S,Z or Y parameter matrices for all frequency points
            % i,j: Port i, port j
            % freqs: The vector of frequencies
            src_visualizer.visualize_frequency_sweep(E,port_1,port_2,freqs)
        end
        
        function visualize_ConvergenceRate(~,E,port_cut,freq_cut)
            % Visualizes the convergence rate of the iterative solver
            % INPUT: E: Tensor of convergence rates for all ports and frequencies
            % port_cut: Port to visualize for
            % freq_cut: Cropping frequency point
            src_visualizer.visualize_ConvergenceRate_(E,port_cut,freq_cut)
        end
        
        function visualize_InterpolatePlot(~,E,sac,cut,port_cut,freq_cut,stepi)
            % Performs an interpolation plot for a slice of an electromagnetic field measurement magnitude
            % INPUT: E: Electromagnetic field measurement
            % sac: String: Sagittal, Coronal or Axial
            % cut: Cropping point
            % port_cut: Port to visualize for
            % freq_cut: Cropping frequency point
            % stepi^3 the number of voxels that each initial voxel is interpolated to
            src_visualizer.visualize_interpolate_plot_(E,sac,cut,port_cut,freq_cut,stepi)
        end
        
        function h = visualize_currents(obj,coil,E,port_num,freq)
            % Visualizes the currents from one port on the coil structure
            % INPUT: COIL: a struct that contains the information for the
            % coil structure
            % INPUT: E: The surface currents on the coil for each port and
            % frequency points
            % port_cut: Port to visualize for
            % freq_cut: Cropping frequency point
            h = src_visualizer.visualize_currents_(obj,coil,E,port_num,freq);
        end
        
        function[Ae] = triangle_area(~,r1,r2,r3)
            % Returns the area of a triangle
            % INPUT: r1,r2,r3: the coordinates of each node
            % OUTPUT: Ae: the area of the triangle
            Ae = src_visualizer.triangle_area_(r1,r2,r3);
        end

    end

end