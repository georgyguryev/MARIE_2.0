classdef COIL_generator
    
    % Class for the generation of the coil structure
    % Author: Ilias I. Giannakopoulos, Cambridge, MA, 2019
    
    methods
        
        function obj = COIL_generator()
            
        end

        function[coil_out] = COIL_OctaCoil(~,a,b,len,t,meshing,shield)
            % Generates an octagonal elliptical birdcage coil with 8 ports
            % INPUT: a: x semi axis 
            % b: y semi axis
            % len: length of the coil
            % t: width
            % meshing: size of the triangular elements
            % shield: struct flag, radius, length, meshing
            % OUTPUT: coil_out: the coil structure
            coil_out = src_coil_generator.COIL_OctaCoil_(a,b,len,t,meshing,shield);
        end
        
        function[coil_out] = COIL_EllipticalTriangularCoil(~,a,b,len,t,meshing,shield)
            % Generates an octagonal elliptical Triangular coil with 8 ports
            % INPUT: a: x semi axis 
            % b: y semi axis
            % len: length of the coil
            % t: width
            % meshing: size of the triangular elements
            % shield: struct flag, radius, length, meshing
            % OUTPUT: coil_out: the coil structure
            coil_out = src_coil_generator.COIL_EllipticalTriangularCoil_(a,b,len,t,meshing,shield);
        end

        function[coil_out] = COIL_OverlapCoil(~,rad1,rad2,over,len,t,meshing,shield)
            % Generates an overlapping coil with 8 ports
            % INPUT: r1: inner radius 
            % r2: outer radius
            % over: overlapping width
            % len: length of the coil
            % t: width
            % meshing: size of the triangular elements
            % shield: struct flag, radius, length, meshing
            % OUTPUT: coil_out: the coil structure
            coil_out = src_coil_generator.COIL_OverlapCoil_(rad1,rad2,over,len,t,meshing,shield);
        end

        function[coil_out] = COIL_SelfCoil(~,rad,dist,len,t,meshing,shield)
            % Generates a self decoupled coil with 8 ports
            % INPUT: rad: radius 
            % dist: distance between the nearest neighbors
            % len: length of the coil
            % t: width
            % meshing: size of the triangular elements
            % shield: struct flag, radius, length, meshing
            % OUTPUT: coil_out: the coil structure
            coil_out = src_coil_generator.COIL_SelfCoil_(rad,dist,len,t,meshing,shield);
        end
        
        function[coil_out] = COIL_SingleRowCoil(~,rad1,rad2,dist,len,z,t,meshing,shield)
            % Generates an elliptical single row coil with 8 ports
            % INPUT: rad: radius 
            % dist: distance between the nearest neighbors
            % len: length of the coil
            % z: displacement in the z direction
            % t: width
            % meshing: size of the triangular elements
            % shield: struct flag, radius, length, meshing
            % OUTPUT: coil_out: the coil structure
            coil_out = src_coil_generator.COIL_SingleRowCoil_(rad1,rad2,dist,len,z,t,meshing,shield);
        end

        function[coil_out] = COIL_TriCoil(~,a,le,meshing,shield)
            % Generates a Wiggins triangular coil with 8 ports
            % INPUT: a: radius 
            % le: length of the coil
            % meshing: size of the triangular elements
            % shield: struct flag, radius, length, meshing
            % OUTPUT: coil_out: the coil structure
            coil_out = src_coil_generator.COIL_TriCoil_(a,le,meshing,shield);
        end

        function[coil_out] = COIL_WigginsCoil(~,r,w,lex,t,meshing,shield)
            % Generates a Wiggins stadium coil with 8 ports
            % INPUT: r: radius of semicircles 
            % w: edge of the rectangle
            % lex: length of the coil
            % t: width
            % meshing: size of the triangular elements
            % shield: struct flag, radius, length, meshing
            % OUTPUT: coil_out: the coil structure
            coil_out = src_coil_generator.COIL_WigginsCoil_(r,w,lex,t,meshing,shield);
        end

        function[coil_out] = COIL_Displacement(~,Cnt,rot,coil_out)
            % Displaces the coil structure
            % INPUT: coil_out: Initial coil structure
            % Cnt: Cartesian displacement
            % rot: rotation in degrees
            % OUTPUT: coil_out: the coil structure
            coil_out = src_coil_generator.COIL_Displacement_(Cnt,rot,coil_out);
        end

        function generate_test_coils(obj, translation)
            
            if nargin < 1
                translation = 0;
            end
            
            % generate half-wave dipole 
            
            if ~exist('2port_asym_dg.txt')
                
                % generate two-port loop
                src_coil_generator.COIL_TwoPortCoil(translation);
            end
            
        end

    end

end