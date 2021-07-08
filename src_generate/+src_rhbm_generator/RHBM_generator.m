classdef RHBM_generator
    
    % Class for the generation of the body model structure
    % Author: Ilias I. Giannakopoulos, Cambridge, MA, 2019
    
    methods
        
        function obj = RHBM_generator()
            
        end

        function[rhbm_out] = RHBM_4comp_cuboid(obj,res,Len,Cnt,e_r,s_e,dens)
            % Generates a four hmogeneous compartment cuboid phantom
            % Len: lengths of edges
            % Cnt: center of the realistic model
            % er: Relative permittivities in each compartment (back left,
            % front left, back right, front right)
            % se: Conductivities in each compartment (back left,
            % front left, back right, front right)
            % dens: Proton densities in each compartment (back left,
            % front left, back right, front right)
            % OUTPUT: rhbm_out: the real
            % OUTPUT: rhbm_out: the realistic body structure
            rhbm_out = src_rhbm_generator.RHBM_4comp_cuboid_(obj,res,Len,Cnt,e_r,s_e,dens);
        end

        function[rhbm_out] = RHBM_cuboid(obj,res,Len,Cnt,e_r,s_e,dens)
            % Generates a homogeneous cuboid phantom
            % INPUT: res: voxel isotropic resolution
            % Len: lengths of edges
            % Cnt: center of the realistic model
            % er: Relative permittivity
            % se: Conductivity
            % dens: Proton density
            % OUTPUT: rhbm_out: the realistic body structure
            rhbm_out = src_rhbm_generator.RHBM_cuboid_(obj,res,Len,Cnt,e_r,s_e,dens);
        end

        function[rhbm_out] = RHBM_cylinder(obj,res,Rad,Len,Cnt,e_r,s_e,dens)
            % Generates a homogeneous cylindrical phantom
            % INPUT: res: voxel isotropic resolution
            % Rad: radius of circles
            % Len: height
            % Cnt: center of the realistic model
            % er: Relative permittivity
            % se: Conductivity
            % dens: Proton density
            % OUTPUT: rhbm_out: the realistic body structure
            rhbm_out = src_rhbm_generator.RHBM_cylinder_(obj,res,Rad,Len,Cnt,e_r,s_e,dens);
        end

        function[rhbm_out] = RHBM_sphere(obj,res,Rad,Cnt,e_r,s_e,dens)
            % Generates a homogeneous spherical phantom
            % INPUT: res: voxel isotropic resolution
            % Rad: radius of sphere
            % Cnt: center of the realistic model
            % er: Relative permittivity
            % se: Conductivity
            % dens: Proton density
            % OUTPUT: rhbm_out: the realistic body structure
            rhbm_out = src_rhbm_generator.RHBM_sphere_(obj,res,Rad,Cnt,e_r,s_e,dens);
        end

        function[rhbm_out] = RHBM_ellipsoid(obj,res,Rad,Cnt,e_r,s_e,dens)
            % Generates a homogeneous ellipsoid phantom
            % INPUT: res: voxel isotropic resolution
            % Rad: lengths of semi-axis
            % Cnt: center of the realistic model
            % er: Relative permittivity
            % se: Conductivity
            % dens: Proton density
            % OUTPUT: rhbm_out: the realistic body structure
            rhbm_out = src_rhbm_generator.RHBM_ellipsoid_(obj,res,Rad,Cnt,e_r,s_e,dens);
        end

        function[rhbm_out] = RHBM_stadium(obj,res,Rad,side,Len,Cnt,e_r,s_e,dens)
            % Generates a homogeneous stadium phantom
            % INPUT: res: voxel isotropic resolution
            % Rad: radius of semicircles
            % side: length of the edge of the rectangle
            % Len: height
            % Cnt: center of the realistic model
            % er: Relative permittivity
            % se: Conductivity
            % dens: Proton density
            % OUTPUT: rhbm_out: the realistic body structure
            rhbm_out = src_rhbm_generator.RHBM_stadium_(obj,res,Rad,side,Len,Cnt,e_r,s_e,dens);
        end

        function[rhbm_out] = RHBM_RHBM(~,s,slice)
            % Crops a Virtual family model and generates a realistic human body model
            % OUTPUT: rhbm_out: the realistic body structure
            % INPUT: s: string with the .mat file of the model, it should
            % contain three variables, r the position, epislon_r the
            % relative perimittivity tensor and sigma_e the conductivity
            % tensor. The code works for Billie, Duke, and Ella models and it needs
            % to be stored with the following name i.e. billie_2_mm_127.74MHz
            % If the code returns an unexpected cut, i.e. legs instead of
            % the head, shift your model in the relevant direction and
            % retry
            % slice: string that picks the desired part of the realistic human body model
            % possible options: 'Full Body','Head','Head and
            % Shoulders','Left Leg', 'Right Leg'
            rhbm_out = src_rhbm_generator.RHBM_RHBM_(s,slice);
        end

        function[rhbm_out] = RHBM_RHBM_pad(~,s,slice,pos,er,se)
            % Crops a Virtual family model and generates a realistic human body
            % model by attaching a hmogeneous dielectric pad on one side of the
            % head
            % INPUT: s: string with the .mat file of the model, it should
            % contain three variables, r the position, epislon_r the
            % relative perimittivity tensor and sigma_e the conductivity
            % tensor. The code works for Billie, Duke, and Ella models and it needs
            % to be stored with the following name i.e. billie_2_mm_127.74MHz
            % If the code returns an unexpected cut, i.e. legs instead of
            % the head, shift your model in the relevant direction and
            % retry
            % slice: string 'Head with Pad'
            % pos: string that defines the position of the pad relevant to
            % the head, 'Left','Right'
            % er: Relative permittivity of the pad
            % se: Conductivity of the pad
            % OUTPUT: rhbm_out: the realistic body structure
            rhbm_out = src_rhbm_generator.RHBM_RHBM_pad_(s,slice,pos,er,se);
        end

        function[rhbm_out] = RHBM_displacement(~,Cnt,rhbm)
            % Displaces the realistic body model
            % INPUT: Cnt: cartesian axis displacement
            % rhbm: realistic body model
            % OUTPUT: rhbm_out: the realistic body structure
            rhbm_out = src_rhbm_generator.RHBM_displacement_(Cnt,rhbm);
        end

        function[rhbm_out] = RHBM_User_Body_Model(obj)
            % Generates a user defined realistic body model
            % OUTPUT: rhbm_out: the realistic body structure
            rhbm_out = src_rhbm_generator.RHBM_User_Body_Model_(obj);
        end
        
        function[Ae] = RHBM_generatedomain(~,res,x,y,z)
            % Returns the area of a triangle
            % INPUT: r1,r2,r3: the coordinates of each node
            % OUTPUT: Ae: the area of the triangle
            Ae = src_rhbm_generator.generatedomain(res,x,y,z);
        end
        
        
        function [] = generate_test_bodies(obj)
            
        end

    end

end