%%    MARIE - MAgnetic Resonance Integral Equation suite
%
%     Main MARIEv2.0 script
%     Run to start GUI
%
%% DATA
%       RHBM  structure
%           name - name of realistic human body model
%           r - mapping of the internal edge number to dof number
%           epsilon_r - voxel dielectric
%           sigma_e  - voxel conductivity
%           rho  - voxel density
%           idxS - indexes of non-air material voxels
%       COIL structure
%           name - name of coil model
%           type - kind of coil model, wire of surface
%           index - mapping of the internal edge number to dof number
%           etod - etod numbering of elements and edges
%           node - coordinates of the nodes 
%           edge - numbering of edges
%           elem - 3 indexes of the nodes defining an element
%           Index_elem - mapping of index to elements
%           port - port definition
%       SOL structure
%           to be used
%       INC structure 
%           name - name of external source
%           Einc - Incident electric field         
%           Hinc - Incident magnetic field
%           r = RHBM.r
% _________________________________________________________________________

clc;
clear all;
close all force;
addpath(genpath(pwd));

MARIEv2;