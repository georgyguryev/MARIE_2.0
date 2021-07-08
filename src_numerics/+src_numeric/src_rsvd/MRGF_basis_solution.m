function [M] = MRGF_basis_solution(M,RHBM,freq,tol,M_fname,GPU_flag)
%
% _________________________________________________________________________
% _________________________________________________________________________
%
%   Function that solve the basis of incident fields in M for a given RHBM
%
%   INPUT:  M         basis with the vectors with incident fields
%           RHBM      realistic human body model struct
%           freq      frequency
%           tol       relative tolerance for the basis truncation
%           M_fname   name for the file of temporary storage of matrix M
%           GPU_flag  1 to use the GPU
%
%   OUTPUT: M         basis vectors with the current solutions
%
% _________________________________________________________________________
%
%
% -------------------------------------------------------------------------
%
%%   This function is part of MARIE
%   MARIE - Magnetic Resonance Integral Equation suite
%           Jorge Fernandez Villena   -- jvillena@mit.edu
%           Athanasios G. Polimeridis -- thanos_p@mit.edu
%           Copyright © 2016
%           RLE Computational Prototyping Group, MIT
% 
%           This software is free and open source
%           Distributed under the GNU-GPLv3 terms
%           For details see MARIE_license.txt
%
% _________________________________________________________________________


% Initialize variables.
if(nargin < 3 )
   fprintf(1, '\n ERROR: not enough arguments\n');
   return
end
if(nargin < 4 || isempty(tol))
   tol = 1e-4;
end
if(nargin < 5 )
   temp_name = [ ];
end
if(nargin < 6 || isempty(GPU_flag))
   GPU_flag = 1;
end

% -------------------------------------------------------------------------
% Initialization of variables
% -------------------------------------------------------------------------

fid = 1;

% -------------------------------------------------------------------------
% Pre-proccess the RHBM and prepare for the VIE solve
% -------------------------------------------------------------------------

tic_i = tic;

[fVIEsolver,Gram,tau3,~,~] = fVIEsolver_Assembly(RHBM,freq,GPU_flag);

tVIE = toc(tic_i);


% -------------------------------------------------------------------------
% Solve the VIE system for each vector in the basis
% -------------------------------------------------------------------------


tic_i = tic;

% allocate space
Ncol = size(M,2);

fprintf(fid, '\n  Solving the system for %d inputs \n', Ncol);
t1 = tic;

for ii=1:Ncol

    % get the incident field excitation and transform into currents
    Jexc = tau3.*(Gram*M(:,ii));
    
    % solve to obtain the total currents
    Jsol = fVIEsolver(Jexc, tol, []);
    
    M(:,ii) = Jsol; % store the solution currents
    
    if ((ii/500) == floor(ii/500))
        fprintf(fid, '\n %d Solves. Elapsed time %g', ii, toc(t1));
        t1 = tic;
    end
    
end

tSolve = toc(tic_i);

clear Jsol; clear Jexc;

fprintf(fid, '\n\n');

fprintf(fid, '\n\n Equivalent current generated for %d vectors. Elapsed time %g \n', Ncol, tSolve);

if (~isempty(M_fname))
    
    % save matrix M to temporary file M_fname
    save(M_fname, 'M', '-v7.3');
    
end

