% vie_assembly
% 
%  Assembles the Toeplitz operators used in J-VIE 
%
%  INPUTS:
%    op   -  'N' or 'K' (char scalar)
%    freq -  frequency of operation (double scalar)
%    res  -  resolution (double, 1 for isotropic, 3 for anisotropic)
%    dims -  dimensions of domain (double, vector of len. 3)
%    pwl  -  (Optional) evaluate PWL (boolean scalar, true for PWL, false for PWC)
%
%  NOTE:
%    * The total number of elements must not exceed 2^32-1.
%    * Anisotropic resolutions have not been tested. Use at your own risk.
%
%  EXAMPLES:
%    * N operator, PWL, anisotropic
%        vie_assembly('N',297.2e6,[1e-3,1e-3,2e-3],[64,64,64])
%    * K operator, PWL, isotropic
%        vie_assembly('K',297.2e6,1e-3,[64,64,64])
%    * N operator, PWC, isotropic
%        vie_assembly('N',297.2e6,1e-3,[64,64,64],false)
%
%  AUTHOR:
%    José Serrallés (jecs@mit.edu)
