function [Jtot] = Jcb2Jtot_mvp(Jcb,PS)
% function [Jtot] = mvp_Jcb2Jext(Jc, Jb, idxS_ext, P)

%% INPUT
%       Jc       - SCOIL surface currents  
%       Jb       - volumetric body currents [Lx x My x Nz x l x q]
%       idxS_ext - indicies of body voxels in the extended domain
%       P        - Projector from surface currents to volumetric
%                  representation currents withing extended domain
%
%
%% OUTPUT
%       Jtot - the resulting volumetric currents after projection
%  
% -------------------------------------------------------------------------
%
% _________________________________________________________________________

    
%check if Jext is set properly
if(isempty(Jcb))
    error('Jc or Jb are empty! \n');
end
    
% total current
Jtot = PS * Jcb;



