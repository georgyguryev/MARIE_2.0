function [Jcb] = Jtot2Jcb_mvp(Jext,PS)
% function [Jc,Jb] = mvp_Jtot2Jcb(Jext,P,S)
% function [Jtot] = mvp_Jcb2Jext(Jc, Jb, idxS_ext, P)

%% INPUT
%
%       Jtot - the resulting volumetric currents after projection
%       P        - Projector from surface currents to volumetric
%                  representation currents withing extended domain
%
%
%% OUTPUT
%       
%       Jc       - SCOIL surface currents  
%       Jb       - volumetric body currents [Lx x My x Nz x l x q]
% -------------------------------------------------------------------------


%check if Jext is set properly
if(isempty(Jext))
    error('Jext is empty! \n');
end

% % find coil and body currents
% Jc = P.' * Jext;
% Jb = S.' * Jext;

Jcb = PS.' * Jext;

