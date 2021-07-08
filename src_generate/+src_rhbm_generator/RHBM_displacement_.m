function [rhbm] = RHBM_displacement_(Cnt,rhbm)

    % Displaces the realistic body model
    % Author: Ilias I. Giannakopoulos, Cambridge, MA, 2019
    
    rhbm.r(:,:,:,1) = rhbm.r(:,:,:,1)+Cnt(1);
    rhbm.r(:,:,:,2) = rhbm.r(:,:,:,2)+Cnt(2);
    rhbm.r(:,:,:,3) = rhbm.r(:,:,:,3)+Cnt(3);
    
end