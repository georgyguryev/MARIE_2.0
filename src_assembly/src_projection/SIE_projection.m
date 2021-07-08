function [P] = SIE_projection(SCOIL, RHBM, f, PWX_type)

% get size of the original VIE problem
[L,M,N,~] = size(RHBM.r);

%% get body dimensions
xb_min = min(RHBM.r(:,1,1,:)); xb_max = max(RHBM.r(:,1,1,:));
yb_min = min(RHBM.r(1,:,1,2)); yb_max = min(RHBM.r(1,:,1,2));
zb_min = min(RHBM.r(1,1,:,3)); zb_max = min(RHBM.r(1,1,:,3));


%% get dimensions of the SCOIL

xc_min = min(SCOIL.node(1,:));  xc_max = max(SCOIL.node(1,:));
yc_min = min(SCOIL.node(2,:));  yc_max = max(SCOIL.node(2,:));
zc_min = min(SCOIL.node(3,:));  zc_max = max(SCOIL.node(3,:));



if ()