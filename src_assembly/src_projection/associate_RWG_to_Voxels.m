function [RWG_to_PWX_centers, RWG_to_PWX_cells, RWG_PWX_list_near, RWG_RWG_list_near] = associate_RWG_to_Voxels(SCOIL, RHBM, PROJ)
% function looks for the closest voxel expansion cell 
% for a given RWG basis

% get frame width (boundary thickness) of near zone
frame_width = PROJ.frame_width;

% 
N_exp_1d     = PROJ.Nexp_1D;
N_near_1d    = 3 * N_exp_1d + 2 * frame_width;
Nexp         = N_exp_1d.^3;
N_near       = N_near_1d.^3;

[L,M,N, ~] = size(RHBM.r);

%% collect data for mapping 

% 
Nsie = max(SCOIL.index);

% init output data
RWG_to_PWX_centers = zeros(Nsie, 3);
RWG_to_PWX_cells   = zeros(Nsie, Nexp);
RWG_PWX_list_near  = zeros(Nsie, N_near);
RWG_RWG_list_near  = [];
RWG_RWG_near = cell(Nsie,1);

% extract coordinates of (terminal + inner) RWG centers
[RWG_cntr] = get_RWG_centers(SCOIL);

%% find nearest center of expansion voxel (for given edge)
for i_sie = 1:Nsie
 
    % get current RWG center
    r_sie = RWG_cntr(i_sie,:).';
    
    % associate current RWG with the nearest voxel 
    [idx_x_ctr, idx_y_ctr, idx_z_ctr]  = find_nearest_voxel(r_sie, RHBM);

    % assign an expansion cell for given RWG function
    RWG_to_PWX_cells(i_sie,:)   = assign_expansion_cell(idx_x_ctr,...
                                                    idx_y_ctr, idx_z_ctr,...
                                                    Nexp, L, M);
    % form a simple near interaction list                                            
    RWG_PWX_list_near(i_sie,:)  = form_expansion_near_list(idx_x_ctr,...
                                                    idx_y_ctr, idx_z_ctr,...
                                                    Nexp, frame_width,...
                                                    L, M, N, RHBM.r);
               

    % store central expansion voxels
    RWG_to_PWX_centers(i_sie,:) = [idx_x_ctr, idx_y_ctr, idx_z_ctr];
    
    %% refactor form_RWG_near_list method  for pure pFFT (without SIE)
    % Get a list of RWG basises that are near to current 
%     Near_RWG = form_RWG_near_list(SCOIL, RHBM, Nexp, i_sie,...
%                                   idx_x_ctr, idx_y_ctr, idx_z_ctr);
                              
%     RWG_RWG_near{i_sie} = form_RWG_near_list(SCOIL, RHBM, Nexp, i_sie,...
%                                   idx_x_ctr, idx_y_ctr, idx_z_ctr);
                              
    % update pair of near interactions
%     RWG_RWG_list_near = [RWG_RWG_list_near; Near_RWG];

end

% % extract unique pairs of near interactions 
% RWG_RWG_list_near = unique(sort(RWG_RWG_list_near,2), 'rows');

% 
% % store resulting list in sparse format
% RWG_RWG_list_near = sparse(RWG_RWG_list_near);







