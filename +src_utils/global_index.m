function [glob_idx] = global_index(idx_x,idx_y,idx_z,L,M)
% function maps 3 indeces to a single global index using Matlab ordering
% i.e. column major ordering

% map 3 indices to a single one
glob_idx = idx_x + (idx_y-1) * L + (idx_z-1) * L * M;