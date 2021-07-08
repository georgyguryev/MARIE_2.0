function [unique_rwg_near] = unique_near_interactions_RWG(rwg_nears)
% function [unique_rwg_near] = unique_near_interactions_RWG(rwg_near)

N_near = size(rwg_nears,1);

% init 
unique_rwg_near = rwg_nears(1,:);

% loop over all nonunique pairs
for i = 2:N_near
    
    flag = 1;
    pair = rwg_nears(i,:);
    
    % loop over unique pairs
    for j = 1:size(unique_rwg_near,1)
        if (isequal(pair, unique_rwg_near(j,:)) || ... 
           isequal(fliplr(pair),unique_rwg_near(j,:)))
            flag = 0;
            break;
        end;
    end;
    
    % add new unique pair to the list
    if 1 == flag
        unique_rwg_near = [unique_rwg_near; pair];
    end;
end;

keyboard;