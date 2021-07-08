function[T_gpu] = to_GPU(T,flag_tucker)

    % 
    [J,K] = size(T);

    switch flag_tucker

        case 0

            T_gpu = gpuArray(T);

        case 1

            T_gpu = struct('G',repmat({[]}, J, K), ...
                          'U1',repmat({[]}, J, K), ...
                          'U2',repmat({[]}, J, K), ...
                          'U3',repmat({[]}, J, K));
            for j=1:J
                for k = 1:K
                    T_gpu(j,k).G = gpuArray(T(j,k).G);
                    T_gpu(j,k).U1 = gpuArray(T(j,k).U1);
                    T_gpu(j,k).U2 = gpuArray(T(j,k).U2);
                    T_gpu(j,k).U3 = gpuArray(T(j,k).U3);
                end
            end

    end

end