% compute inductance of the two-port coil
loop.L = 0.17;
loop.W = 0.055;
loop.T = 0.01;

mu_0 = 4 * pi * 1e-7;

scaling = mu_0 / (4 * pi );
scaling = 1e-7;


% define magnetic induction 

I_cd_1 = @(x,y) (loop.L - y) / sqrt((loop.W - x).^2 + (loop.L - y).^2);
I_cd_2 = @(x,y) (loop.L + y) / sqrt((loop.W - x).^2 + (loop.L + y).^2);


I_bc_1 = @(x,y) (loop.W - x) / sqrt((loop.W - x).^2 + (loop.L - y).^2);
I_bc_2 = @(x,y) (loop.W + x) / sqrt((loop.W + x).^2 + (loop.L - y).^2);


B_cd = @(x,y) scaling * (I_cd_1(x,y) + I_cd_2(x,y)) / (loop.W - x);
B_bc = @(x,y) scaling * (I_bc_1(x,y) + I_bc_2(x,y)) / (loop.L - y);


%%



orders = [1:32,64];

L = zeros(length(orders),1);

for order = 1:length(orders)

    [w, x] = gauss_1d(orders(order));
    
    %%
    
    frame_W = loop.W - loop.T / 4;
    frame_L = loop.L - loop.T / 4;
    
    
    w_x = frame_W * w;
    w_y = frame_L * w;
    
    xx = x * frame_W;
    yy = x * frame_L;
    
    %%
    
    
    for i = 1:length(x)
        
        x_i = xx(i);
        
        for j = 1:length(x)
            
            y_j = yy(j);
            
            L(order) = L(order) + 2 * w_x(i) * w_y(j) * (B_cd(x_i, y_j) + B_bc(x_i, y_j));
            
        end
    end
    
end

%% 

close all;


figure(); plot(L);
