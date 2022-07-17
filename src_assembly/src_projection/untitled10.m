clc;
clear all;

N_max = 500;

C = zeros(2*N_max,5);

t_start_pf = tic;

parfor i = 1:N_max
    
    for j = 1:N_max
        C(i,:) = C(i,:)+ ones(1,5) * ( exp(i)^2 * 3 - sin(i));
    end;
end;

t_end_pf = toc(t_start_pf);

t_start_f = tic;
for i = 1:N_max
    for j = 1:N_max
        C(i,:) = C(i,:) +  ones(1,5) * (exp(i)^2 * 3 - sin(i));
    end;
end;

t_end_f = toc(t_start_f);


disp(t_end_pf);
disp(t_end_f);
