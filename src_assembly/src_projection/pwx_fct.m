function [weight] = pwx_fct(z1, z2, z3, b_basis)
% function [weight] = pwx_fct(z1, z2, z3, b_basis)
% function returns coefficients for PWC/PWL basis integration

switch b_basis
    case 1
        weight = 1;
    case 2
        weight = z1 / 2;
    case 3
        weight = z2 / 2;
    case 4
        weight = z3 / 2;
end