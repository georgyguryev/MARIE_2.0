function output = apply_preconditioner(Prec, vec, i_perm)
% function allows to apply two types of preconditioners:
%
% 1) if Prec is a matrix, the preconditioner is treated as not inverted and
% therefore is applied using backslash: (Prec \ residual) or (Prec
% \residual(i_perm))
% 2) if Prec is a vector, it is treated as the inverted preconditioner and
% is applied to the residual vector through ".*": (Prec .* residual) 

[m,n] = size(Prec);

if m==n
    if isempty(i_perm)
%         output = Prec \ vec;
        output = Prec * vec;
    else
        output = Prec \ vec(i_perm);
    end
else 
    output = Prec .* vec;
end