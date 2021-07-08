
% Clear variables (created before) in the memory namespace 
% clear;

% Alias for object files extension
if ispc
    str_obj_ext = '.obj';
elseif isunix || ismac
    str_obj_ext = '.o';
end

% in mac it adds \t\n to strjoin() so it does not create a proper string!

%% Compilation
% define compilation command with proper parameters
str_mex_command_with_opts = 'mex "-I../include" -c -outdir "../" ';

cd './utils';
str_src_utils = 'GL_1D.cpp';
eval([str_mex_command_with_opts str_src_utils]);

cd '../WS_EA'
str_src_ea = '';
str_src_ea = strcat(str_src_ea, ' A_functions_ws_ea.cpp');
str_src_ea = strcat(str_src_ea, ' N_functions_ws_ea.cpp');
str_src_ea = strcat(str_src_ea, ' PSI_limits_ws_ea.cpp');
str_src_ea = strcat(str_src_ea, ' Simplex_ws_ea.cpp');
str_src_ea = strcat(str_src_ea, ' THETA_limits_ws_ea.cpp');
str_src_ea = strcat(str_src_ea, ' quadric_ws_ea.cpp');
eval([str_mex_command_with_opts str_src_ea]);

cd '../WS_ST'
str_src_st = '';
str_src_st = strcat(str_src_st, ' a_functions_ws_st.cpp');
str_src_st = strcat(str_src_st, ' lambda_limits_ws_st.cpp');
str_src_st = strcat(str_src_st, ' n_functions_ws_st.cpp');
str_src_st = strcat(str_src_st, ' psi_limits_ws_st.cpp');
str_src_st = strcat(str_src_st, ' quadric_ws_st.cpp');
str_src_st = strcat(str_src_st, ' subtriangles_ws_st.cpp');
str_src_st = strcat(str_src_st, ' u_limits_ws_st.cpp');
eval([str_mex_command_with_opts str_src_st]);

cd '../WS_VA'
str_src_va = '';
str_src_va = strcat(str_src_va, ' rho_limits_ws_va.cpp');
str_src_va = strcat(str_src_va, ' theta_p_limits_ws_va.cpp');
str_src_va = strcat(str_src_va, ' theta_q_limits_ws_va.cpp');
str_src_va = strcat(str_src_va, ' quadric_ws_va.cpp');
eval([str_mex_command_with_opts str_src_va]);

cd '../'

%%  Creating a string with obj files to be linked in mex

str_obj = '';
str_obj = strcat(str_obj, strcat(' GL_1D', str_obj_ext));

str_obj = strcat(str_obj, strcat(' A_functions_ws_ea', str_obj_ext));
str_obj = strcat(str_obj, strcat(' N_functions_ws_ea', str_obj_ext));
str_obj = strcat(str_obj, strcat(' PSI_limits_ws_ea', str_obj_ext));
str_obj = strcat(str_obj, strcat(' Simplex_ws_ea', str_obj_ext));
str_obj = strcat(str_obj, strcat(' THETA_limits_ws_ea'), str_obj_ext);
str_obj = strcat(str_obj, strcat(' quadric_ws_ea', str_obj_ext));

str_obj = strcat(str_obj, strcat(' a_functions_ws_st', str_obj_ext));
str_obj = strcat(str_obj, strcat(' lambda_limits_ws_st', str_obj_ext));
str_obj = strcat(str_obj, strcat(' n_functions_ws_st', str_obj_ext));
str_obj = strcat(str_obj, strcat(' psi_limits_ws_st', str_obj_ext));
str_obj = strcat(str_obj, strcat(' quadric_ws_st', str_obj_ext));
str_obj = strcat(str_obj, strcat(' subtriangles_ws_st', str_obj_ext));
str_obj = strcat(str_obj, strcat(' u_limits_ws_st', str_obj_ext));

str_obj = strcat(str_obj, strcat(' rho_limits_ws_va', str_obj_ext));
str_obj = strcat(str_obj, strcat(' theta_p_limits_ws_va', str_obj_ext));
str_obj = strcat(str_obj, strcat(' theta_q_limits_ws_va', str_obj_ext));
str_obj = strcat(str_obj, strcat(' quadric_ws_va', str_obj_ext));

%% EA mex
str_mex_ea = '';
str_mex_ea = strcat(str_mex_ea, 'mex -output solve_ea ');
str_mex_ea = strcat(str_mex_ea, str_obj);
str_mex_ea = strcat(str_mex_ea, ' solve_ea_mex.cpp create_EA.cpp Kernels.cpp');
eval(str_mex_ea);

%% ST mex
str_mex_st = '';
str_mex_st = strcat(str_mex_st, 'mex -output solve_st ');
str_mex_st = strcat(str_mex_st, str_obj);
str_mex_st = strcat(str_mex_st, ' solve_st_mex.cpp create_ST.cpp Kernels.cpp');
eval(str_mex_st);

%% ST mex
str_mex_va = '';
str_mex_va = strcat(str_mex_va, 'mex -output solve_va ');
str_mex_va = strcat(str_mex_va, str_obj);
str_mex_va = strcat(str_mex_va, ' solve_va_mex.cpp create_VA.cpp Kernels.cpp');
eval(str_mex_va);

% Remove object files after compilation
delete(strcat('*',str_obj_ext));

%% End of the script



% list = strjoin(cellstr(ls('*.obj')));
% eval(['mex -output solve_ea ' list ' solve_ea_mex.cpp create_EA.cpp Kernels.cpp'])
% eval(['mex -output solve_va ' list ' solve_va_mex.cpp create_VA.cpp Kernels.cpp'])
% eval(['mex -output solve_st ' list ' solve_st_mex.cpp create_ST.cpp Kernels.cpp'])


% aa = {'aa1','ab2','ac3'...
%       'ad4', 'af5'};
% bb = {'bb1', 'bc2', 'bc3'};
% aa_bb = horzcat( aa , bb);
% % nem = len(aa_bb);
% for i = len(aa_bb)
%    disp(aa_bb(i)) 
% end

%% End of the file
