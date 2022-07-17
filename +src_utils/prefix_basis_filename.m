function pref_File = prefix_basis_filename(task_settings)

file = task_settings.basis.Filename;
solver_mode = lower(task_settings.vsie.Solver_mode);

% label basis type with prefix
switch solver_mode
    case 'explicit'
        pref_File = strcat('Coupled_basis_', file);
    case 'coil_implicit'
        pref_File = strcat('Perturbation_basis_', file);
    case 'tissue_implicit'
        pref_File = strcat('Coupled_basis_', file);
end

