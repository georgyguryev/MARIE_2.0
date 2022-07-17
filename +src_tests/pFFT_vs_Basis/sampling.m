path_to_ini = './+src_tests/pFFT_vs_Basis/ini_files/';

ini_paths = {'Bcage_Cylinder_0.1R_0.3H_res_5mm',...
            'Bcage_Cylinder_0.145R_0.3H_res_5mm',...
            'Triangular_coil_Cylinder_0.1R_0.3H_res_5mm',...
            'Skyra_coil_RHBM_5mm'}; 
        
task_fnames = {'Bcage_cylinder_5mm_Dir.ini',...
               'Bcage_cylinder_5mm_pFFT.ini',...
               'Bcage_cylinder_5mm_Basis.ini',...
               'Bcage_cylinder_5mm_Dir_decoupled.ini',...
               'Bcage_cylinder_5mm_pFFT_decoupled.ini',...
               'Bcage_cylinder_5mm_Basis_decoupled.ini'};

           
for i_ini_file = 3:length(ini_paths)
    
    ini_path = ini_paths{i_ini_file};
    
    %% loop over tasks, keep a record of iterations and solutions
    
    SOL_marie2 = cell(length(task_fnames),1);
    
    %%
    for i_task = 1:length(task_fnames)
        
        task_fname = task_fnames{i_task};
        
        task_cur_path = fullfile(path_to_ini, ini_path);
        
        % get full file name with a path
        task_specs = fullfile(task_cur_path, task_fname);
        
        % create simulator
        simulator = src_task.EM_simulator();
        
        % construct appropriate task from task_filename with a factory method
        simulator.load_task(task_specs);
        
        Zbc = src_coupling.assemble_coupling_matrix(simulator.task_specs_,...
                                                    simulator.task_runner_.dims,...
                                                    simulator.task_runner_.scatterer,...
                                                    simulator.task_runner_.coil,... 
                                                    simulator.task_runner_.task_settings_.general.Freq_min);
        
        
        % run simulation task
        simulator.run_task();
        
        %% Return results of simulation
        
        SOL_marie2{i_task} = simulator.get_solution();
        
    end
    
    
end
   