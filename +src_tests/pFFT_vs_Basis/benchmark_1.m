 % MARIE 2.0 example

% addpath(genpath('../'));

% task_filename

task_path  = './+src_tests/pFFT_vs_Basis/ini_files/';

ini_paths = {'Bcage_Cylinder_0.1R_0.3H_res_5mm',...
            'Bcage_Cylinder_0.145R_0.3H_res_5mm',...
            'Bcage_RHBM_res_5mm',...
            'Triangular_coil_Cylinder_0.1R_0.3H_res_5mm',...
            'Triangular_coil_RHBM_res_5mm',...
            'Skyra_coil_RHBM_5mm'}; 
           
for i_ini_file = 2:3%length(ini_paths)
    
%     task_fnames = {'Bcage_cylinder_5mm_pFFT.ini',...
%                'Bcage_cylinder_5mm_pFFT_decoupled.ini'};
%            
           
%     task_fnames = {'Bcage_cylinder_5mm_Dir.ini',...
%                'Bcage_cylinder_5mm_pFFT.ini',...
%                'Bcage_cylinder_5mm_Basis.ini',...
%                'Bcage_cylinder_5mm_Dir_decoupled.ini',...
%                'Bcage_cylinder_5mm_pFFT_decoupled.ini',...
%                'Bcage_cylinder_5mm_Basis_decoupled.ini'};       
    

%     task_fnames = {'Bcage_cylinder_5mm_Basis.ini',...
%                    'Bcage_cylinder_5mm_Basis_decoupled.ini'};

    ini_path = ini_paths{i_ini_file};
    
    %% loop over tasks, keep a record of iterations and solutions
%     
%     SOL_marie2 = cell(length(task_fnames),1);
%     
%     %%
%     for i_task =  1:length(task_fnames)
%         
%         task_fname = task_fnames{i_task};
%         
%         task_cur_path = fullfile(task_path, ini_path);
%         
%         % get full file name with a path
%         task_specs = fullfile(task_cur_path, task_fname);
%         
%         % create simulator
%         simulator = src_task.EM_simulator();
%         
%         % construct appropriate task from task_filename with a factory method
%         simulator.load_task(task_specs);
%         
%         % run simulation task
%         simulator.run_task();
%         
%         %% Return results of simulation
%         
%         SOL_marie2{i_task} = simulator.get_solution();
%         
%         clear simulator;
%         
%     end
%     %%
%     sv_filename = sprintf('./+src_tests/pFFT_vs_Basis/PWC/sweep/1e-3/%s_GMRES_DR.mat', ini_path);
%     
%     % save results
%     save(sv_filename, 'SOL_marie2', '-v7.3');
%     
%     clear SOL_marie2;
    
    %% loop over tasks, keep a record of iterations and solutions
    

%     task_fnames = {'Bcage_cylinder_5mm_Dir_TFQMR.ini',...
%         'Bcage_cylinder_5mm_pFFT_TFQMR.ini',...
%         'Bcage_cylinder_5mm_Basis_TFQMR.ini',...
%         'Bcage_cylinder_5mm_Dir_decoupled_TFQMR.ini',...
%         'Bcage_cylinder_5mm_pFFT_decoupled_TFQMR.ini',...
%         'Bcage_cylinder_5mm_Basis_decoupled_TFQMR.ini'};
    
    task_fnames = {'Bcage_cylinder_5mm_Dir_decoupled.ini',...
                   'Bcage_cylinder_5mm_pFFT_decoupled_TFQMR.ini'};
    
    
%     task_fnames = {'Bcage_cylinder_5mm_Basis_TFQMR.ini',...
%                    'Bcage_cylinder_5mm_Basis_decoupled_TFQMR.ini'};

    SOL_marie2 = cell(length(task_fnames),1);
    
    %%
    for i_task = 1:length(task_fnames)
        
        task_fname = task_fnames{i_task};
        
        task_cur_path = fullfile(task_path, ini_path);
        
        % get full file name with a path
        task_specs = fullfile(task_cur_path, task_fname);
        
        % create simulator
        simulator = src_task.EM_simulator();
        
        % construct appropriate task from task_filename with a factory method
        simulator.load_task(task_specs);
        
        % run simulation task
        simulator.run_task();
        
        %% Return results of simulation
        
        SOL_marie2{i_task} = simulator.get_solution();
        
    end
    
    sv_filename = sprintf('./+src_tests/pFFT_vs_Basis/PWC/sweep/1e-3/%s_fine.mat', ini_path);
    
    % save results
    save(sv_filename, 'SOL_marie2', '-v7.3');
       
    
    clear SOL_marie2;
end
 
