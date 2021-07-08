function [ini] = create_ini_vsie(rhbm,coil,inp)

    rhbm_name = rhbm.name;
    coil_name = coil.name;
    ini = strcat('./data/inputs/user/',rhbm_name,'_',coil_name,'.ini');
    fid = fopen(ini, 'W');
    
    bodyname = strcat(rhbm_name,'.mat');
    coilnamemat = strcat('user/',coil_name,'.mat');
    coilnamemsh = strcat('msh_files/',coil_name,'.msh');
    coilnametxt = strcat('lumped_elements_files/',extractBefore(coil_name,'_'),'.txt');
    solutionname = strcat('./data/solutions/SOLUTION_',rhbm_name,'_',coil_name,'.mat');

    fprintf(fid, '[General] \n');
    fprintf(fid, 'Task = VSIE ;\n');
    fprintf(fid, 'Precision = DOUBLE ;\n');  %Could be float in the future
    fprintf(fid, 'Freq_min = %e ;\n',inp.freq.freq_min);
    fprintf(fid, 'Freq_max = %e ;\n',inp.freq.freq_max);
    fprintf(fid, 'N_freqs = %d ;\n',inp.freq.range);
    
    if strcmp(inp.PU,'CPU')
        temp = 0;
    elseif strcmp(inp.PU,'GPU')
        temp = 1;
    end
    
    fprintf(fid, 'GPU_flag = %d ;\n',temp);
    fprintf(fid,'\n');
    fprintf(fid, '[VIE] \n');
    fprintf(fid, 'BodyModel = %s ;\n',bodyname);
    
    if strcmp(inp.PWX,'PWC')
        temp = 0;
    elseif strcmp(inp.PWX,'PWL')
        temp = 1;
    end
    
    fprintf(fid, 'PWX = %d \n',temp);
    if strcmp(inp.VIE_solver.name,'G')
        fprintf(fid, 'Iterative_method = GMRES ;\n');
        fprintf(fid, 'SolTolerance = %e ;\n', inp.VIE_solver.tol);
        fprintf(fid, 'MaxIter = %d ;\n', inp.VIE_solver.outer_iter);
        fprintf(fid, 'Restart = %d ;\n', inp.VIE_solver.inner_iter);
    elseif  strcmp(inp.VIE_solver.name,'B')
        fprintf(fid, 'Iterative_method = BICGSTAB ;\n');
        fprintf(fid, 'SolTolerance = %e ;\n', inp.VIE_solver.tol);
        fprintf(fid, 'MaxIter = %d ;\n', inp.VIE_solver.outer_iter);
    end
    fprintf(fid, 'SolverHybridModel = %s ;\n',inp.PU);
    fprintf(fid, 'Tucker = %d ;\n',inp.VIE_tucker.flag);
    fprintf(fid, 'TuckerTolerance = %e ;\n',inp.VIE_tucker.tol);
    fprintf(fid, 'Np_quad_far = %d ;\n',inp.Np_quad_far);
    fprintf(fid, 'Np_quad_medium = %d ;\n',inp.Np_quad_medium);
    fprintf(fid, 'Np_quad_near = %d ;\n',inp.Np_quad_near);
    fprintf(fid, '\n');
    
    fprintf(fid, '[SIE] ;\n');
    if strcmp(inp.SIE_solver.name,'G')
        fprintf(fid, 'InvertMatrix = GMRES ;\n');
        fprintf(fid, 'SolTolerance = %e ;\n', inp.VIE_solver.tol);
        fprintf(fid, 'MaxIter = %d ;\n', inp.VIE_solver.outer_iter);
        fprintf(fid, 'Restart = %d ;\n', inp.VIE_solver.inner_iter);
    elseif  strcmp(inp.SIE_solver.name,'B')
        fprintf(fid, 'InvertMatrix = BICGSTAB ;\n');
        fprintf(fid, 'SolTolerance = %e ;\n', inp.VIE_solver.tol);
        fprintf(fid, 'MaxIter = %d ;\n', inp.VIE_solver.outer_iter);
    elseif strcmp(inp.SIE_solver.name,'L')
        fprintf(fid, 'InvertMatrix = LU ;\n');
    end
    fprintf(fid, 'CoilMeshFile = %s ;\n',coilnamemsh);
    fprintf(fid, 'CoilPortsFile = %s ;\n',coilnametxt);
    fprintf(fid, 'CoilModel = %s ;\n',coilnamemat);
    fprintf(fid, 'Coil_Path = %s ;\n','./data/coils');
    fprintf(fid, '\n');
    
    if strcmp(inp.VSIE_coupling,'Coupled')
        temp_vsie = 'Unnested';
    else
        temp_vsie = 'Nested';
    end
    
    fprintf(fid, '[VSIE] \n');
    fprintf(fid, 'VSIE_Solver = %s ;\n',temp_vsie);
    fprintf(fid, 'Coupling_Mode = Implicit ;\n');
    fprintf(fid, 'N_expansion_voxels = %d ;\n',inp.VSIE_expansion_voxels);
    fprintf(fid, 'N_near_voxels = 1 ;\n');
    fprintf(fid, 'Proj_assembly_mode = Sphere ;\n');
    fprintf(fid, 'Near_boundary_width = 1 ;\n');
    fprintf(fid, 'Precond_type = Symmetric ;\n');
    fprintf(fid, 'Iterative_method = GMRES_DR ;\n');
    fprintf(fid, 'Tolerance = %e ;\n', inp.VSIE_solver.tol);
    fprintf(fid, 'MaxIter = %d ;\n', inp.VSIE_solver.outer_iter);
    fprintf(fid, 'Restart = %d ;\n', inp.VSIE_solver.inner_iter);
    fprintf(fid, 'Ritz = %d \n', inp.VSIE_solver.inner_iter);
    fprintf(fid, 'Np_quad_coup_sie = %d ;\n', inp.Np_quad_coup_sie);
    fprintf(fid, 'Np_quad_coup_vie = %d ;\n', inp.Np_quad_coup_vie);
    fprintf(fid, '\n');
    
    fprintf(fid, '[Output] \n');
    fprintf(fid, 'Solution = %s ;\n',solutionname);
    
    fclose(fid);

end

