function [op_out] = getOPERATORS(domain,freq,op,l,tucker,tucker_tol)
    % function that retrieves the operators N and K for JM-VIE
    % 
    % INPUT
    % r         3-D Domain
    % f         working frequency
    % 
    % op        flag for choosing the operator
    %                N: only operator N
    %                K: only operator K 
    % l         flag for choosing the basis functions
    %                1: for PWC
    %                4: for PWL
    % tucker    flag for compressed operators using the HOSVD
    % 
    % OUTPUT
    % op_out    discrete operator
    % _________________________________________________________________________
    
    
    
    % -------------------------------------------------------------------------
    % Parallelization
    % -------------------------------------------------------------------------
    try
        if verLessThan('matlab', '8.3.0')
            isOpen = (matlabpool('size') > 0); %#ok<*DPOOL>
            if (~isOpen)
                mycluster = parcluster;
                matlabpool('local',mycluster.NumWorkers);
            end
        end
    catch
        fprintf(1, '\n\n WARNING: matlabpool could not open, running on sequential');
        fprintf(1, '\n            computational performance may be compromised');
        fprintf(1, '\n            consider checking your MATLAB parallel capabilities');
        fprintf(1, '\n\n');
        pause(1);
    end

    % -------------------------------------------------------------------------
    % Initialization
    % -------------------------------------------------------------------------    
    
    % set default tucker tolerance if none provided
    if nargin < 6
        tucker_tol = 1e-16;
    end
    
    
    fid = 1;
    emc = src_utils.EM_utils(freq);
    res = domain.res;
    dom_dims = size(domain.x_tensor);
    r = domain.r;
    
    % Pick PWC/PWL flag
    if l == 4
        pwx      = 10;
        pwx_flag = true;
    else
        pwx      = 1;
        pwx_flag = false;
    end
    
    % Pick FFT+Circ flag
    if l == 1 && ~tucker
        flag_fft_circ = 1;
    elseif l == 4 && ~tucker
        flag_fft_circ = 2;
    elseif l == 1 && tucker
        flag_fft_circ = 3;
    elseif l == 4 && tucker
        flag_fft_circ = 4;
    else
        warning('Unavailable FFT+Circulan form');
        return
    end
    
    
%% Assembly of VIE operators with Toeplitz structure

    % if vie assembly was sucessfully compiled - use mexed version
    if (3 == exist('vie_assembly', 'file'))
    tic
        % run mex-ed VIE assembly, reshape output operator to be consistent
        % with the matlab ordering
        op_out = vie_assembly(op, freq, res, dom_dims, pwx_flag);
        op_out = permute(op_out, [3 4 5 1 2]);
    toc
    else
    tic    
        % Retrieve Operator
        if strcmp(op ,'N')
            
            op_out = assembly_N_(r,emc,res,pwx,l);
            
        elseif strcmp(op, 'K')
            
            op_out = assembly_K_(r,emc,res,pwx,l);
            
        else
            warning('Unavailable Operator');
            return
        end
    toc
    end

    %% Transform Toeplitz to Circulant VIE operator and compute its FFT
    
    switch flag_fft_circ
        
        case 1
            
            op_out = assembly_fft_circ_pwc_(op_out);
            
        case 2
            
            op_out = assembly_fft_circ_pwl_(op_out);
            
        case 3
            
            op_out = assembly_fft_circ_tucker_pwc_(op_out,tucker_tol);
            
        case 4
            
            op_out = assembly_fft_circ_tucker_pwl_(op_out,tucker_tol);
    end
    
    infocir = whos('op_out');
    fprintf(fid, '\n Frequency:           %3.3f MHz', freq/1e6);
    fprintf(fid, '\n Domain Dimensions:   %dx%dx%d', dom_dims);
    fprintf(fid, '\n Memory space:        %.6f MB', infocir.bytes/(1024*1024));
    fprintf(fid, '\n ----------------------------------------------------------\n ');

end