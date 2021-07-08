classdef Task_settings < handle

    methods

        % define constructor
        function obj = Task_settings(task_settings)

            obj.general = task_settings.general;
            obj.sie     = task_settings.sie;
            obj.vie     = task_settings.vie;
            obj.vsie    = task_settings.vsie;
            obj.basis   = task_settings.basis;
            
            % check input settings
            obj.validate_task_settings_();

        end

        % define interface methods
        function task_type = task_type(obj)
            task_type = obj.general.Task;
        end

        % get frequency info
        function [freqs, N_freqs] = get_working_frequencies(obj)

            % get number of working frequencies
            N_freqs  = obj.general.N_freqs;
            freq_min = obj.general.Freq_min;
            freq_max = obj.general.Freq_max;

            if N_freqs == 1
                freqs = freq_min;
            else
                df    = (freq_max - freq_min) / (N_freqs - 1);
                freqs = freq_min:df:freq_max;
            end

        end


        function [coil_msh, coil_txt] = get_coil_filenames(obj)
            % path to files
            coil_path = obj.sie.Coil_path;

            % full .msh and .txt file names
            coil_msh = fullfile(coil_path, obj.sie.Mesh_file);
            coil_txt = fullfile(coil_path, obj.sie.Ports_file);

        end

        function body_fname = get_body_filename(obj)
            body_path  = obj.vie.Body_path;
            body_fname = fullfile(body_path, obj.vie.Body_file);
        end

    end


    methods (Access = private)

        function [valid] = validate_task_settings_(obj)

            % validate frequencies
            assert(obj.frequencies_are_valid_(), 'Sanity check failed for frequency properties');

            valid = true;

        end

        function [valid] = frequencies_are_valid_(obj)

            % check if freq_min is provided
            assert(isnumeric(obj.general.Freq_min) && obj.general.Freq_min > 0, 'Variable freq_min was not provided');

            % check if N_freqs is specified properly
            if (~isnumeric(obj.general.N_freqs) || obj.general.N_freqs < 0 || ...
                      floor(obj.general.N_freqs) ~= obj.general.N_freqs)

                obj.general.N_freqs = 1;

            end

            % if freq_max is provided check that freq_min <= freq_max;
            if ~isempty(obj.general.Freq_max) && isnumeric(obj.general.Freq_max)

                assert(obj.general.Freq_max >= obj.general.Freq_min);

            else

                obj.general.Freq_max = obj.general.Freq_min;
                obj.general.N_freqs  = 1;

            end

            valid = true;

        end

    end

    % define list of properties
    properties (SetAccess = immutable)

        general
        sie
        vie
        vsie
        basis
        output

    end
end
