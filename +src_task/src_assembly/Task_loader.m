classdef Task_loader < handle
% class src_task.Task_loader parses .txt Task file and creates the
% appropriate TaskRunner based on the splecified simulation type
% (either SIE, VIE or VSIE)

    % =================================================================== %

    properties (Access = private)
        specs_file
        task_settings
    end


    % =================================================================== %

    methods

        % function defines constructor
        function obj = Task_loader(specs_file)

            obj.specs_file      = specs_file;
            obj.task_settings = obj.parse_task_spec_();
        end

        % --------------------------------------------------------------- %
        function [task_runner] = create_task_runner(obj)
            % method creates a specific Task object

            switch obj.task_settings.task_type()
                case 'SIE'
                    task_runner = TaskRunner_SIE(obj.task_settings);
                case 'VIE'
                    task_runner = TaskRunner_VIE(obj.task_settings);
                case 'VSIE'
                    task_runner = TaskRunner_VSIE(obj.task_settings);
            end
        end

        % --------------------------------------------------------------- %
        function [task_settings] = get_task_settings(obj)
            % method returns task settings

            task_settings = obj.task_settings;
        end

        % --------------------------------------------------------------- %
        function [specs_fname] = get_specs_filename(obj)
            % method returns specs file name
            specs_fname = obj.specs_file;
        end

    end

    % =================================================================== %
    % define private methods
    methods (Access = private)

        function [task_settings] = parse_task_spec_(obj)

            % instantiate ini file parser
            ini_parser = src_utils.IniConfig();

            % parse and preprocess file
            ini_parser.ReadFile(obj.specs_file);
            setting_types = ini_parser.GetSections();

            % combine keys and values; form settings structure
            for i = 1:length(setting_types)

                setting_type = setting_types{i};

                [keys, ~] = ini_parser.GetKeys(setting_type);
                values    = ini_parser.GetValues(setting_type, keys);

                % check if some values represent numbers
                dbl_values = str2double(values);

                % find idx of non-nan values
                idx = find(~isnan(dbl_values));

                if nnz(dbl_values) > 0

                    cell_dbl_values = num2cell(dbl_values,length(values));

                    values(idx) = cell_dbl_values(idx);
                end

                field_name = lower(setting_type(2:end-1));

                task_settings.(field_name) = cell2struct(values, keys);
            end

            % create an instance of Task_settings
            task_settings = Task_settings(task_settings);

        end
    end
end
