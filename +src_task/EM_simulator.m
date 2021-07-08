classdef  (Sealed = true) EM_simulator < handle

    % ================================================================== %
    methods
        
        function obj = EM_simulator()
            
            obj.title.Y = 'Admittance matrix, dB';
            obj.title.Z = 'Impedance matrix, dB';
            obj.title.S = 'Scattering matrix, dB';

        end

        % -------------------------------------------------------------- %
        function [status] = load_task(obj, specs_file)
            % method loads simulation task form the specs file

            % load and parse specs file
            obj.task_loader_ = Task_loader(specs_file);
            
            % init settings
            obj.task_specs_  = obj.task_loader_.get_task_settings();

            % set appropriate task runner
            obj.task_runner_ = obj.task_loader_.create_task_runner();
            
            % call a constructor 
            obj.visualizer_ = src_visualizer.Visualizer();

            status = true;
        end

        % -------------------------------------------------------------- %
        function [status] = run_task(obj)
            % method runs the silulation, defined in task_specs_

            % assert that task_specs are set properly
            assert(~isempty(obj.task_specs_), 'User must define Task first');

            % run task_runner_.run()
            obj.task_runner_.run();

            % update solution
            obj.task_solution_ = obj.task_runner_.get_solution();

            status = true;
        end

        % -------------------------------------------------------------- %
        function  update_coil_position(obj) 
            obj.task_runner_.update_coil_position();
        end
        % -------------------------------------------------------------- %
        
        function show_coil(obj)
            % method visualizes the simulation task

            % call visualization method from task runner
            obj.visualizer_.visualize_coil(obj.task_runner_.coil);
        end
        
        % -------------------------------------------------------------- %
        
        function show_coil_normals(obj)
            
            obj.visualizer_.visualize_normals(obj.task_runner_.coil);
        end
        
        % -------------------------------------------------------------- %
        function show_body(obj)
            
            dom_vie = obj.task_runner_.scatterer.dom_vie;
            id = obj.task_runner_.scatterer.index_vie.S_1d;
            xd = dom_vie.x_tensor;
            yd = dom_vie.y_tensor;
            zd = dom_vie.z_tensor;
            obj.visualizer_.visualize_body(xd,yd,zd,id);
        end

        % -------------------------------------------------------------- %
        function show_coil_and_body(obj)
            
            dom_vie = obj.task_runner_.scatterer.dom_vie;
            id = obj.task_runner_.scatterer.index_vie.S_1d;
            xd = dom_vie.x_tensor;
            yd = dom_vie.y_tensor;
            zd = dom_vie.z_tensor;
            obj.visualizer_.visualize_coil_and_body(xd,yd,zd,id,obj.task_runner_.coil);
        end
        
        % -------------------------------------------------------------- %
        
        function show_frequency_sweep_Sparam(obj, location, port_1, port_2)
            
            % get S parameters
            switch location 
                case 'coil'
                    S_param = obj.task_solution_.Network_param_.coil.S;
                case 'src'
                    S_param = obj.task_solution_.Network_param_.src.S;
            end
            
            % get working frequencies
            [freqs,~] = obj.task_specs_.get_working_frequencies();
            
            % call visualizer function 
            obj.visualizer_.show_frequency_sweep(S_param, port_1, port_2, freqs);
        end
        
        
        % -------------------------------------------------------------- %
        
        
        function show_frequency_sweep_Zparam(obj, location, port_1, port_2)
            
            % get S parameters
            switch location 
                case 'coil'
                    Z_param = obj.task_solution_.Network_param_.coil.Z;
                case 'src'
                    Z_param = obj.task_solution_.Network_param_.src.Z;
            end
            
            % get working frequencies
            [freqs,~] = obj.task_specs_.get_working_frequencies();
            
            % call visualizer function 
            obj.visualizer_.show_frequency_sweep(Z_param, port_1, port_2, freqs);
        end
        
         % -------------------------------------------------------------- %
        
        function show_frequency_sweep_Yparam(obj, location, port_1, port_2)
            
            % get S parameters
            switch location 
                case 'coil'
                    Y_param = obj.task_solution_.Network_param_.coil.Y;
                case 'src'
                    Y_param = obj.task_solution_.Network_param_.src.Y;
            end
            
            % get working frequencies
            [freqs,~] = obj.task_specs_.get_working_frequencies();
            
            % call visualizer function 
            obj.visualizer_.show_frequency_sweep(Y_param, port_1, port_2, freqs);
        end       
        % -------------------------------------------------------------- %
        
        function show_Sparam_matrix(obj, location, freq)
            
            % get S parameters
            switch location 
                case 'coil'
                    S_param = obj.task_solution_.Network_param_.coil.S;
                case 'src'
                    S_param = obj.task_solution_.Network_param_.src.S;
            end
                        
            % get working frequencies
            [freqs,~] = obj.task_specs_.get_working_frequencies();
            
            [~,freq_num] = find(freq == freqs);
            
            % call visualizer if there exists a result for a given
            % frequency
            if ~isempty(freq_num)
                obj.visualizer_.visualize_Network_param_matrix(S_param, freq_num, obj.title.S);
            else
                warning('No data for the specified frequency!');
            end
            
        end
        
        % -------------------------------------------------------------- %
        
        function show_Zparam_matrix(obj, location, freq)
            
            % get Z parameters
            switch location 
                case 'coil'
                    Z_param = obj.task_solution_.Network_param_.coil.Z;
                case 'src'
                    Z_param = obj.task_solution_.Network_param_.src.Z;
            end
                        
            % get working frequencies
            [freqs,~] = obj.task_specs_.get_working_frequencies();
            
            [~,freq_num] = find(freq == freqs);
            
            % call visualizer if there exists a result for a given
            % frequency
            if ~isempty(freq_num)
                obj.visualizer_.visualize_Network_param_matrix(Z_param, freq_num, obj.title.Z);
            else
                warning('No data for the specified frequency!');
            end
            
        end

        % -------------------------------------------------------------- %
        
        function show_Yparam_matrix(obj, location, freq)
            
            % get Y parameters
            switch location 
                case 'coil'
                    Y_param = obj.task_solution_.Network_param_.coil.Y;
                case 'src'
                    Y_param = obj.task_solution_.Network_param_.src.Y;
            end
            
            Y_title = 'Scattering matrix, dB';
                        
            % get working frequencies
            [freqs,~] = obj.task_specs_.get_working_frequencies();
            
            [~,freq_num] = find(freq == freqs);
            
            % call visualizer if there exists a result for a given
            % frequency
            if ~isempty(freq_num)
                obj.visualizer_.visualize_Network_param_matrix(Y_param, freq_num, obj.title.Y);
            else
                warning('No data for the specified frequency!');
            end
            
        end
        

        % -------------------------------------------------------------- %
        
        function show_relative_permittivity(obj, sac, cut)
            
            % get material properties (complex relative permittivity) 
            relative_permittivity = obj.task_runner_.scatterer.prop_vie.Mr;
            
            % check if data is empty
            if isempty(relative_permittivity)
                warning("Can not visualize relative permittivity: Tensor is empty! ");
                return
            end
            
            % call visualizer
            obj.visualizer_.visualize_properties(relative_permittivity, sac, cut);
        end
        
         % -------------------------------------------------------------- %
        
        function show_conductivity(obj, sac, cut)
            
            % get conductivity
            conductivity = obj.task_runner_.scatterer.prop_vie.sigma_e;
            
            % check if data is empty
            if isempty(conductivity)
                warning("Can not visualize conductivity: Tensor is empty! ");
                return
            end
            
            % call visualizer
            obj.visualizer_.visualize_properties(conductivity, sac, cut);
        end

        % -------------------------------------------------------------- %
        
        function show_density(obj, sac, cut)
            
            % get conductivity
            density = obj.task_runner_.scatterer.prop_vie.rho;
            
            % check if data is empty
            if isempty(density)
                warning("Can not visualize density: Tensor is empty! ");
                return
            end
            
            % call visualizer
            obj.visualizer_.visualize_properties(density, sac, cut);
        end

        % -------------------------------------------------------------- %
        
        function h = show_total_electric_field(obj, sac, cut, port, freq)
            
            % get total electric field
            E_tot = obj.task_solution_.get_E_total;
            
            % get working frequencies
            [freqs,~] = obj.task_specs_.get_working_frequencies();
            
            [~,freq_num] = find(freq == freqs);
            
            % check if data is empty
            if isempty(E_tot)
                warning("Can not visualize total electric field: Tensor is empty! ");
                return
            end
            
            % call visualizer if there exists a result for a given
            % frequency
            if ~isempty(freq_num)
                h = obj.visualizer_.visualize_fields(E_tot, sac, cut, port, freq_num);
            else
                warning('No data for the specified frequency!');
                return
            end
        end

        % -------------------------------------------------------------- %
        
        function h = show_total_magnetic_field(obj, sac, cut, port, freq)
            
            % get total magnetic field
            H_tot = obj.task_solution_.get_H_total;
            
            % get working frequencies
            [freqs,~] = obj.task_specs_.get_working_frequencies();
            
            [~,freq_num] = find(freq == freqs);
            
            % check if data is empty
            if isempty(H_tot)
                warning("Can not visualize total magnetic field: Tensor is empty! ");
                return
            end
            
            % call visualizer if there exists a result for a given
            % frequency
            if ~isempty(freq_num)
               h = obj.visualizer_.visualize_fields(H_tot, sac, cut, port, freq_num);
            else
                warning('No data for the specified frequency!');
                return
            end
        end

        % -------------------------------------------------------------- %
        
        function h = show_polarization_currents(obj, sac, cut, port, freq)
            
            % get polarization currents field
            Jb = obj.task_solution_.get_plrz_currents_tensor();
            
            % get working frequencies
            [freqs,~] = obj.task_specs_.get_working_frequencies();
            
            [~,freq_num] = find(freq == freqs);
            
            % check if data is empty
            if isempty(Jb)
                warning("Can not visualize polarization currents: Tensor is empty! ");
                return
            end
            
            % call visualizer if there exists a result for a given
            % frequency
            if ~isempty(freq_num)
                h = obj.visualizer_.visualize_fields(Jb, sac, cut, port, freq_num);
            else
                warning('No data for the specified frequency!');
                return
            end
        end

        % -------------------------------------------------------------- %
        
        function h = show_surface_currents(obj, port, freq)
            
            % get surface currents
            Jc = obj.task_solution_.get_surface_currents;
            
            % get working frequencies
            [freqs,~] = obj.task_specs_.get_working_frequencies();
            
            [~,freq_num] = find(freq == freqs);
            
            % check if data is empty
            if isempty(Jc)
                warning("Can not visualize surface currents: Vector is empty! ");
                return
            end
            
            % call visualizer if there exists a result for a given
            % frequency
            if ~isempty(freq_num)
                h = obj.visualizer_.visualize_currents(obj.task_runner_.coil, Jc, port, freq_num);
            else
                warning('No data for the specified frequency!');
                return
            end
        end

        % -------------------------------------------------------------- %
        
        function show_absorbed_power(obj, sac, cut, port, freq)
            
            % get local absorbed power
            P_abs = obj.task_solution_.get_P_abs;
            
            % get working frequencies
            [freqs,~] = obj.task_specs_.get_working_frequencies();
            
            [~,freq_num] = find(freq == freqs);
            
            % check if data is empty
            if isempty(P_abs)
                warning("Can not visualize local absorbed power: Tensor is empty! ");
                return
            end
            
            % call visualizer if there exists a result for a given
            % frequency
            if ~isempty(freq_num)
                obj.visualizer_.visualize_fields(P_abs, sac, cut, port, freq_num);
            else
                warning('No data for the specified frequency!');
                return
            end
        end

        % -------------------------------------------------------------- %
        
        function show_B1_plus(obj, sac, cut, port, freq)
            
            % get b1 plus
            B1_p = obj.task_solution_.get_B1_p;
            
            % get working frequencies
            [freqs,~] = obj.task_specs_.get_working_frequencies();
            
            [~,freq_num] = find(freq == freqs);
            
            % check if data is empty
            if isempty(B1_p)
                warning("Can not visualize B1 plus: Tensor is empty! ");
                return
            end
            
            % call visualizer if there exists a result for a given
            % frequency
            if ~isempty(freq_num)
                obj.visualizer_.visualize_fields(B1_p, sac, cut, port, freq_num);
            else
                warning('No data for the specified frequency!');
                return
            end
        end

        % -------------------------------------------------------------- %
        
        function show_B1_minus(obj, sac, cut, port, freq)
            
            % get b1 minus
            B1_m = obj.task_solution_.get_B1_p;
            
            % get working frequencies
            [freqs,~] = obj.task_specs_.get_working_frequencies();
            
            [~,freq_num] = find(freq == freqs);
            
            % check if data is empty
            if isempty(B1_m)
                warning("Can not visualize B1 minus: Tensor is empty! ");
                return
            end
            
            % call visualizer if there exists a result for a given
            % frequency
            if ~isempty(freq_num)
                obj.visualizer_.visualize_fields(B1_m, sac, cut, port, freq_num);
            else
                warning('No data for the specified frequency!');
                return
            end
        end


        % -------------------------------------------------------------- %

        function [task_solution] = get_solution(obj)
            % function returns the task_solution_ object to a user

            task_solution = obj.task_solution_;
        end

    end

    % =================================================================== %

    properties %(Access = private)

        % object parses specs file and constructs task_specs_ & task_runner_
        task_loader_ = [];

        % object stores task specifications and settings
        task_specs_  = [];

        % solver engine
        task_runner_ = [];

        % object stores the solutions
        task_solution_  = [];
                
        % object visualizes 
        visualizer_ = [];
        
        title = struct();

        
    end
    

end
