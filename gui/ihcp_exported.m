classdef ihcp_exported < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        ihcpUIFigure                   matlab.ui.Figure
        FileMenu                       matlab.ui.container.Menu
        OpenMenu                       matlab.ui.container.Menu
        SetcurrentvaluesasdefaultMenu  matlab.ui.container.Menu
        ResettodefaultMenu             matlab.ui.container.Menu
        MattLenahan2023Menu            matlab.ui.container.Menu
        ThermocouplepositionsPanel     matlab.ui.container.Panel
        T_Cu3pos                       matlab.ui.control.NumericEditField
        T_Cu3EditFieldLabel            matlab.ui.control.Label
        T_Inco2pos                     matlab.ui.control.NumericEditField
        T_Inco2EditFieldLabel          matlab.ui.control.Label
        T_Inco1pos                     matlab.ui.control.NumericEditField
        T_Inco1EditFieldLabel          matlab.ui.control.Label
        T_Cu2pos                       matlab.ui.control.NumericEditField
        T_Cu2EditFieldLabel            matlab.ui.control.Label
        GeometryPanel                  matlab.ui.container.Panel
        lengthCu2                      matlab.ui.control.NumericEditField
        LengthCu2EditFieldLabel        matlab.ui.control.Label
        lengthInco2                    matlab.ui.control.NumericEditField
        LengthInco2EditFieldLabel      matlab.ui.control.Label
        lengthH25                      matlab.ui.control.NumericEditField
        LengthH25EditFieldLabel        matlab.ui.control.Label
        lengthInco1                    matlab.ui.control.NumericEditField
        LengthInco1EditFieldLabel      matlab.ui.control.Label
        lengthCu1                      matlab.ui.control.NumericEditField
        LengthCu1EditFieldLabel        matlab.ui.control.Label
        BrowseButton                   matlab.ui.control.Button
        ParametersPanel                matlab.ui.container.Panel
        plotCheck                      matlab.ui.control.CheckBox
        dt                             matlab.ui.control.NumericEditField
        dxLabel_2                      matlab.ui.control.Label
        dx                             matlab.ui.control.NumericEditField
        dxLabel                        matlab.ui.control.Label
        epsilonEditField               matlab.ui.control.NumericEditField
        epsilonEditFieldLabel          matlab.ui.control.Label
        maxIterSpinner                 matlab.ui.control.Spinner
        maxIterSpinnerLabel            matlab.ui.control.Label
        RTOLerrorEditField             matlab.ui.control.NumericEditField
        RTOLerrorEditFieldLabel        matlab.ui.control.Label
        RTOLhEditField                 matlab.ui.control.NumericEditField
        RTOLhEditFieldLabel            matlab.ui.control.Label
        hInitialEditField              matlab.ui.control.NumericEditField
        hInitialEditFieldLabel         matlab.ui.control.Label
        rSpinner                       matlab.ui.control.Spinner
        rSpinnerLabel                  matlab.ui.control.Label
        RunButton                      matlab.ui.control.Button
        selectedfile                   matlab.ui.control.EditField
        SelectedfileEditFieldLabel     matlab.ui.control.Label
    end

    
    methods (Access = private)

        function getfile(app)
            % fix for lost focus
            f = figure('Position', [-100 -100 0 0]); % off-screen
            % file picker
            [file, path] = uigetfile({'*.*'; '*.csv'; '*.dat'; ...
                '*.txt'}, "Open...");
            delete(f);
            % if file selected i.e. user doesn't cancel
            if file ~= 0
                app.selectedfile.Value = fullfile(path, file);
                app.RunButton.Enable = "on";
                % prepare
                prepdata(fullfile(path, file))
            end
        end

        function savestate(app)

            state.T_Cu2pos.Value = app.T_Cu2pos.Value;
            state.T_Inco1pos.Value = app.T_Inco1pos.Value;
            state.T_Inco2pos.Value = app.T_Inco2pos.Value;
            state.T_Cu3pos.Value = app.T_Cu3pos.Value;

            state.lengthCu1.Value = app.lengthCu1.Value;
            state.lengthInco1.Value = app.lengthInco1.Value;
            state.lengthH25.Value = app.lengthH25.Value;
            state.lengthInco2.Value = app.lengthInco2.Value;
            state.lengthCu2.Value = app.lengthCu2.Value;

            state.dx.Value = app.dx.Value;
            state.dt.Value = app.dt.Value;

            state.rSpinner.Value = app.rSpinner.Value;
            state.epsilonEditField.Value = app.epsilonEditField.Value;
            state.hInitialEditField.Value = app.hInitialEditField.Value;
            state.RTOLhEditField.Value = app.RTOLhEditField.Value;
            state.RTOLerrorEditField.Value = app.RTOLerrorEditField.Value;
            state.maxIterSpinner.Value = app.maxIterSpinner.Value;

            state.plotCheck.Value = app.plotCheck.Value;

            save("defaultState.mat", "state")
        end

        function loadstate(app)
            load("defaultState.mat", "state")

            app.T_Cu2pos.Value = state.T_Cu2pos.Value;
            app.T_Inco1pos.Value = state.T_Inco1pos.Value;
            app.T_Inco2pos.Value = state.T_Inco2pos.Value;
            app.T_Cu3pos.Value = state.T_Cu3pos.Value;

            app.lengthCu1.Value = state.lengthCu1.Value;
            app.lengthInco1.Value = state.lengthInco1.Value;
            app.lengthH25.Value = state.lengthH25.Value;
            app.lengthInco2.Value = state.lengthInco2.Value;
            app.lengthCu2.Value = state.lengthCu2.Value;

            app.dt.Value = state.dt.Value;
            app.dx.Value = state.dx.Value;

            app.rSpinner.Value = state.rSpinner.Value;
            app.epsilonEditField.Value = state.epsilonEditField.Value;
            app.hInitialEditField.Value = state.hInitialEditField.Value;
            app.RTOLhEditField.Value = state.RTOLhEditField.Value;
            app.RTOLerrorEditField.Value = state.RTOLerrorEditField.Value;
            app.maxIterSpinner.Value = state.maxIterSpinner.Value;

            app.plotCheck.Value = state.plotCheck.Value;
        end

    end


    % Callbacks that handle component events
    methods (Access = private)

        % Code that executes after component creation
        function startupFcn(app)
            % button icon
            app.RunButton.Icon = fullfile(matlabroot,'ui','install',...
                'installer_login','images','toolstrip','Run_24.png');
            try loadstate(app)
            catch
                uialert(app.ihcpUIFigure, "No defaults found.", "Warning", Icon="warning")
            end
        end

        % Menu selected function: OpenMenu
        function OpenMenuSelected(app, event)
            getfile(app)
        end

        % Button pushed function: RunButton
        function RunButtonPushed(app, event)
            % gather all entered data
            % filename
            filename = app.selectedfile.Value;
           
            % geometry
            geometry = cell(6, 1);
            geometry{1} = [app.lengthCu1.Value app.lengthInco1.Value ...
                app.lengthH25.Value app.lengthInco2.Value ...
                app.lengthCu2.Value];
            geometry{2} = sum(geometry{1}, "all");
            geometry{3} = app.dx.Value;
            geometry{4} = round(geometry{1}./geometry{3});
            if mod(geometry{1}, geometry{3})
                uialert(app.ihcpUIFigure, "Length(s) not divisible by dx", "Error", Icon="error")
                return
            end
            geometry{5} = sum(geometry{4});
            positions = [app.T_Cu2pos.Value app.T_Inco1pos.Value ...
                app.T_Inco2pos.Value app.T_Cu3pos.Value];
            geometry{6} = round(.5+(positions)/geometry{3});
            if any(mod(.5+(positions), geometry{3}))
                uialert(app.ihcpUIFigure, "TC position(s) not in node centre(s)", "Error", Icon="error")
                return
            end

            % materials
            materials = cat(3, materiallookup("Copper"), ...
                materiallookup("Inconel 718"), ...
                materiallookup("Haynes 25"), ...
                materiallookup("Inconel 718"), ...
                materiallookup("Copper"));

            % parameters
            parameters.r = app.rSpinner.Value;
            parameters.dt = app.dt.Value;
            parameters.epsilon = app.epsilonEditField.Value;
            parameters.hInitial = app.hInitialEditField.Value;
            parameters.maxiter = app.maxIterSpinner.Value;
            parameters.RTOLh = app.RTOLhEditField.Value;
            parameters.RTOLerror = app.RTOLerrorEditField.Value;

            inverse(app, filename, geometry, materials, parameters)

        end

        % Button pushed function: BrowseButton
        function BrowseButtonPushed(app, event)
            getfile(app)
        end

        % Menu selected function: SetcurrentvaluesasdefaultMenu
        function SetcurrentvaluesasdefaultMenuSelected(app, event)
            savestate(app)
        end

        % Menu selected function: ResettodefaultMenu
        function ResettodefaultMenuSelected(app, event)
            loadstate(app)
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create ihcpUIFigure and hide until all components are created
            app.ihcpUIFigure = uifigure('Visible', 'off');
            app.ihcpUIFigure.Position = [100 100 470 396];
            app.ihcpUIFigure.Name = 'ihcp';

            % Create FileMenu
            app.FileMenu = uimenu(app.ihcpUIFigure);
            app.FileMenu.Text = 'File';

            % Create OpenMenu
            app.OpenMenu = uimenu(app.FileMenu);
            app.OpenMenu.MenuSelectedFcn = createCallbackFcn(app, @OpenMenuSelected, true);
            app.OpenMenu.Accelerator = 'o';
            app.OpenMenu.Text = 'Open...';

            % Create SetcurrentvaluesasdefaultMenu
            app.SetcurrentvaluesasdefaultMenu = uimenu(app.FileMenu);
            app.SetcurrentvaluesasdefaultMenu.MenuSelectedFcn = createCallbackFcn(app, @SetcurrentvaluesasdefaultMenuSelected, true);
            app.SetcurrentvaluesasdefaultMenu.Accelerator = 'd';
            app.SetcurrentvaluesasdefaultMenu.Text = 'Set current values as default';

            % Create ResettodefaultMenu
            app.ResettodefaultMenu = uimenu(app.FileMenu);
            app.ResettodefaultMenu.MenuSelectedFcn = createCallbackFcn(app, @ResettodefaultMenuSelected, true);
            app.ResettodefaultMenu.Accelerator = 'r';
            app.ResettodefaultMenu.Text = 'Reset to default';

            % Create MattLenahan2023Menu
            app.MattLenahan2023Menu = uimenu(app.FileMenu);
            app.MattLenahan2023Menu.ForegroundColor = [0.149 0.149 0.149];
            app.MattLenahan2023Menu.Enable = 'off';
            app.MattLenahan2023Menu.Text = 'Matt Lenahan 2023';

            % Create SelectedfileEditFieldLabel
            app.SelectedfileEditFieldLabel = uilabel(app.ihcpUIFigure);
            app.SelectedfileEditFieldLabel.Position = [21 367 70 20];
            app.SelectedfileEditFieldLabel.Text = 'Selected file';

            % Create selectedfile
            app.selectedfile = uieditfield(app.ihcpUIFigure, 'text');
            app.selectedfile.Editable = 'off';
            app.selectedfile.Placeholder = '<none>';
            app.selectedfile.Position = [21 347 350 20];

            % Create RunButton
            app.RunButton = uibutton(app.ihcpUIFigure, 'push');
            app.RunButton.ButtonPushedFcn = createCallbackFcn(app, @RunButtonPushed, true);
            app.RunButton.BackgroundColor = [0.9608 0.9608 0.9608];
            app.RunButton.Enable = 'off';
            app.RunButton.Position = [21 17 430 40];
            app.RunButton.Text = 'Run';

            % Create ParametersPanel
            app.ParametersPanel = uipanel(app.ihcpUIFigure);
            app.ParametersPanel.Title = 'Parameters';
            app.ParametersPanel.FontWeight = 'bold';
            app.ParametersPanel.Position = [241 67 210 270];

            % Create rSpinnerLabel
            app.rSpinnerLabel = uilabel(app.ParametersPanel);
            app.rSpinnerLabel.HorizontalAlignment = 'right';
            app.rSpinnerLabel.Position = [51 218 40 21];
            app.rSpinnerLabel.Text = 'r';

            % Create rSpinner
            app.rSpinner = uispinner(app.ParametersPanel);
            app.rSpinner.Limits = [1 50];
            app.rSpinner.Position = [101 218 100 20];
            app.rSpinner.Value = 5;

            % Create hInitialEditFieldLabel
            app.hInitialEditFieldLabel = uilabel(app.ParametersPanel);
            app.hInitialEditFieldLabel.HorizontalAlignment = 'right';
            app.hInitialEditFieldLabel.Position = [41 139 50 20];
            app.hInitialEditFieldLabel.Text = 'hInitial';

            % Create hInitialEditField
            app.hInitialEditField = uieditfield(app.ParametersPanel, 'numeric');
            app.hInitialEditField.Limits = [1 Inf];
            app.hInitialEditField.Position = [101 137 100 22];
            app.hInitialEditField.Value = 1000;

            % Create RTOLhEditFieldLabel
            app.RTOLhEditFieldLabel = uilabel(app.ParametersPanel);
            app.RTOLhEditFieldLabel.HorizontalAlignment = 'right';
            app.RTOLhEditFieldLabel.Position = [41 176 53 22];
            app.RTOLhEditFieldLabel.Text = 'RTOLh';

            % Create RTOLhEditField
            app.RTOLhEditField = uieditfield(app.ParametersPanel, 'numeric');
            app.RTOLhEditField.Limits = [0 Inf];
            app.RTOLhEditField.Position = [101 178 100 20];
            app.RTOLhEditField.Value = 0.001;

            % Create RTOLerrorEditFieldLabel
            app.RTOLerrorEditFieldLabel = uilabel(app.ParametersPanel);
            app.RTOLerrorEditFieldLabel.HorizontalAlignment = 'right';
            app.RTOLerrorEditFieldLabel.Position = [31 159 62 20];
            app.RTOLerrorEditFieldLabel.Text = 'RTOLerror';

            % Create RTOLerrorEditField
            app.RTOLerrorEditField = uieditfield(app.ParametersPanel, 'numeric');
            app.RTOLerrorEditField.Limits = [0 Inf];
            app.RTOLerrorEditField.Position = [101 157 100 22];
            app.RTOLerrorEditField.Value = 1e-06;

            % Create maxIterSpinnerLabel
            app.maxIterSpinnerLabel = uilabel(app.ParametersPanel);
            app.maxIterSpinnerLabel.HorizontalAlignment = 'right';
            app.maxIterSpinnerLabel.Position = [51 198 45 20];
            app.maxIterSpinnerLabel.Text = 'maxIter';

            % Create maxIterSpinner
            app.maxIterSpinner = uispinner(app.ParametersPanel);
            app.maxIterSpinner.Limits = [1 50];
            app.maxIterSpinner.Position = [101 198 100 20];
            app.maxIterSpinner.Value = 20;

            % Create epsilonEditFieldLabel
            app.epsilonEditFieldLabel = uilabel(app.ParametersPanel);
            app.epsilonEditFieldLabel.HorizontalAlignment = 'right';
            app.epsilonEditFieldLabel.Position = [51 117 43 22];
            app.epsilonEditFieldLabel.Text = 'epsilon';

            % Create epsilonEditField
            app.epsilonEditField = uieditfield(app.ParametersPanel, 'numeric');
            app.epsilonEditField.Limits = [0 Inf];
            app.epsilonEditField.Position = [101 119 100 20];
            app.epsilonEditField.Value = 0.01;

            % Create dxLabel
            app.dxLabel = uilabel(app.ParametersPanel);
            app.dxLabel.HorizontalAlignment = 'right';
            app.dxLabel.Position = [61 90 30 19];
            app.dxLabel.Text = 'dx';

            % Create dx
            app.dx = uieditfield(app.ParametersPanel, 'numeric');
            app.dx.Limits = [0 Inf];
            app.dx.Position = [101 89 101 20];
            app.dx.Value = 0.000125;

            % Create dxLabel_2
            app.dxLabel_2 = uilabel(app.ParametersPanel);
            app.dxLabel_2.HorizontalAlignment = 'right';
            app.dxLabel_2.Position = [61 69 30 19];
            app.dxLabel_2.Text = 'dt';

            % Create dt
            app.dt = uieditfield(app.ParametersPanel, 'numeric');
            app.dt.Limits = [0 Inf];
            app.dt.Position = [101 70 101 20];
            app.dt.Value = 0.00025;

            % Create plotCheck
            app.plotCheck = uicheckbox(app.ParametersPanel);
            app.plotCheck.Text = 'Generate plot';
            app.plotCheck.Position = [101 39 100 20];

            % Create BrowseButton
            app.BrowseButton = uibutton(app.ihcpUIFigure, 'push');
            app.BrowseButton.ButtonPushedFcn = createCallbackFcn(app, @BrowseButtonPushed, true);
            app.BrowseButton.Position = [371 347 80 20];
            app.BrowseButton.Text = 'Browse...';

            % Create GeometryPanel
            app.GeometryPanel = uipanel(app.ihcpUIFigure);
            app.GeometryPanel.Title = 'Geometry';
            app.GeometryPanel.FontWeight = 'bold';
            app.GeometryPanel.Position = [21 197 210 140];

            % Create LengthCu1EditFieldLabel
            app.LengthCu1EditFieldLabel = uilabel(app.GeometryPanel);
            app.LengthCu1EditFieldLabel.HorizontalAlignment = 'right';
            app.LengthCu1EditFieldLabel.Position = [17 89 74 22];
            app.LengthCu1EditFieldLabel.Text = 'Length Cu 1';

            % Create lengthCu1
            app.lengthCu1 = uieditfield(app.GeometryPanel, 'numeric');
            app.lengthCu1.Position = [101 88 100 22];

            % Create LengthInco1EditFieldLabel
            app.LengthInco1EditFieldLabel = uilabel(app.GeometryPanel);
            app.LengthInco1EditFieldLabel.HorizontalAlignment = 'right';
            app.LengthInco1EditFieldLabel.Position = [6 69 85 22];
            app.LengthInco1EditFieldLabel.Text = 'Length Inco 1';

            % Create lengthInco1
            app.lengthInco1 = uieditfield(app.GeometryPanel, 'numeric');
            app.lengthInco1.Position = [101 69 100 22];

            % Create LengthH25EditFieldLabel
            app.LengthH25EditFieldLabel = uilabel(app.GeometryPanel);
            app.LengthH25EditFieldLabel.HorizontalAlignment = 'right';
            app.LengthH25EditFieldLabel.Position = [11 49 80 22];
            app.LengthH25EditFieldLabel.Text = 'Length H25 ';

            % Create lengthH25
            app.lengthH25 = uieditfield(app.GeometryPanel, 'numeric');
            app.lengthH25.Position = [101 49 100 22];

            % Create LengthInco2EditFieldLabel
            app.LengthInco2EditFieldLabel = uilabel(app.GeometryPanel);
            app.LengthInco2EditFieldLabel.HorizontalAlignment = 'right';
            app.LengthInco2EditFieldLabel.Position = [1 29 90 22];
            app.LengthInco2EditFieldLabel.Text = 'Length Inco 2 ';

            % Create lengthInco2
            app.lengthInco2 = uieditfield(app.GeometryPanel, 'numeric');
            app.lengthInco2.Position = [101 28 100 22];

            % Create LengthCu2EditFieldLabel
            app.LengthCu2EditFieldLabel = uilabel(app.GeometryPanel);
            app.LengthCu2EditFieldLabel.HorizontalAlignment = 'right';
            app.LengthCu2EditFieldLabel.Position = [11 10 80 22];
            app.LengthCu2EditFieldLabel.Text = 'Length Cu 2 ';

            % Create lengthCu2
            app.lengthCu2 = uieditfield(app.GeometryPanel, 'numeric');
            app.lengthCu2.Position = [101 9 100 22];

            % Create ThermocouplepositionsPanel
            app.ThermocouplepositionsPanel = uipanel(app.ihcpUIFigure);
            app.ThermocouplepositionsPanel.Tooltip = {'Referenced to start of first copper section input above'};
            app.ThermocouplepositionsPanel.Title = 'Thermocouple positions';
            app.ThermocouplepositionsPanel.FontWeight = 'bold';
            app.ThermocouplepositionsPanel.Position = [21 67 210 120];

            % Create T_Cu2EditFieldLabel
            app.T_Cu2EditFieldLabel = uilabel(app.ThermocouplepositionsPanel);
            app.T_Cu2EditFieldLabel.HorizontalAlignment = 'right';
            app.T_Cu2EditFieldLabel.Position = [51 67 41 22];
            app.T_Cu2EditFieldLabel.Text = 'T_Cu2';

            % Create T_Cu2pos
            app.T_Cu2pos = uieditfield(app.ThermocouplepositionsPanel, 'numeric');
            app.T_Cu2pos.Position = [101 68 100 21];

            % Create T_Inco1EditFieldLabel
            app.T_Inco1EditFieldLabel = uilabel(app.ThermocouplepositionsPanel);
            app.T_Inco1EditFieldLabel.HorizontalAlignment = 'right';
            app.T_Inco1EditFieldLabel.Position = [41 47 48 22];
            app.T_Inco1EditFieldLabel.Text = 'T_Inco1';

            % Create T_Inco1pos
            app.T_Inco1pos = uieditfield(app.ThermocouplepositionsPanel, 'numeric');
            app.T_Inco1pos.Position = [101 49 100 20];

            % Create T_Inco2EditFieldLabel
            app.T_Inco2EditFieldLabel = uilabel(app.ThermocouplepositionsPanel);
            app.T_Inco2EditFieldLabel.HorizontalAlignment = 'right';
            app.T_Inco2EditFieldLabel.Position = [41 28 48 22];
            app.T_Inco2EditFieldLabel.Text = 'T_Inco2';

            % Create T_Inco2pos
            app.T_Inco2pos = uieditfield(app.ThermocouplepositionsPanel, 'numeric');
            app.T_Inco2pos.Position = [101 28 100 22];

            % Create T_Cu3EditFieldLabel
            app.T_Cu3EditFieldLabel = uilabel(app.ThermocouplepositionsPanel);
            app.T_Cu3EditFieldLabel.HorizontalAlignment = 'right';
            app.T_Cu3EditFieldLabel.Position = [51 7 41 22];
            app.T_Cu3EditFieldLabel.Text = 'T_Cu3';

            % Create T_Cu3pos
            app.T_Cu3pos = uieditfield(app.ThermocouplepositionsPanel, 'numeric');
            app.T_Cu3pos.Position = [101 8 100 21];

            % Show the figure after all components are created
            app.ihcpUIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = ihcp_exported

            runningApp = getRunningApp(app);

            % Check for running singleton app
            if isempty(runningApp)

                % Create UIFigure and components
                createComponents(app)

                % Register the app with App Designer
                registerApp(app, app.ihcpUIFigure)

                % Execute the startup function
                runStartupFcn(app, @startupFcn)
            else

                % Focus the running singleton app
                figure(runningApp.ihcpUIFigure)

                app = runningApp;
            end

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.ihcpUIFigure)
        end
    end
end