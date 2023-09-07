classdef ihcp_exported < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        ihcpUIFigure                matlab.ui.Figure
        FileMenu                    matlab.ui.container.Menu
        OpenMenu                    matlab.ui.container.Menu
        PropertiesMenu              matlab.ui.container.Menu
        EnforceMenu                 matlab.ui.container.Menu
        SamematerialMenu            matlab.ui.container.Menu
        samematerialUC              matlab.ui.container.Menu
        samematerialCD              matlab.ui.container.Menu
        ThermocouplesPanel          matlab.ui.container.Panel
        centreTC                    matlab.ui.control.CheckBox
        HeattransfercoefficientsButtonGroup  matlab.ui.container.ButtonGroup
        AverageButton               matlab.ui.control.RadioButton
        IndividualButton            matlab.ui.control.RadioButton
        ML23Label                   matlab.ui.control.Label
        BrowseButton                matlab.ui.control.Button
        plotCheck                   matlab.ui.control.CheckBox
        ParametersPanel             matlab.ui.container.Panel
        maxIterSpinner              matlab.ui.control.Spinner
        maxIterSpinnerLabel         matlab.ui.control.Label
        RTOLerrorEditField          matlab.ui.control.NumericEditField
        RTOLerrorEditFieldLabel     matlab.ui.control.Label
        RTOLqhEditField             matlab.ui.control.NumericEditField
        RTOLqhEditFieldLabel        matlab.ui.control.Label
        hInitialEditField           matlab.ui.control.NumericEditField
        hInitialEditFieldLabel      matlab.ui.control.Label
        qInitialEditField           matlab.ui.control.NumericEditField
        qInitialEditFieldLabel      matlab.ui.control.Label
        rSpinner                    matlab.ui.control.Spinner
        rSpinnerLabel               matlab.ui.control.Label
        DownstreamPanel             matlab.ui.container.Panel
        down2length                 matlab.ui.control.NumericEditField
        LengthEditField_7Label_2    matlab.ui.control.Label
        down2mat                    matlab.ui.control.DropDown
        MaterialDropDown_7Label_2   matlab.ui.control.Label
        downconv                    matlab.ui.control.Lamp
        ConvergentLamp_3Label       matlab.ui.control.Label
        downdt                      matlab.ui.control.NumericEditField
        dxLabel_6                   matlab.ui.control.Label
        downdx                      matlab.ui.control.NumericEditField
        dxLabel_5                   matlab.ui.control.Label
        down1length                 matlab.ui.control.NumericEditField
        LengthEditField_7Label      matlab.ui.control.Label
        down1mat                    matlab.ui.control.DropDown
        MaterialDropDown_7Label     matlab.ui.control.Label
        Section2Label_5             matlab.ui.control.Label
        Section1Label_4             matlab.ui.control.Label
        CentralPanel                matlab.ui.container.Panel
        centTC                      matlab.ui.control.NumericEditField
        CentralTCEditFieldLabel     matlab.ui.control.Label
        cent3length                 matlab.ui.control.NumericEditField
        LengthEditFieldLabel_5      matlab.ui.control.Label
        cent3mat                    matlab.ui.control.DropDown
        MaterialDropDownLabel_5     matlab.ui.control.Label
        cent2length                 matlab.ui.control.NumericEditField
        LengthEditFieldLabel_4      matlab.ui.control.Label
        cent2mat                    matlab.ui.control.DropDown
        MaterialDropDownLabel_4     matlab.ui.control.Label
        cent1length                 matlab.ui.control.NumericEditField
        LengthEditFieldLabel_3      matlab.ui.control.Label
        cent1mat                    matlab.ui.control.DropDown
        MaterialDropDownLabel_3     matlab.ui.control.Label
        centconv                    matlab.ui.control.Lamp
        ConvergentLamp_2Label       matlab.ui.control.Label
        centdt                      matlab.ui.control.NumericEditField
        dxLabel_4                   matlab.ui.control.Label
        centdx                      matlab.ui.control.NumericEditField
        dxLabel_3                   matlab.ui.control.Label
        Section3Label               matlab.ui.control.Label
        Section2Label_4             matlab.ui.control.Label
        Section1Label_3             matlab.ui.control.Label
        UpstreamPanel               matlab.ui.container.Panel
        up1length                   matlab.ui.control.NumericEditField
        LengthEditFieldLabel        matlab.ui.control.Label
        up2length                   matlab.ui.control.NumericEditField
        LengthEditFieldLabel_2      matlab.ui.control.Label
        up2mat                      matlab.ui.control.DropDown
        MaterialDropDownLabel_2     matlab.ui.control.Label
        upconv                      matlab.ui.control.Lamp
        ConvergentLampLabel         matlab.ui.control.Label
        updt                        matlab.ui.control.NumericEditField
        dxLabel_2                   matlab.ui.control.Label
        updx                        matlab.ui.control.NumericEditField
        dxLabel                     matlab.ui.control.Label
        up1mat                      matlab.ui.control.DropDown
        MaterialDropDownLabel       matlab.ui.control.Label
        Section2Label               matlab.ui.control.Label
        Section1Label               matlab.ui.control.Label
        RunButton                   matlab.ui.control.Button
        selectedfile                matlab.ui.control.EditField
        SelectedfileEditFieldLabel  matlab.ui.control.Label
    end

    
    methods (Access = private)

        function checkconvergence(app)
            % lookup materials
            mat1u = materiallookup(app.up1mat.Value);
            mat2u = materiallookup(app.up2mat.Value);
            mat1c = materiallookup(app.cent1mat.Value);
            mat2c = materiallookup(app.cent2mat.Value);
            mat3c = materiallookup(app.cent3mat.Value);
            mat1d = materiallookup(app.down1mat.Value);
            mat2d = materiallookup(app.down2mat.Value);

            T = 20;
            % calculate all tau
            % upstream
            k1u = mat1u(1, 1) + T*mat1u(1, 2) + T^2*mat1u(1, 3);
            c_p1u = mat1u(2, 1) + T*mat1u(2, 2) + T^2*mat1u(2, 3);
            rho1u = mat1u(3, 1) + T*mat1u(3, 2) + T^2*mat1u(3, 3);
            tau1u = k1u*rho1u^-1*c_p1u^-1*app.updt.Value...
                *app.updx.Value^-2;
            
            k2u = mat2u(1, 1) + T*mat2u(1, 2) + T^2*mat2u(1, 3);
            c_p2u = mat2u(2, 1) + T*mat2u(2, 2) + T^2*mat2u(2, 3);
            rho2u = mat2u(3, 1) + T*mat2u(3, 2) + T^2*mat2u(3, 3);
            tau2u = k2u*rho2u^-1*c_p2u^-1*app.updt.Value...
                *app.updx.Value^-2;

            if tau1u > .5 || tau2u > .5
                app.upconv.Color = [1, 0, 0]; % red
            else
                app.upconv.Color = [0, 1, 0]; % green
            end

            % central
            k1c = mat1c(1, 1) + T*mat1c(1, 2) + T^2*mat1c(1, 3);
            c_p1c = mat1c(2, 1) + T*mat1c(2, 2) + T^2*mat1c(2, 3);
            rho1c = mat1c(3, 1) + T*mat1c(3, 2) + T^2*mat1c(3, 3);
            tau1c = k1c*rho1c^-1*c_p1c^-1*app.centdt.Value...
                *app.centdx.Value^-2;

            k2c = mat2c(1, 1) + T*mat2c(1, 2) + T^2*mat2c(1, 3);
            c_p2c = mat2c(2, 1) + T*mat2c(2, 2) + T^2*mat2c(2, 3);
            rho2c = mat2c(3, 1) + T*mat2c(3, 2) + T^2*mat2c(3, 3);
            tau2c = k2c*rho2c^-1*c_p2c^-1*app.centdt.Value...
                *app.centdx.Value^-2;

            k3c = mat3c(1, 1) + T*mat3c(1, 2) + T^2*mat3c(1, 3);
            c_p3c = mat3c(2, 1) + T*mat3c(2, 2) + T^2*mat3c(2, 3);
            rho3c = mat3c(3, 1) + T*mat3c(3, 2) + T^2*mat3c(3, 3);
            tau3c = k3c*rho3c^-1*c_p3c^-1*app.centdt.Value...
                *app.centdx.Value^-2;

            if tau1c > .5 || tau2c > .5 || tau3c > .5
                app.centconv.Color = [1, 0, 0];
            else
                app.centconv.Color = [0, 1, 0];
            end

            % downstream
            k1d = mat1d(1, 1) + T*mat1d(1, 2) + T^2*mat1d(1, 3);
            c_p1d = mat1d(2, 1) + T*mat1d(2, 2) + T^2*mat1d(2, 3);
            rho1d = mat1d(3, 1) + T*mat1d(3, 2) + T^2*mat1d(3, 3);
            tau1d = k1d*rho1d^-1*c_p1d^-1*app.downdt.Value...
                *app.downdx.Value^-2;
            
            k2d = mat2d(1, 1) + T*mat2d(1, 2) + T^2*mat2d(1, 3);
            c_p2d = mat2d(2, 1) + T*mat2d(2, 2) + T^2*mat2d(2, 3);
            rho2d = mat2d(3, 1) + T*mat2d(3, 2) + T^2*mat2d(3, 3);
            tau2d = k2d*rho2d^-1*c_p2d^-1*app.downdt.Value...
                *app.downdx.Value^-2;

            if tau1d > .5 || tau2d > .5
                app.downconv.Color = [1, 0, 0];
            else
                app.downconv.Color = [0, 1, 0];
            end

            if any([tau1u tau1c tau1d tau2u tau2c tau2d tau3c] > 0.5)
                app.RunButton.Enable = "off";
            else
                app.RunButton.Enable = "on";
            end

        end

        function getfile(app)
            % fix for lost focus
            f = figure('Position', [-100 -100 0 0]); % off-screen
            % file picker
            [file, path] = uigetfile({'*.*'; '*.csv'; '*.dat'; ...
                '*.txt'}, "Open...");
            delete(f);
            if file ~= 0
                app.selectedfile.Value = fullfile(path, file);
                app.RunButton.Enable = "on";
            end
        end

    end


    % Callbacks that handle component events
    methods (Access = private)

        % Code that executes after component creation
        function startupFcn(app)
            % enforce same material between sections
            samematerialUCSelected(app)
            samematerialCDSelected(app)
            % button icon
            app.RunButton.Icon = fullfile(matlabroot,'ui','install',...
                'installer_login','images','toolstrip','Run_24.png');
        end

        % Menu selected function: OpenMenu
        function OpenMenuSelected(app, event)
            getfile(app)
        end

        % Value changed function: updt, updx
        function updxValueChanged(app, event)
            % if upstream value is changed, downstream changes too
            app.downdx.Value = app.updx.Value;
            app.downdt.Value = app.updt.Value;
            % re-check convergence
            checkconvergence(app)
        end

        % Value changed function: centdt, centdx
        function centdxValueChanged(app, event)
            % re-check convergence
            checkconvergence(app)
        end

        % Button pushed function: RunButton
        function RunButtonPushed(app, event)
            % gather all entered data
            % filename
            filename = app.selectedfile.Value;
            
            % geometry
            % upstream
            geom.up = cell(3, 1);
            geom.up{1} = app.updx.Value;
            geom.up{2} = [app.up1length.Value app.up2length.Value];
            geom.up{3} = [geom.up{2}(1)/geom.up{1},...
                geom.up{2}(2)/geom.up{1}];
            % central
            geom.cent = cell(4, 1);
            geom.cent{1} = app.centdx.Value;
            geom.cent{2} = [app.cent1length.Value app.cent2length.Value...
                app.cent3length.Value];
            geom.cent{3} = [geom.cent{2}(1)/geom.cent{1}... 
                geom.cent{2}(2)/geom.cent{1} geom.cent{2}(3)/geom.cent{1}];
            geom.cent{4} = app.centTC.Value;
            % downstream
            geom.down = cell(3, 1);
            geom.down{1} = app.downdx.Value;
            geom.down{2} = [app.down1length.Value app.down2length.Value];
            geom.down{3} = [geom.down{2}(1)/geom.down{1}...
                geom.down{2}(2)/geom.down{1}];
            
            % materials
            mat.up = cat(3, materiallookup(app.up1mat.Value),...
                materiallookup(app.up2mat.Value));
            mat.cent = cat(3, materiallookup(app.cent1mat.Value),...
                materiallookup(app.cent2mat.Value),...
                materiallookup(app.cent3mat.Value));
            mat.down = cat(3, materiallookup(app.down1mat.Value),...
                materiallookup(app.down2mat.Value));

            % parameters
            param.r = app.rSpinner.Value;
            param.qInitial = app.qInitialEditField.Value;
            param.hInitial = app.hInitialEditField.Value;
            param.maxiter = app.maxIterSpinner.Value;
            param.RTOLqh = app.RTOLqhEditField.Value;
            param.RTOLerror = app.RTOLerrorEditField.Value;

            % run upstreamdownstream.m
            upstreamdownstream(app, filename, geom, mat, param,...
                app.updt.Value)
            % run central.m
            central(app, filename, geom, mat, param, app.centdt.Value)
        end

        % Menu selected function: samematerialUC
        function samematerialUCSelected(app, event)
            if app.samematerialUC.Checked == "on"
                % toggle
                app.samematerialUC.Checked = "off";
                % same material not enforced, box editable
                app.cent1mat.Enable = "on";
            else
                % toggle
                app.samematerialUC.Checked = "on";
                % same material enforced, not editable & change value
                app.cent1mat.Enable = "off";
                app.cent1mat.Value = app.up2mat.Value;
            end
        end

        % Menu selected function: samematerialCD
        function samematerialCDSelected(app, event)
            if app.samematerialCD.Checked == "on"
                % toggle
                app.samematerialCD.Checked = "off";
                % same material not enforced, box enabled
                app.cent3mat.Enable = "on";
            else
                % toggle
                app.samematerialCD.Checked = "on";
                % same material enforced, box disabled
                app.cent3mat.Enable = "off";
                % update value
                app.cent3mat.Value = app.down1mat.Value;
            end
        end

        % Value changed function: up2mat
        function up2matValueChanged(app, event)
            if app.samematerialUC.Checked == "on"
                app.cent1mat.Value = app.up2mat.Value;
            end
            checkconvergence(app)
        end

        % Value changed function: down1mat
        function down1matValueChanged(app, event)
            if app.samematerialCD.Checked == "on"
                app.cent3mat.Value = app.down1mat.Value;
            end
            checkconvergence(app)
        end

        % Value changed function: cent1mat, cent2mat, cent3mat, down2mat, 
        % ...and 1 other component
        function up1matValueChanged(app, event)
            checkconvergence(app)
        end

        % Button pushed function: BrowseButton
        function BrowseButtonPushed(app, event)
            getfile(app)
        end

        % Value changed function: centTC
        function centTCValueChanged(app, event)
            if app.centTC.Value > app.cent2length.Value
                app.centTC.Value = app.cent2length.Value;
            end
        end

        % Value changed function: centreTC
        function centreTCValueChanged(app, event)
            if app.centreTC.Value == 1
                app.IndividualButton.Enable = "on";
                app.centTC.Enable = "on";
            else    
                app.IndividualButton.Enable = "off";
                app.AverageButton.Value = 1;
                app.centTC.Enable = "off";
            end
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create ihcpUIFigure and hide until all components are created
            app.ihcpUIFigure = uifigure('Visible', 'off');
            app.ihcpUIFigure.Position = [100 100 659 750];
            app.ihcpUIFigure.Name = 'ihcp';
            app.ihcpUIFigure.Resize = 'off';

            % Create FileMenu
            app.FileMenu = uimenu(app.ihcpUIFigure);
            app.FileMenu.Text = 'File';

            % Create OpenMenu
            app.OpenMenu = uimenu(app.FileMenu);
            app.OpenMenu.MenuSelectedFcn = createCallbackFcn(app, @OpenMenuSelected, true);
            app.OpenMenu.Accelerator = 'o';
            app.OpenMenu.Text = 'Open...';

            % Create PropertiesMenu
            app.PropertiesMenu = uimenu(app.ihcpUIFigure);
            app.PropertiesMenu.Text = 'Properties';

            % Create EnforceMenu
            app.EnforceMenu = uimenu(app.PropertiesMenu);
            app.EnforceMenu.Text = 'Enforce...';

            % Create SamematerialMenu
            app.SamematerialMenu = uimenu(app.EnforceMenu);
            app.SamematerialMenu.Text = 'Same material...';

            % Create samematerialUC
            app.samematerialUC = uimenu(app.SamematerialMenu);
            app.samematerialUC.MenuSelectedFcn = createCallbackFcn(app, @samematerialUCSelected, true);
            app.samematerialUC.Text = 'Upstream-central';

            % Create samematerialCD
            app.samematerialCD = uimenu(app.SamematerialMenu);
            app.samematerialCD.MenuSelectedFcn = createCallbackFcn(app, @samematerialCDSelected, true);
            app.samematerialCD.Text = 'Central-downstream';

            % Create SelectedfileEditFieldLabel
            app.SelectedfileEditFieldLabel = uilabel(app.ihcpUIFigure);
            app.SelectedfileEditFieldLabel.Position = [21 711 70 20];
            app.SelectedfileEditFieldLabel.Text = 'Selected file';

            % Create selectedfile
            app.selectedfile = uieditfield(app.ihcpUIFigure, 'text');
            app.selectedfile.Editable = 'off';
            app.selectedfile.Placeholder = '<none>';
            app.selectedfile.Position = [101 711 450 20];

            % Create RunButton
            app.RunButton = uibutton(app.ihcpUIFigure, 'push');
            app.RunButton.ButtonPushedFcn = createCallbackFcn(app, @RunButtonPushed, true);
            app.RunButton.Enable = 'off';
            app.RunButton.Position = [461 41 180 90];
            app.RunButton.Text = 'Run';

            % Create UpstreamPanel
            app.UpstreamPanel = uipanel(app.ihcpUIFigure);
            app.UpstreamPanel.Title = 'Upstream';
            app.UpstreamPanel.FontWeight = 'bold';
            app.UpstreamPanel.Position = [21 521 430 170];

            % Create Section1Label
            app.Section1Label = uilabel(app.UpstreamPanel);
            app.Section1Label.Position = [81 119 60 20];
            app.Section1Label.Text = 'Section 1';

            % Create Section2Label
            app.Section2Label = uilabel(app.UpstreamPanel);
            app.Section2Label.Position = [282 117 60 22];
            app.Section2Label.Text = 'Section 2';

            % Create MaterialDropDownLabel
            app.MaterialDropDownLabel = uilabel(app.UpstreamPanel);
            app.MaterialDropDownLabel.HorizontalAlignment = 'right';
            app.MaterialDropDownLabel.Position = [11 88 58 21];
            app.MaterialDropDownLabel.Text = 'Material';

            % Create up1mat
            app.up1mat = uidropdown(app.UpstreamPanel);
            app.up1mat.Items = {'Copper', 'Inconel 718', 'Haynes 25'};
            app.up1mat.ValueChangedFcn = createCallbackFcn(app, @up1matValueChanged, true);
            app.up1mat.Position = [81 89 120 20];
            app.up1mat.Value = 'Copper';

            % Create dxLabel
            app.dxLabel = uilabel(app.UpstreamPanel);
            app.dxLabel.HorizontalAlignment = 'right';
            app.dxLabel.Position = [42 19 30 19];
            app.dxLabel.Text = 'dx';

            % Create updx
            app.updx = uieditfield(app.UpstreamPanel, 'numeric');
            app.updx.Limits = [0 Inf];
            app.updx.ValueChangedFcn = createCallbackFcn(app, @updxValueChanged, true);
            app.updx.Position = [82 19 79 20];
            app.updx.Value = 0.00025;

            % Create dxLabel_2
            app.dxLabel_2 = uilabel(app.UpstreamPanel);
            app.dxLabel_2.HorizontalAlignment = 'right';
            app.dxLabel_2.Position = [152 19 30 19];
            app.dxLabel_2.Text = 'dt';

            % Create updt
            app.updt = uieditfield(app.UpstreamPanel, 'numeric');
            app.updt.Limits = [0 Inf];
            app.updt.ValueChangedFcn = createCallbackFcn(app, @updxValueChanged, true);
            app.updt.Position = [192 19 79 20];
            app.updt.Value = 0.00025;

            % Create ConvergentLampLabel
            app.ConvergentLampLabel = uilabel(app.UpstreamPanel);
            app.ConvergentLampLabel.HorizontalAlignment = 'right';
            app.ConvergentLampLabel.Position = [301 19 67 20];
            app.ConvergentLampLabel.Text = 'Convergent';

            % Create upconv
            app.upconv = uilamp(app.UpstreamPanel);
            app.upconv.Position = [381 24 10 10];

            % Create MaterialDropDownLabel_2
            app.MaterialDropDownLabel_2 = uilabel(app.UpstreamPanel);
            app.MaterialDropDownLabel_2.HorizontalAlignment = 'right';
            app.MaterialDropDownLabel_2.Position = [210 88 58 21];
            app.MaterialDropDownLabel_2.Text = 'Material';

            % Create up2mat
            app.up2mat = uidropdown(app.UpstreamPanel);
            app.up2mat.Items = {'Copper', 'Inconel 718', 'Haynes 25'};
            app.up2mat.ValueChangedFcn = createCallbackFcn(app, @up2matValueChanged, true);
            app.up2mat.Position = [280 89 120 20];
            app.up2mat.Value = 'Inconel 718';

            % Create LengthEditFieldLabel_2
            app.LengthEditFieldLabel_2 = uilabel(app.UpstreamPanel);
            app.LengthEditFieldLabel_2.HorizontalAlignment = 'right';
            app.LengthEditFieldLabel_2.Position = [211 59 60 20];
            app.LengthEditFieldLabel_2.Text = 'Length';

            % Create up2length
            app.up2length = uieditfield(app.UpstreamPanel, 'numeric');
            app.up2length.Limits = [0 Inf];
            app.up2length.Position = [280 59 120 20];
            app.up2length.Value = 0.0065;

            % Create LengthEditFieldLabel
            app.LengthEditFieldLabel = uilabel(app.UpstreamPanel);
            app.LengthEditFieldLabel.HorizontalAlignment = 'right';
            app.LengthEditFieldLabel.Position = [12 59 60 20];
            app.LengthEditFieldLabel.Text = 'Length';

            % Create up1length
            app.up1length = uieditfield(app.UpstreamPanel, 'numeric');
            app.up1length.Limits = [0 Inf];
            app.up1length.Position = [81 59 120 20];
            app.up1length.Value = 0.0015;

            % Create CentralPanel
            app.CentralPanel = uipanel(app.ihcpUIFigure);
            app.CentralPanel.Title = 'Central';
            app.CentralPanel.FontWeight = 'bold';
            app.CentralPanel.Position = [21 321 620 190];

            % Create Section1Label_3
            app.Section1Label_3 = uilabel(app.CentralPanel);
            app.Section1Label_3.Position = [80 138 60 22];
            app.Section1Label_3.Text = 'Section 1';

            % Create Section2Label_4
            app.Section2Label_4 = uilabel(app.CentralPanel);
            app.Section2Label_4.Position = [280 138 60 22];
            app.Section2Label_4.Text = 'Section 2';

            % Create Section3Label
            app.Section3Label = uilabel(app.CentralPanel);
            app.Section3Label.Position = [482 137 60 22];
            app.Section3Label.Text = 'Section 3';

            % Create dxLabel_3
            app.dxLabel_3 = uilabel(app.CentralPanel);
            app.dxLabel_3.HorizontalAlignment = 'right';
            app.dxLabel_3.Position = [41 20 30 19];
            app.dxLabel_3.Text = 'dx';

            % Create centdx
            app.centdx = uieditfield(app.CentralPanel, 'numeric');
            app.centdx.Limits = [0 Inf];
            app.centdx.ValueChangedFcn = createCallbackFcn(app, @centdxValueChanged, true);
            app.centdx.Position = [81 20 81 20];
            app.centdx.Value = 0.0001;

            % Create dxLabel_4
            app.dxLabel_4 = uilabel(app.CentralPanel);
            app.dxLabel_4.HorizontalAlignment = 'right';
            app.dxLabel_4.Position = [151 20 30 19];
            app.dxLabel_4.Text = 'dt';

            % Create centdt
            app.centdt = uieditfield(app.CentralPanel, 'numeric');
            app.centdt.Limits = [0 Inf];
            app.centdt.ValueChangedFcn = createCallbackFcn(app, @centdxValueChanged, true);
            app.centdt.Position = [191 20 81 20];
            app.centdt.Value = 0.00125;

            % Create ConvergentLamp_2Label
            app.ConvergentLamp_2Label = uilabel(app.CentralPanel);
            app.ConvergentLamp_2Label.HorizontalAlignment = 'right';
            app.ConvergentLamp_2Label.Position = [302 19 67 20];
            app.ConvergentLamp_2Label.Text = 'Convergent';

            % Create centconv
            app.centconv = uilamp(app.CentralPanel);
            app.centconv.Position = [382 24 10 10];

            % Create MaterialDropDownLabel_3
            app.MaterialDropDownLabel_3 = uilabel(app.CentralPanel);
            app.MaterialDropDownLabel_3.HorizontalAlignment = 'right';
            app.MaterialDropDownLabel_3.Position = [10 108 58 21];
            app.MaterialDropDownLabel_3.Text = 'Material';

            % Create cent1mat
            app.cent1mat = uidropdown(app.CentralPanel);
            app.cent1mat.Items = {'Copper', 'Inconel 718', 'Haynes 25'};
            app.cent1mat.ValueChangedFcn = createCallbackFcn(app, @up1matValueChanged, true);
            app.cent1mat.Position = [80 109 120 20];
            app.cent1mat.Value = 'Inconel 718';

            % Create LengthEditFieldLabel_3
            app.LengthEditFieldLabel_3 = uilabel(app.CentralPanel);
            app.LengthEditFieldLabel_3.HorizontalAlignment = 'right';
            app.LengthEditFieldLabel_3.Position = [11 79 60 20];
            app.LengthEditFieldLabel_3.Text = 'Length';

            % Create cent1length
            app.cent1length = uieditfield(app.CentralPanel, 'numeric');
            app.cent1length.Limits = [0 Inf];
            app.cent1length.Position = [80 79 120 20];
            app.cent1length.Value = 0.0015;

            % Create MaterialDropDownLabel_4
            app.MaterialDropDownLabel_4 = uilabel(app.CentralPanel);
            app.MaterialDropDownLabel_4.HorizontalAlignment = 'right';
            app.MaterialDropDownLabel_4.Position = [210 108 58 21];
            app.MaterialDropDownLabel_4.Text = 'Material';

            % Create cent2mat
            app.cent2mat = uidropdown(app.CentralPanel);
            app.cent2mat.Items = {'Copper', 'Inconel 718', 'Haynes 25'};
            app.cent2mat.ValueChangedFcn = createCallbackFcn(app, @up1matValueChanged, true);
            app.cent2mat.Position = [280 109 120 20];
            app.cent2mat.Value = 'Haynes 25';

            % Create LengthEditFieldLabel_4
            app.LengthEditFieldLabel_4 = uilabel(app.CentralPanel);
            app.LengthEditFieldLabel_4.HorizontalAlignment = 'right';
            app.LengthEditFieldLabel_4.Position = [211 79 60 20];
            app.LengthEditFieldLabel_4.Text = 'Length';

            % Create cent2length
            app.cent2length = uieditfield(app.CentralPanel, 'numeric');
            app.cent2length.Limits = [0 Inf];
            app.cent2length.Position = [281 79 119 20];
            app.cent2length.Value = 0.005;

            % Create MaterialDropDownLabel_5
            app.MaterialDropDownLabel_5 = uilabel(app.CentralPanel);
            app.MaterialDropDownLabel_5.HorizontalAlignment = 'right';
            app.MaterialDropDownLabel_5.Position = [409 109 58 21];
            app.MaterialDropDownLabel_5.Text = 'Material';

            % Create cent3mat
            app.cent3mat = uidropdown(app.CentralPanel);
            app.cent3mat.Items = {'Copper', 'Inconel 718', 'Haynes 25'};
            app.cent3mat.ValueChangedFcn = createCallbackFcn(app, @up1matValueChanged, true);
            app.cent3mat.Position = [479 110 120 20];
            app.cent3mat.Value = 'Inconel 718';

            % Create LengthEditFieldLabel_5
            app.LengthEditFieldLabel_5 = uilabel(app.CentralPanel);
            app.LengthEditFieldLabel_5.HorizontalAlignment = 'right';
            app.LengthEditFieldLabel_5.Position = [410 80 60 20];
            app.LengthEditFieldLabel_5.Text = 'Length';

            % Create cent3length
            app.cent3length = uieditfield(app.CentralPanel, 'numeric');
            app.cent3length.Limits = [0 Inf];
            app.cent3length.Position = [479 80 120 20];
            app.cent3length.Value = 0.0015;

            % Create CentralTCEditFieldLabel
            app.CentralTCEditFieldLabel = uilabel(app.CentralPanel);
            app.CentralTCEditFieldLabel.HorizontalAlignment = 'right';
            app.CentralTCEditFieldLabel.Position = [211 49 60 20];
            app.CentralTCEditFieldLabel.Text = 'Central TC';

            % Create centTC
            app.centTC = uieditfield(app.CentralPanel, 'numeric');
            app.centTC.Limits = [0 Inf];
            app.centTC.ValueChangedFcn = createCallbackFcn(app, @centTCValueChanged, true);
            app.centTC.Tooltip = {'Distance into the central piece of material, from the left edge'};
            app.centTC.Position = [281 49 120 20];
            app.centTC.Value = 0.0015;

            % Create DownstreamPanel
            app.DownstreamPanel = uipanel(app.ihcpUIFigure);
            app.DownstreamPanel.Title = 'Downstream';
            app.DownstreamPanel.FontWeight = 'bold';
            app.DownstreamPanel.Position = [21 141 430 170];

            % Create Section1Label_4
            app.Section1Label_4 = uilabel(app.DownstreamPanel);
            app.Section1Label_4.Position = [82 119 60 20];
            app.Section1Label_4.Text = 'Section 1';

            % Create Section2Label_5
            app.Section2Label_5 = uilabel(app.DownstreamPanel);
            app.Section2Label_5.Position = [282 117 60 22];
            app.Section2Label_5.Text = 'Section 2';

            % Create MaterialDropDown_7Label
            app.MaterialDropDown_7Label = uilabel(app.DownstreamPanel);
            app.MaterialDropDown_7Label.HorizontalAlignment = 'right';
            app.MaterialDropDown_7Label.Position = [12 89 58 20];
            app.MaterialDropDown_7Label.Text = 'Material';

            % Create down1mat
            app.down1mat = uidropdown(app.DownstreamPanel);
            app.down1mat.Items = {'Copper', 'Inconel 718', 'Haynes 25'};
            app.down1mat.ValueChangedFcn = createCallbackFcn(app, @down1matValueChanged, true);
            app.down1mat.Position = [82 87 120 22];
            app.down1mat.Value = 'Inconel 718';

            % Create LengthEditField_7Label
            app.LengthEditField_7Label = uilabel(app.DownstreamPanel);
            app.LengthEditField_7Label.HorizontalAlignment = 'right';
            app.LengthEditField_7Label.Position = [12 59 60 20];
            app.LengthEditField_7Label.Text = 'Length';

            % Create down1length
            app.down1length = uieditfield(app.DownstreamPanel, 'numeric');
            app.down1length.Limits = [0 Inf];
            app.down1length.Position = [82 59 120 20];
            app.down1length.Value = 0.0065;

            % Create dxLabel_5
            app.dxLabel_5 = uilabel(app.DownstreamPanel);
            app.dxLabel_5.HorizontalAlignment = 'right';
            app.dxLabel_5.Enable = 'off';
            app.dxLabel_5.Position = [42 19 30 19];
            app.dxLabel_5.Text = 'dx';

            % Create downdx
            app.downdx = uieditfield(app.DownstreamPanel, 'numeric');
            app.downdx.Limits = [0 Inf];
            app.downdx.Enable = 'off';
            app.downdx.Position = [82 19 79 20];
            app.downdx.Value = 0.00025;

            % Create dxLabel_6
            app.dxLabel_6 = uilabel(app.DownstreamPanel);
            app.dxLabel_6.HorizontalAlignment = 'right';
            app.dxLabel_6.Enable = 'off';
            app.dxLabel_6.Position = [152 19 30 19];
            app.dxLabel_6.Text = 'dt';

            % Create downdt
            app.downdt = uieditfield(app.DownstreamPanel, 'numeric');
            app.downdt.Limits = [0 Inf];
            app.downdt.Enable = 'off';
            app.downdt.Position = [192 19 79 20];
            app.downdt.Value = 0.00025;

            % Create ConvergentLamp_3Label
            app.ConvergentLamp_3Label = uilabel(app.DownstreamPanel);
            app.ConvergentLamp_3Label.HorizontalAlignment = 'right';
            app.ConvergentLamp_3Label.Position = [301 19 67 20];
            app.ConvergentLamp_3Label.Text = 'Convergent';

            % Create downconv
            app.downconv = uilamp(app.DownstreamPanel);
            app.downconv.Position = [381 24 10 10];

            % Create MaterialDropDown_7Label_2
            app.MaterialDropDown_7Label_2 = uilabel(app.DownstreamPanel);
            app.MaterialDropDown_7Label_2.HorizontalAlignment = 'right';
            app.MaterialDropDown_7Label_2.Position = [211 89 58 20];
            app.MaterialDropDown_7Label_2.Text = 'Material';

            % Create down2mat
            app.down2mat = uidropdown(app.DownstreamPanel);
            app.down2mat.Items = {'Copper', 'Inconel 718', 'Haynes 25'};
            app.down2mat.ValueChangedFcn = createCallbackFcn(app, @up1matValueChanged, true);
            app.down2mat.Position = [281 87 120 22];
            app.down2mat.Value = 'Copper';

            % Create LengthEditField_7Label_2
            app.LengthEditField_7Label_2 = uilabel(app.DownstreamPanel);
            app.LengthEditField_7Label_2.HorizontalAlignment = 'right';
            app.LengthEditField_7Label_2.Position = [211 59 60 20];
            app.LengthEditField_7Label_2.Text = 'Length';

            % Create down2length
            app.down2length = uieditfield(app.DownstreamPanel, 'numeric');
            app.down2length.Limits = [0 Inf];
            app.down2length.Position = [281 59 120 20];
            app.down2length.Value = 0.0015;

            % Create ParametersPanel
            app.ParametersPanel = uipanel(app.ihcpUIFigure);
            app.ParametersPanel.Title = 'Parameters';
            app.ParametersPanel.FontWeight = 'bold';
            app.ParametersPanel.Position = [21 21 430 110];

            % Create rSpinnerLabel
            app.rSpinnerLabel = uilabel(app.ParametersPanel);
            app.rSpinnerLabel.HorizontalAlignment = 'right';
            app.rSpinnerLabel.Position = [11 49 40 20];
            app.rSpinnerLabel.Text = 'r';

            % Create rSpinner
            app.rSpinner = uispinner(app.ParametersPanel);
            app.rSpinner.Limits = [1 50];
            app.rSpinner.Position = [60 49 61 20];
            app.rSpinner.Value = 5;

            % Create qInitialEditFieldLabel
            app.qInitialEditFieldLabel = uilabel(app.ParametersPanel);
            app.qInitialEditFieldLabel.HorizontalAlignment = 'right';
            app.qInitialEditFieldLabel.Position = [141 49 50 20];
            app.qInitialEditFieldLabel.Text = 'qInitial';

            % Create qInitialEditField
            app.qInitialEditField = uieditfield(app.ParametersPanel, 'numeric');
            app.qInitialEditField.Limits = [1 Inf];
            app.qInitialEditField.Position = [201 47 59 22];
            app.qInitialEditField.Value = 1000;

            % Create hInitialEditFieldLabel
            app.hInitialEditFieldLabel = uilabel(app.ParametersPanel);
            app.hInitialEditFieldLabel.HorizontalAlignment = 'right';
            app.hInitialEditFieldLabel.Position = [291 49 50 20];
            app.hInitialEditFieldLabel.Text = 'hInitial';

            % Create hInitialEditField
            app.hInitialEditField = uieditfield(app.ParametersPanel, 'numeric');
            app.hInitialEditField.Limits = [1 Inf];
            app.hInitialEditField.Position = [351 49 60 20];
            app.hInitialEditField.Value = 1000;

            % Create RTOLqhEditFieldLabel
            app.RTOLqhEditFieldLabel = uilabel(app.ParametersPanel);
            app.RTOLqhEditFieldLabel.HorizontalAlignment = 'right';
            app.RTOLqhEditFieldLabel.Position = [141 17 53 22];
            app.RTOLqhEditFieldLabel.Text = 'RTOLq/h';

            % Create RTOLqhEditField
            app.RTOLqhEditField = uieditfield(app.ParametersPanel, 'numeric');
            app.RTOLqhEditField.Limits = [0 Inf];
            app.RTOLqhEditField.Position = [201 19 60 20];
            app.RTOLqhEditField.Value = 0.001;

            % Create RTOLerrorEditFieldLabel
            app.RTOLerrorEditFieldLabel = uilabel(app.ParametersPanel);
            app.RTOLerrorEditFieldLabel.HorizontalAlignment = 'right';
            app.RTOLerrorEditFieldLabel.Position = [281 19 62 20];
            app.RTOLerrorEditFieldLabel.Text = 'RTOLerror';

            % Create RTOLerrorEditField
            app.RTOLerrorEditField = uieditfield(app.ParametersPanel, 'numeric');
            app.RTOLerrorEditField.Limits = [0 Inf];
            app.RTOLerrorEditField.Position = [350 19 60 20];
            app.RTOLerrorEditField.Value = 1e-06;

            % Create maxIterSpinnerLabel
            app.maxIterSpinnerLabel = uilabel(app.ParametersPanel);
            app.maxIterSpinnerLabel.HorizontalAlignment = 'right';
            app.maxIterSpinnerLabel.Position = [6 17 45 22];
            app.maxIterSpinnerLabel.Text = 'maxIter';

            % Create maxIterSpinner
            app.maxIterSpinner = uispinner(app.ParametersPanel);
            app.maxIterSpinner.Limits = [1 50];
            app.maxIterSpinner.Position = [60 19 61 20];
            app.maxIterSpinner.Value = 20;

            % Create plotCheck
            app.plotCheck = uicheckbox(app.ihcpUIFigure);
            app.plotCheck.Text = 'Generate plot on exit';
            app.plotCheck.Position = [462 18 133 22];

            % Create BrowseButton
            app.BrowseButton = uibutton(app.ihcpUIFigure, 'push');
            app.BrowseButton.ButtonPushedFcn = createCallbackFcn(app, @BrowseButtonPushed, true);
            app.BrowseButton.Position = [561 711 80 20];
            app.BrowseButton.Text = 'Browse...';

            % Create ML23Label
            app.ML23Label = uilabel(app.ihcpUIFigure);
            app.ML23Label.HorizontalAlignment = 'right';
            app.ML23Label.FontColor = [0.651 0.651 0.651];
            app.ML23Label.Position = [600 18 41 22];
            app.ML23Label.Text = 'ML 23';

            % Create HeattransfercoefficientsButtonGroup
            app.HeattransfercoefficientsButtonGroup = uibuttongroup(app.ihcpUIFigure);
            app.HeattransfercoefficientsButtonGroup.Title = 'Heat transfer coefficients';
            app.HeattransfercoefficientsButtonGroup.FontWeight = 'bold';
            app.HeattransfercoefficientsButtonGroup.Position = [462 141 179 169];

            % Create IndividualButton
            app.IndividualButton = uiradiobutton(app.HeattransfercoefficientsButtonGroup);
            app.IndividualButton.Text = 'Individual';
            app.IndividualButton.Position = [11 123 73 22];
            app.IndividualButton.Value = true;

            % Create AverageButton
            app.AverageButton = uiradiobutton(app.HeattransfercoefficientsButtonGroup);
            app.AverageButton.Text = 'Average';
            app.AverageButton.Position = [11 101 66 22];

            % Create ThermocouplesPanel
            app.ThermocouplesPanel = uipanel(app.ihcpUIFigure);
            app.ThermocouplesPanel.Title = 'Thermocouples';
            app.ThermocouplesPanel.FontWeight = 'bold';
            app.ThermocouplesPanel.Position = [461 521 180 170];

            % Create centreTC
            app.centreTC = uicheckbox(app.ThermocouplesPanel);
            app.centreTC.ValueChangedFcn = createCallbackFcn(app, @centreTCValueChanged, true);
            app.centreTC.Text = 'Use TC in central piece';
            app.centreTC.Position = [11 117 147 22];
            app.centreTC.Value = true;

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