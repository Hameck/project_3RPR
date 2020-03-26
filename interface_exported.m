classdef interface_exported < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        RPRAnalysistoolByAMazouniHMartinezUIFigure  matlab.ui.Figure
        TabGroup                       matlab.ui.container.TabGroup
        ParametersTab                  matlab.ui.container.Tab
        AregularcylindricalworkspacewithadiameterofLabel  matlab.ui.control.Label
        ProblemParametersandconstraintsLabel  matlab.ui.control.Label
        workspaceDiameter              matlab.ui.control.EditField
        ConstraintsLabel               matlab.ui.control.Label
        WorkspaceLabel                 matlab.ui.control.Label
        mLabel_17                      matlab.ui.control.Label
        condEnabled                    matlab.ui.control.CheckBox
        lowerLimitEnabled              matlab.ui.control.CheckBox
        highLimitEnabled               matlab.ui.control.CheckBox
        Label                          matlab.ui.control.Label
        conditionnement                matlab.ui.control.EditField
        Label_2                        matlab.ui.control.Label
        lowerLimitValue                matlab.ui.control.EditField
        Label_3                        matlab.ui.control.Label
        higherLimitValue               matlab.ui.control.EditField
        Label_5                        matlab.ui.control.Label
        WorkspaceSamplingLabel         matlab.ui.control.Label
        MobileplatformerotationrangeLabel  matlab.ui.control.Label
        minPhi                         matlab.ui.control.Spinner
        maxPhi                         matlab.ui.control.Spinner
        toLabel                        matlab.ui.control.Label
        Label_6                        matlab.ui.control.Label
        SampletheworkspaceintoLabel    matlab.ui.control.Label
        cerclesandLabel                matlab.ui.control.Label
        equallyspacedsudivisionsLabel  matlab.ui.control.Label
        ncercles                       matlab.ui.control.Spinner
        ndivs                          matlab.ui.control.Spinner
        samplingPreview                matlab.ui.control.UIAxes
        OptimisationTab                matlab.ui.container.Tab
        OptimisationparametersLabel    matlab.ui.control.Label
        ThesystemconfigurationXcontainsDdLabel  matlab.ui.control.Label
        LowerboundsPanel               matlab.ui.container.Panel
        DLabel                         matlab.ui.control.Label
        mLabel                         matlab.ui.control.Label
        lb_D                           matlab.ui.control.NumericEditField
        mLabel_2                       matlab.ui.control.Label
        lb_d                           matlab.ui.control.NumericEditField
        dLabel                         matlab.ui.control.Label
        UpperboundsPanel               matlab.ui.container.Panel
        DLabel_2                       matlab.ui.control.Label
        mLabel_5                       matlab.ui.control.Label
        ub_D                           matlab.ui.control.NumericEditField
        mLabel_6                       matlab.ui.control.Label
        ub_d                           matlab.ui.control.NumericEditField
        dLabel_2                       matlab.ui.control.Label
        InitialStatePanel              matlab.ui.container.Panel
        DLabel_3                       matlab.ui.control.Label
        mLabel_9                       matlab.ui.control.Label
        D0                             matlab.ui.control.NumericEditField
        mLabel_10                      matlab.ui.control.Label
        d0                             matlab.ui.control.NumericEditField
        dLabel_3                       matlab.ui.control.Label
        MaxiterationsLabel             matlab.ui.control.Label
        MaxfunctionevaluationLabel     matlab.ui.control.Label
        MaxIterations                  matlab.ui.control.NumericEditField
        MaxEvals                       matlab.ui.control.NumericEditField
        LogIterations                  matlab.ui.control.CheckBox
        Label_7                        matlab.ui.control.Label
        RunOptimisationButton          matlab.ui.control.Button
        ExitStatus                     matlab.ui.control.Lamp
        ExitFlag                       matlab.ui.control.Label
        ResultsTab                     matlab.ui.container.Tab
        TabGroup2                      matlab.ui.container.TabGroup
        OptimalConfigTab               matlab.ui.container.Tab
        OptiConfigAx                   matlab.ui.control.UIAxes
        InitialConfigTab               matlab.ui.container.Tab
        InitConfigAx                   matlab.ui.control.UIAxes
        XmSliderLabel                  matlab.ui.control.Label
        XmSlider                       matlab.ui.control.Slider
        YmSliderLabel                  matlab.ui.control.Label
        YmSlider                       matlab.ui.control.Slider
        PhiKnobLabel                   matlab.ui.control.Label
        PhiKnob                        matlab.ui.control.Knob
        OptimaleConfigurationPanel     matlab.ui.container.Panel
        DLabel_4                       matlab.ui.control.Label
        mLabel_13                      matlab.ui.control.Label
        X_D                            matlab.ui.control.NumericEditField
        mLabel_14                      matlab.ui.control.Label
        X_d                            matlab.ui.control.NumericEditField
        dLabel_4                       matlab.ui.control.Label
        CurrentPositionLabel           matlab.ui.control.Label
        XLabel                         matlab.ui.control.Label
        currentX                       matlab.ui.control.NumericEditField
        mLabel_15                      matlab.ui.control.Label
        mLabel_16                      matlab.ui.control.Label
        currentY                       matlab.ui.control.NumericEditField
        YLabel                         matlab.ui.control.Label
        Label_8                        matlab.ui.control.Label
        currentPhi                     matlab.ui.control.NumericEditField
        OptimisationresultsLabel       matlab.ui.control.Label
        PhiLabel                       matlab.ui.control.Label
        GoToButton                     matlab.ui.control.Button
        ValidationTab                  matlab.ui.container.Tab
        ValidationAxis                 matlab.ui.control.UIAxes
        AnalyseSolutionButton          matlab.ui.control.Button
        ThistakesafewminutesPleasewaitLabel  matlab.ui.control.Label
    end

    
    properties (Access = private)
        runningOptimisation = false; % Description
        pointsList = []; % Description
        initConfig = [10 1]; % Description
        optimalConfig = [10 1]; % Description
    end
    
    methods (Access = private)
        function [points_x, points_y] = generateCerclePoints(~, diameter, ndivs)
            points_x = zeros([1 ndivs]);
            points_y = zeros([1 ndivs]);
            angles = linspace(0,2*pi,ndivs);
            for i=1:ndivs
                points_x(i) = (diameter/2) * cos(angles(i));
                points_y(i) = (diameter/2) * sin(angles(i));
            end
        end
        function [points_x, points_y] = generatePoints(app, diameter, n_cercles, ndivs)
            points_x = [];
            points_y = [];
            r = linspace(0,diameter/2, n_cercles);
            for i=1:n_cercles
              [xs, ys] = generateCerclePoints(app, r(i)*2, ndivs);
              points_x = [points_x, xs];
              points_y = [points_y, ys];
            end
            % Store pointsList
            app.pointsList = [];
            for i=1:length(points_x)
                app.pointsList = [app.pointsList, [points_x(i);points_y(i)]];
            end
        end
        function updateSamplingPreview(app)
            n_cercles = app.ncercles.Value + 1;
            n_divs = app.ndivs.Value + 1;
            diameter = str2double(app.workspaceDiameter.Value);
            [points_x, points_y] = generatePoints(app, diameter, n_cercles, n_divs);
            cla(app.samplingPreview);
            hold(app.samplingPreview,'on');
            rectangle(app.samplingPreview, 'Position',[-diameter/2,-diameter/2,diameter,diameter],'Curvature',[1 1]);
            plot(app.samplingPreview, points_x, points_y,'r*');
            axis(app.samplingPreview, 'equal');
            hold(app.samplingPreview,'off');
        end
        
        function toggleOptimisation(app, status)
            app.runningOptimisation = status;
            if status
                % Start
                app.Running.Enable = 'on';
                app.Running.Color = '0.00,1.00,0.00';
            else
                % Stop optimization
                app.Running.Color = '0.90,0.90,0.90';
                app.Running.Enable = 'off';
            end
        end
        
        function runOptimisation(app)
            % global parameters
            global conditionnement pointsList minPhi maxPhi maxJointLimit minJointLimit maxJointLimitEnabled minJointLimitEnabled condEnabled;
            
            conditionnement = str2double(app.conditionnement.Value);
            pointsList = app.pointsList;
            minPhi = app.minPhi.Value;
            maxPhi = app.maxPhi.Value;
            maxJointLimit = str2double(app.higherLimitValue.Value);
            minJointLimit = str2double(app.lowerLimitValue.Value);
            
            maxJointLimitEnabled = app.highLimitEnabled.Value;
            minJointLimitEnabled = app.lowerLimitEnabled.Value;
            condEnabled = app.condEnabled.Value;
            
            % initial conditions
            x0 = app.initConfig; % D d
            lb =[app.lb_D.Value; app.lb_d.Value];
            ub =[app.ub_D.Value; app.ub_d.Value];
            
            display = 'off';
            if app.LogIterations.Value
                display = 'iter';
            end
            app.ExitFlag.Text = "Running...";
            app.ExitStatus.Color = "0.94,0.94,0.94";
            app.ExitStatus.Enable = "off";
            options = optimset('Display',display, 'MaxIter', app.MaxIterations.Value,'MaxFunEvals', app.MaxEvals.Value, 'TolFun',1e-18, 'TolX',1e-18);
            [X,~,EXITFLAG] = fmincon('objectiveFun', x0, [], [], [], [], lb, ub , 'nonlcong', options);
            setExitFlag(app, EXITFLAG);
            app.optimalConfig = X;
            app.X_D.Value = X(1);
            app.X_d.Value = X(2);
            drawOptiConfig(app, app.optimalConfig);
        end
        
        function setExitFlag(app, flag)
            app.ExitStatus.Enable = "on";
            switch flag
                case 1
                    app.ExitFlag.Text = "First order optimality conditions satisfied.";
                    app.ExitStatus.Color = "0.00,1.00,0.00";
                case 0
                    app.ExitFlag.Text = "Too many function evaluations or iterations.";
                    app.ExitStatus.Color = "0.85,0.33,0.10";
                case -1
                    app.ExitFlag.Text = "Stopped by output/plot function.";
                    app.ExitStatus.Color = "0.93,0.69,0.13";
                case -2
                    app.ExitFlag.Text = "No feasible point found.";
                    app.ExitStatus.Color = "0.85,0.33,0.10";
                case 2
                    app.ExitFlag.Text = "Local minimum possible. Constraints satisfied. Change in X too small.";
                    app.ExitStatus.Color = "0.00,1.00,0.00";
                case 3
                    app.ExitFlag.Text = "Change in objective function too small.";
                    app.ExitStatus.Color = "0.85,0.33,0.10";
                case 4
                    app.ExitFlag.Text = "Computed search direction too small.";
                    app.ExitStatus.Color = "0.85,0.33,0.10";
                case 5
                    app.ExitFlag.Text = "Predicted change in objective function too small.";
                    app.ExitStatus.Color = "0.85,0.33,0.10";
                case -3
                    app.ExitFlag.Text = "Problem seems unbounded.";
                    app.ExitStatus.Color = "0.85,0.33,0.10";
                otherwise
                    app.ExitFlag.Text = "Unknown flag";
                    app.ExitStatus.Color = "0.85,0.33,0.10";
            end
        end
        
        function drawInitConfig(app, x)
            % Draws the robot using the given configuration x
            D = x(1);
            d = x(2);
    
            A_1 = [0 0]';
            A_2 = [D 0]';
            A_3 = [D*1/2 D*sqrt(3)/2]';
            
            x = app.XmSlider.Value;
            y = app.YmSlider.Value;
            phi = deg2rad(app.PhiKnob.Value);
            
            % inverse kinematics
            B_1 = [x+d*cos(phi-5*pi/6) y+d*sin(phi-5*pi/6)]';
            B_2 = [x+d*cos(phi-pi/6) y+d*sin(phi-pi/6)]';
            B_3 = [x+d*cos(phi+pi/2) y+d*sin(phi+pi/2)]';
    
            P = [x y]';
            
            % Plot init
            cla(app.InitConfigAx);
            hold(app.InitConfigAx, "on");
            plot(app.InitConfigAx, [A_1(1) B_1(1)],[A_1(2) B_1(2)], '-o', 'linewidth', 2);
            plot(app.InitConfigAx, [A_2(1) B_2(1)],[A_2(2) B_2(2)], '-o', 'linewidth', 2);
            plot(app.InitConfigAx, [A_3(1) B_3(1)],[A_3(2) B_3(2)], '-o', 'linewidth', 2);
            plot(app.InitConfigAx, [B_1(1) B_2(1) B_3(1) B_1(1)],[B_1(2) B_2(2) B_3(2) B_1(2)], '-o', 'linewidth', 2);
            plot(app.InitConfigAx, [P(1)],[P(2)], '-+', 'linewidth', 1);
            grid(app.InitConfigAx, "on");
            xlim(app.OptiConfigAx, [0 D]);
            ylim(app.OptiConfigAx, [0 D]);
            axis(app.OptiConfigAx, "equal");
            hold(app.InitConfigAx, "off");
        end
        
        function drawOptiConfig(app, x)
            % Draws the robot using the given configuration x
            D = x(1);
            d = x(2);
    
            A_1 = [0 0]';
            A_2 = [D 0]';
            A_3 = [D*1/2 D*sqrt(3)/2]';
            
            x = app.XmSlider.Value;
            y = app.YmSlider.Value;
            phi = deg2rad(app.PhiKnob.Value);
            
            % inverse kinematics
            B_1 = [x+d*cos(phi-5*pi/6) y+d*sin(phi-5*pi/6)]';
            B_2 = [x+d*cos(phi-pi/6) y+d*sin(phi-pi/6)]';
            B_3 = [x+d*cos(phi+pi/2) y+d*sin(phi+pi/2)]';
    
            P = [x y]';
            
            % Plot init
            cla(app.OptiConfigAx);
            hold(app.OptiConfigAx, "on");
            plot(app.OptiConfigAx, [A_1(1) B_1(1)],[A_1(2) B_1(2)], '-o', 'linewidth', 2);
            plot(app.OptiConfigAx, [A_2(1) B_2(1)],[A_2(2) B_2(2)], '-o', 'linewidth', 2);
            plot(app.OptiConfigAx, [A_3(1) B_3(1)],[A_3(2) B_3(2)], '-o', 'linewidth', 2);
            plot(app.OptiConfigAx, [B_1(1) B_2(1) B_3(1) B_1(1)],[B_1(2) B_2(2) B_3(2) B_1(2)], '-o', 'linewidth', 2);
            plot(app.OptiConfigAx, [P(1)],[P(2)], '-+', 'linewidth', 1);
            grid(app.OptiConfigAx, "on");
            xlim(app.OptiConfigAx, [0 D]);
            ylim(app.OptiConfigAx, [0 D]);
            axis(app.OptiConfigAx, "equal");
            hold(app.OptiConfigAx, "off");
        end
    end
    

    % Callbacks that handle component events
    methods (Access = private)

        % Code that executes after component creation
        function startupFcn(app)
            % Call update sampling function
            updateSamplingPreview(app);
            drawInitConfig(app, app.initConfig);
            drawOptiConfig(app, app.optimalConfig);
            movegui(app.RPRAnalysistoolByAMazouniHMartinezUIFigure,'center');
        end

        % Value changed function: ncercles
        function ncerclesValueChanged(app, event)
            % Call update sampling function
            updateSamplingPreview(app);
        end

        % Value changed function: ndivs
        function ndivsValueChanged(app, event)
            % Call update sampling function
            updateSamplingPreview(app);           
        end

        % Value changed function: workspaceDiameter
        function workspaceDiameterValueChanged(app, event)
            % Call update sampling function
            updateSamplingPreview(app);
        end

        % Button pushed function: RunOptimisationButton
        function RunOptimisationButtonPushed(app, event)
            % if app.runningOptimisation
            %     Stop optimization
            %     toggleOptimisation(app, false)
            % else
            %     Start
            %     toggleOptimisation(app, true)
            %     
            % end
            runOptimisation(app)
        end

        % Value changed function: D0
        function D0ValueChanged(app, event)
            app.initConfig = [app.D0.Value; app.d0.Value];
        end

        % Value changed function: d0
        function d0ValueChanged(app, event)
            app.initConfig = [app.D0.Value; app.d0.Value];
        end

        % Value changed function: XmSlider
        function XmSliderValueChanged(app, event)
            drawInitConfig(app, app.initConfig);
            drawOptiConfig(app, app.optimalConfig);
            app.currentX.Value = app.XmSlider.Value;
        end

        % Value changed function: YmSlider
        function YmSliderValueChanged(app, event)
            drawInitConfig(app, app.initConfig);
            drawOptiConfig(app, app.optimalConfig);
            app.currentY.Value = app.YmSlider.Value;
        end

        % Value changed function: PhiKnob
        function PhiKnobValueChanged(app, event)
            drawInitConfig(app, app.initConfig);
            drawOptiConfig(app, app.optimalConfig);
            app.currentPhi.Value = app.PhiKnob.Value;
        end

        % Button pushed function: GoToButton
        function GoToButtonPushed(app, event)
            app.XmSlider.Value = app.currentX.Value;
            app.YmSlider.Value = app.currentY.Value;
            app.PhiKnob.Value = app.currentPhi.Value;
            drawInitConfig(app, app.initConfig);
            drawOptiConfig(app, app.optimalConfig);
        end

        % Button pushed function: AnalyseSolutionButton
        function AnalyseSolutionButtonPushed(app, event)
            % Draws the robot using the given configuration x
            conditionnement = str2double(app.conditionnement.Value);
            minPhi = app.minPhi.Value;
            maxPhi = app.maxPhi.Value;
            x = app.optimalConfig;
            D = x(1);
            d = x(2);
            
            A_1 = [0 0]';
            A_2 = [D 0]';
            A_3 = [D*1/2 D*sqrt(3)/2]';

            pointsList = [];
            N = 15;
            taille = 1.5;
            
            for i=-(taille/2):(taille/N):(taille/2)
                for j=-(taille/2):(taille/N):(taille/2)
                    pointsList = [pointsList, [i+D/2;j+D*sqrt(3)/6]];
                end
            end

            % end-effector position
            x = D/2;
            y = D*sqrt(3)/6;
            phi = deg2rad((maxPhi-minPhi)/2);
            
            % inverse kinematics
            B_1 = [x+d*cos(phi-5*pi/6) y+d*sin(phi-5*pi/6)]';
            B_2 = [x+d*cos(phi-pi/6) y+d*sin(phi-pi/6)]';
            B_3 = [x+d*cos(phi+pi/2) y+d*sin(phi+pi/2)]';
            
            cla(app.ValidationAxis);
            plot(app.ValidationAxis, [A_1(1) A_2(1) A_3(1) A_1(1)],[A_1(2) A_2(2) A_3(2) A_1(2)], "--ok", "linewidth", 3); hold(app.ValidationAxis,'on');
            plot(app.ValidationAxis, [A_1(1) B_1(1)],[A_1(2) B_1(2)], '-o', 'linewidth', 3); hold(app.ValidationAxis,'on');
            plot(app.ValidationAxis, [A_2(1) B_2(1)],[A_2(2) B_2(2)], '-o', 'linewidth', 3); hold(app.ValidationAxis,'on');
            plot(app.ValidationAxis, [A_3(1) B_3(1)],[A_3(2) B_3(2)], '-o', 'linewidth', 3); hold(app.ValidationAxis,'on');
            plot(app.ValidationAxis, [B_1(1) B_2(1) B_3(1) B_1(1)],[B_1(2) B_2(2) B_3(2) B_1(2)], '-o', 'linewidth', 3);hold(app.ValidationAxis,'on');
            
            for i=1:length(pointsList)
              for phi=deg2rad(minPhi):pi/36:deg2rad(maxPhi)
                
                % end-effector position
                x = pointsList(1,i);
                y = pointsList(2,i);

                % inverse kinematics
                B_1 = [x+d*cos(phi-5*pi/6) y+d*sin(phi-5*pi/6)]';
                B_2 = [x+d*cos(phi-pi/6) y+d*sin(phi-pi/6)]';
                B_3 = [x+d*cos(phi+pi/2) y+d*sin(phi+pi/2)]';

                P = [x y]';

                rho_1 = norm(B_1 - A_1);
                rho_2 = norm(B_2 - A_2);
                rho_3 = norm(B_3 - A_3);


                % jacobian
                E = [0, -1; 1, 0];
                v_1 = (B_1-A_1)/norm(B_1-A_1);
                v_2 = (B_2-A_2)/norm(B_2-A_2);
                v_3 = (B_3-A_3)/norm(B_3-A_3);

                % paralel Jacobian
                A = [(v_1)',-(v_1)'*E*(P-B_1); ...
                     (v_2)',-(v_2)'*E*(P-B_2); ...
                     (v_3)',-(v_3)'*E*(P-B_3)];

                M = isnan(A);

                if(det(A) == 0 || any(M(:)))
                  J = [];
                  if(det(A) == 0)
                    disp('Paralel singularity');
                  else
                    disp('any M');  
                  end
                else
                  % serial Jacobian 
                  B = eye(3);

                  % condition
                  J = inv(A)*B;
                end

                if((conditionnement*cond(J)-1)>0)
                    plot(app.ValidationAxis, P(1),P(2), '-+r', 'linewidth', 2);hold(app.ValidationAxis,'on');
                else
                    plot(app.ValidationAxis, P(1),P(2), '-+k', 'linewidth', 2);hold(app.ValidationAxis,'on');
                end
              end
            end
            rectangle(app.ValidationAxis, 'Position',[D/2-0.5, D*sqrt(3)/6-0.5, 1, 1],'Curvature',1,'LineWidth',3);
            title(app.ValidationAxis, strcat(['D = ',num2str(D),' and d = ',num2str(d)]));
            axis(app.ValidationAxis, 'equal');
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create RPRAnalysistoolByAMazouniHMartinezUIFigure and hide until all components are created
            app.RPRAnalysistoolByAMazouniHMartinezUIFigure = uifigure('Visible', 'off');
            app.RPRAnalysistoolByAMazouniHMartinezUIFigure.Position = [100 100 928 551];
            app.RPRAnalysistoolByAMazouniHMartinezUIFigure.Name = '3-RPR Analysis tool By A. Mazouni, H. Martinez';
            app.RPRAnalysistoolByAMazouniHMartinezUIFigure.Resize = 'off';

            % Create TabGroup
            app.TabGroup = uitabgroup(app.RPRAnalysistoolByAMazouniHMartinezUIFigure);
            app.TabGroup.TabLocation = 'bottom';
            app.TabGroup.Position = [1 -1 927 553];

            % Create ParametersTab
            app.ParametersTab = uitab(app.TabGroup);
            app.ParametersTab.Title = 'Parameters';

            % Create AregularcylindricalworkspacewithadiameterofLabel
            app.AregularcylindricalworkspacewithadiameterofLabel = uilabel(app.ParametersTab);
            app.AregularcylindricalworkspacewithadiameterofLabel.Position = [52 317 275 22];
            app.AregularcylindricalworkspacewithadiameterofLabel.Text = 'A regular cylindrical workspace with a diameter of ';

            % Create ProblemParametersandconstraintsLabel
            app.ProblemParametersandconstraintsLabel = uilabel(app.ParametersTab);
            app.ProblemParametersandconstraintsLabel.FontSize = 20;
            app.ProblemParametersandconstraintsLabel.FontWeight = 'bold';
            app.ProblemParametersandconstraintsLabel.Position = [53 481 359 25];
            app.ProblemParametersandconstraintsLabel.Text = 'Problem Parameters and constraints';

            % Create workspaceDiameter
            app.workspaceDiameter = uieditfield(app.ParametersTab, 'text');
            app.workspaceDiameter.ValueChangedFcn = createCallbackFcn(app, @workspaceDiameterValueChanged, true);
            app.workspaceDiameter.HorizontalAlignment = 'center';
            app.workspaceDiameter.Position = [326 315 68 22];
            app.workspaceDiameter.Value = '1';

            % Create ConstraintsLabel
            app.ConstraintsLabel = uilabel(app.ParametersTab);
            app.ConstraintsLabel.FontSize = 16;
            app.ConstraintsLabel.FontWeight = 'bold';
            app.ConstraintsLabel.Position = [53 280 100 22];
            app.ConstraintsLabel.Text = 'Constraints:';

            % Create WorkspaceLabel
            app.WorkspaceLabel = uilabel(app.ParametersTab);
            app.WorkspaceLabel.FontSize = 16;
            app.WorkspaceLabel.FontWeight = 'bold';
            app.WorkspaceLabel.Position = [53 350 96 22];
            app.WorkspaceLabel.Text = 'Workspace:';

            % Create mLabel_17
            app.mLabel_17 = uilabel(app.ParametersTab);
            app.mLabel_17.Position = [402 315 25 22];
            app.mLabel_17.Text = 'm';

            % Create condEnabled
            app.condEnabled = uicheckbox(app.ParametersTab);
            app.condEnabled.Text = '';
            app.condEnabled.Position = [53 248 25 22];
            app.condEnabled.Value = true;

            % Create lowerLimitEnabled
            app.lowerLimitEnabled = uicheckbox(app.ParametersTab);
            app.lowerLimitEnabled.Text = '';
            app.lowerLimitEnabled.Position = [53 206 25 22];
            app.lowerLimitEnabled.Value = true;

            % Create highLimitEnabled
            app.highLimitEnabled = uicheckbox(app.ParametersTab);
            app.highLimitEnabled.Text = '';
            app.highLimitEnabled.Position = [53 165 25 22];
            app.highLimitEnabled.Value = true;

            % Create Label
            app.Label = uilabel(app.ParametersTab);
            app.Label.Position = [77 246 401 22];
            app.Label.Text = 'The inverse conditioning of the kinematic Jacobians must be greater than';

            % Create conditionnement
            app.conditionnement = uieditfield(app.ParametersTab, 'text');
            app.conditionnement.HorizontalAlignment = 'center';
            app.conditionnement.Position = [472 244 68 22];
            app.conditionnement.Value = '0.1';

            % Create Label_2
            app.Label_2 = uilabel(app.ParametersTab);
            app.Label_2.Position = [77 204 401 22];
            app.Label_2.Text = 'The minimum prismatic joint movement is limited to';

            % Create lowerLimitValue
            app.lowerLimitValue = uieditfield(app.ParametersTab, 'text');
            app.lowerLimitValue.HorizontalAlignment = 'center';
            app.lowerLimitValue.Position = [361 202 68 22];
            app.lowerLimitValue.Value = '0.5';

            % Create Label_3
            app.Label_3 = uilabel(app.ParametersTab);
            app.Label_3.Position = [77 163 283 22];
            app.Label_3.Text = 'The maximum prismatic joint movement is limited to';

            % Create higherLimitValue
            app.higherLimitValue = uieditfield(app.ParametersTab, 'text');
            app.higherLimitValue.HorizontalAlignment = 'center';
            app.higherLimitValue.Position = [364 162 68 22];
            app.higherLimitValue.Value = '2.5';

            % Create Label_5
            app.Label_5 = uilabel(app.ParametersTab);
            app.Label_5.BackgroundColor = [0.8 0.8 0.8];
            app.Label_5.HorizontalAlignment = 'center';
            app.Label_5.FontSize = 17;
            app.Label_5.Position = [53 389 818 71];
            app.Label_5.Text = {'This application optimizes the design of a 3-RPR parallel robot, according to certain constraints.'; 'Unchecked constraints are not taken into account in the optimization.'};

            % Create WorkspaceSamplingLabel
            app.WorkspaceSamplingLabel = uilabel(app.ParametersTab);
            app.WorkspaceSamplingLabel.FontSize = 16;
            app.WorkspaceSamplingLabel.FontWeight = 'bold';
            app.WorkspaceSamplingLabel.Position = [53 121 173 22];
            app.WorkspaceSamplingLabel.Text = 'Workspace Sampling:';

            % Create MobileplatformerotationrangeLabel
            app.MobileplatformerotationrangeLabel = uilabel(app.ParametersTab);
            app.MobileplatformerotationrangeLabel.Position = [54 72 172 22];
            app.MobileplatformerotationrangeLabel.Text = 'Mobile platforme rotation range';

            % Create minPhi
            app.minPhi = uispinner(app.ParametersTab);
            app.minPhi.Position = [229 71 62 22];
            app.minPhi.Value = -65;

            % Create maxPhi
            app.maxPhi = uispinner(app.ParametersTab);
            app.maxPhi.Position = [326 70 62 22];
            app.maxPhi.Value = -5;

            % Create toLabel
            app.toLabel = uilabel(app.ParametersTab);
            app.toLabel.Position = [299 71 29 22];
            app.toLabel.Text = '° to';

            % Create Label_6
            app.Label_6 = uilabel(app.ParametersTab);
            app.Label_6.Position = [394 71 29 22];
            app.Label_6.Text = '°';

            % Create SampletheworkspaceintoLabel
            app.SampletheworkspaceintoLabel = uilabel(app.ParametersTab);
            app.SampletheworkspaceintoLabel.Position = [53 35 153 22];
            app.SampletheworkspaceintoLabel.Text = 'Sample the workspace into ';

            % Create cerclesandLabel
            app.cerclesandLabel = uilabel(app.ParametersTab);
            app.cerclesandLabel.Position = [279 35 70 22];
            app.cerclesandLabel.Text = 'cercles and ';

            % Create equallyspacedsudivisionsLabel
            app.equallyspacedsudivisionsLabel = uilabel(app.ParametersTab);
            app.equallyspacedsudivisionsLabel.Position = [421 35 151 22];
            app.equallyspacedsudivisionsLabel.Text = 'equally spaced sudivisions.';

            % Create ncercles
            app.ncercles = uispinner(app.ParametersTab);
            app.ncercles.Limits = [1 Inf];
            app.ncercles.ValueChangedFcn = createCallbackFcn(app, @ncerclesValueChanged, true);
            app.ncercles.Position = [205 35 62 22];
            app.ncercles.Value = 2;

            % Create ndivs
            app.ndivs = uispinner(app.ParametersTab);
            app.ndivs.Limits = [4 Inf];
            app.ndivs.ValueChangedFcn = createCallbackFcn(app, @ndivsValueChanged, true);
            app.ndivs.Position = [348 35 62 22];
            app.ndivs.Value = 5;

            % Create samplingPreview
            app.samplingPreview = uiaxes(app.ParametersTab);
            title(app.samplingPreview, 'Sampling Preview')
            xlabel(app.samplingPreview, 'X (m)')
            ylabel(app.samplingPreview, 'Y(m)')
            app.samplingPreview.Position = [600 55 317 309];

            % Create OptimisationTab
            app.OptimisationTab = uitab(app.TabGroup);
            app.OptimisationTab.Title = 'Optimisation';

            % Create OptimisationparametersLabel
            app.OptimisationparametersLabel = uilabel(app.OptimisationTab);
            app.OptimisationparametersLabel.FontSize = 24;
            app.OptimisationparametersLabel.FontWeight = 'bold';
            app.OptimisationparametersLabel.Position = [53 477 291 29];
            app.OptimisationparametersLabel.Text = 'Optimisation parameters';

            % Create ThesystemconfigurationXcontainsDdLabel
            app.ThesystemconfigurationXcontainsDdLabel = uilabel(app.OptimisationTab);
            app.ThesystemconfigurationXcontainsDdLabel.FontSize = 16;
            app.ThesystemconfigurationXcontainsDdLabel.Position = [53 408 310 22];
            app.ThesystemconfigurationXcontainsDdLabel.Text = 'The system configuration X contains:  D, d';

            % Create LowerboundsPanel
            app.LowerboundsPanel = uipanel(app.OptimisationTab);
            app.LowerboundsPanel.TitlePosition = 'centertop';
            app.LowerboundsPanel.Title = 'Lower bounds:';
            app.LowerboundsPanel.FontWeight = 'bold';
            app.LowerboundsPanel.FontSize = 16;
            app.LowerboundsPanel.Position = [98 223 188 154];

            % Create DLabel
            app.DLabel = uilabel(app.LowerboundsPanel);
            app.DLabel.Position = [40 63 28 22];
            app.DLabel.Text = 'D = ';

            % Create mLabel
            app.mLabel = uilabel(app.LowerboundsPanel);
            app.mLabel.Position = [131 63 25 22];
            app.mLabel.Text = 'm';

            % Create lb_D
            app.lb_D = uieditfield(app.LowerboundsPanel, 'numeric');
            app.lb_D.Position = [71 63 49 22];
            app.lb_D.Value = 1;

            % Create mLabel_2
            app.mLabel_2 = uilabel(app.LowerboundsPanel);
            app.mLabel_2.Position = [131 30 25 22];
            app.mLabel_2.Text = 'm';

            % Create lb_d
            app.lb_d = uieditfield(app.LowerboundsPanel, 'numeric');
            app.lb_d.Position = [71 30 49 22];
            app.lb_d.Value = 0.1;

            % Create dLabel
            app.dLabel = uilabel(app.LowerboundsPanel);
            app.dLabel.Position = [40 30 25 22];
            app.dLabel.Text = 'd =';

            % Create UpperboundsPanel
            app.UpperboundsPanel = uipanel(app.OptimisationTab);
            app.UpperboundsPanel.TitlePosition = 'centertop';
            app.UpperboundsPanel.Title = 'Upper bounds:';
            app.UpperboundsPanel.FontWeight = 'bold';
            app.UpperboundsPanel.FontSize = 16;
            app.UpperboundsPanel.Position = [372 223 179 155];

            % Create DLabel_2
            app.DLabel_2 = uilabel(app.UpperboundsPanel);
            app.DLabel_2.Position = [40 64 28 22];
            app.DLabel_2.Text = 'D = ';

            % Create mLabel_5
            app.mLabel_5 = uilabel(app.UpperboundsPanel);
            app.mLabel_5.Position = [131 64 25 22];
            app.mLabel_5.Text = 'm';

            % Create ub_D
            app.ub_D = uieditfield(app.UpperboundsPanel, 'numeric');
            app.ub_D.Position = [71 64 49 22];
            app.ub_D.Value = 15;

            % Create mLabel_6
            app.mLabel_6 = uilabel(app.UpperboundsPanel);
            app.mLabel_6.Position = [131 31 25 22];
            app.mLabel_6.Text = 'm';

            % Create ub_d
            app.ub_d = uieditfield(app.UpperboundsPanel, 'numeric');
            app.ub_d.Position = [71 31 49 22];
            app.ub_d.Value = 5;

            % Create dLabel_2
            app.dLabel_2 = uilabel(app.UpperboundsPanel);
            app.dLabel_2.Position = [40 31 25 22];
            app.dLabel_2.Text = 'd =';

            % Create InitialStatePanel
            app.InitialStatePanel = uipanel(app.OptimisationTab);
            app.InitialStatePanel.TitlePosition = 'centertop';
            app.InitialStatePanel.Title = 'Initial State:';
            app.InitialStatePanel.FontWeight = 'bold';
            app.InitialStatePanel.FontSize = 16;
            app.InitialStatePanel.Position = [642 223 179 154];

            % Create DLabel_3
            app.DLabel_3 = uilabel(app.InitialStatePanel);
            app.DLabel_3.Position = [40 63 28 22];
            app.DLabel_3.Text = 'D = ';

            % Create mLabel_9
            app.mLabel_9 = uilabel(app.InitialStatePanel);
            app.mLabel_9.Position = [131 63 25 22];
            app.mLabel_9.Text = 'm';

            % Create D0
            app.D0 = uieditfield(app.InitialStatePanel, 'numeric');
            app.D0.ValueChangedFcn = createCallbackFcn(app, @D0ValueChanged, true);
            app.D0.Position = [71 63 49 22];
            app.D0.Value = 5;

            % Create mLabel_10
            app.mLabel_10 = uilabel(app.InitialStatePanel);
            app.mLabel_10.Position = [131 30 25 22];
            app.mLabel_10.Text = 'm';

            % Create d0
            app.d0 = uieditfield(app.InitialStatePanel, 'numeric');
            app.d0.ValueChangedFcn = createCallbackFcn(app, @d0ValueChanged, true);
            app.d0.Position = [71 30 49 22];
            app.d0.Value = 1;

            % Create dLabel_3
            app.dLabel_3 = uilabel(app.InitialStatePanel);
            app.dLabel_3.Position = [40 30 25 22];
            app.dLabel_3.Text = 'd =';

            % Create MaxiterationsLabel
            app.MaxiterationsLabel = uilabel(app.OptimisationTab);
            app.MaxiterationsLabel.FontSize = 16;
            app.MaxiterationsLabel.Position = [53 147 110 22];
            app.MaxiterationsLabel.Text = 'Max iterations:';

            % Create MaxfunctionevaluationLabel
            app.MaxfunctionevaluationLabel = uilabel(app.OptimisationTab);
            app.MaxfunctionevaluationLabel.FontSize = 16;
            app.MaxfunctionevaluationLabel.Position = [253 146 178 22];
            app.MaxfunctionevaluationLabel.Text = 'Max function evaluation:';

            % Create MaxIterations
            app.MaxIterations = uieditfield(app.OptimisationTab, 'numeric');
            app.MaxIterations.Position = [165 147 50 22];
            app.MaxIterations.Value = 4000;

            % Create MaxEvals
            app.MaxEvals = uieditfield(app.OptimisationTab, 'numeric');
            app.MaxEvals.Position = [430 146 56 22];
            app.MaxEvals.Value = 600;

            % Create LogIterations
            app.LogIterations = uicheckbox(app.OptimisationTab);
            app.LogIterations.Text = '';
            app.LogIterations.FontSize = 16;
            app.LogIterations.Position = [53 107 26 22];
            app.LogIterations.Value = true;

            % Create Label_7
            app.Label_7 = uilabel(app.OptimisationTab);
            app.Label_7.FontSize = 16;
            app.Label_7.Position = [77 107 102 22];
            app.Label_7.Text = 'Log iterations';

            % Create RunOptimisationButton
            app.RunOptimisationButton = uibutton(app.OptimisationTab, 'push');
            app.RunOptimisationButton.ButtonPushedFcn = createCallbackFcn(app, @RunOptimisationButtonPushed, true);
            app.RunOptimisationButton.FontSize = 16;
            app.RunOptimisationButton.FontWeight = 'bold';
            app.RunOptimisationButton.Position = [604 126 217 43];
            app.RunOptimisationButton.Text = 'Run Optimisation';

            % Create ExitStatus
            app.ExitStatus = uilamp(app.OptimisationTab);
            app.ExitStatus.Enable = 'off';
            app.ExitStatus.Position = [616 136 20 20];
            app.ExitStatus.Color = [0.9412 0.9412 0.9412];

            % Create ExitFlag
            app.ExitFlag = uilabel(app.OptimisationTab);
            app.ExitFlag.Position = [607 93 263 22];
            app.ExitFlag.Text = 'Press "Run Optimisation" ...';

            % Create ResultsTab
            app.ResultsTab = uitab(app.TabGroup);
            app.ResultsTab.Title = 'Results';

            % Create TabGroup2
            app.TabGroup2 = uitabgroup(app.ResultsTab);
            app.TabGroup2.Position = [49 60 435 410];

            % Create OptimalConfigTab
            app.OptimalConfigTab = uitab(app.TabGroup2);
            app.OptimalConfigTab.Title = 'Optimal Config';

            % Create OptiConfigAx
            app.OptiConfigAx = uiaxes(app.OptimalConfigTab);
            title(app.OptiConfigAx, 'Title')
            xlabel(app.OptiConfigAx, 'X')
            ylabel(app.OptiConfigAx, 'Y')
            app.OptiConfigAx.PlotBoxAspectRatio = [1.18575851393189 1 1];
            app.OptiConfigAx.XGrid = 'on';
            app.OptiConfigAx.YGrid = 'on';
            app.OptiConfigAx.Position = [1 1 433 381];

            % Create InitialConfigTab
            app.InitialConfigTab = uitab(app.TabGroup2);
            app.InitialConfigTab.Title = 'Initial Config';

            % Create InitConfigAx
            app.InitConfigAx = uiaxes(app.InitialConfigTab);
            title(app.InitConfigAx, 'Title')
            xlabel(app.InitConfigAx, 'X')
            ylabel(app.InitConfigAx, 'Y')
            app.InitConfigAx.Position = [1 1 433 381];

            % Create XmSliderLabel
            app.XmSliderLabel = uilabel(app.ResultsTab);
            app.XmSliderLabel.HorizontalAlignment = 'right';
            app.XmSliderLabel.Position = [626 447 40 22];
            app.XmSliderLabel.Text = 'X (m)';

            % Create XmSlider
            app.XmSlider = uislider(app.ResultsTab);
            app.XmSlider.Limits = [0 10];
            app.XmSlider.ValueChangedFcn = createCallbackFcn(app, @XmSliderValueChanged, true);
            app.XmSlider.Position = [687 456 150 3];
            app.XmSlider.Value = 4;

            % Create YmSliderLabel
            app.YmSliderLabel = uilabel(app.ResultsTab);
            app.YmSliderLabel.HorizontalAlignment = 'right';
            app.YmSliderLabel.Position = [553 230 40 22];
            app.YmSliderLabel.Text = 'Y (m)';

            % Create YmSlider
            app.YmSlider = uislider(app.ResultsTab);
            app.YmSlider.Limits = [0 10];
            app.YmSlider.Orientation = 'vertical';
            app.YmSlider.ValueChangedFcn = createCallbackFcn(app, @YmSliderValueChanged, true);
            app.YmSlider.Position = [614 239 3 150];
            app.YmSlider.Value = 4;

            % Create PhiKnobLabel
            app.PhiKnobLabel = uilabel(app.ResultsTab);
            app.PhiKnobLabel.HorizontalAlignment = 'center';
            app.PhiKnobLabel.Position = [744 242 36 22];
            app.PhiKnobLabel.Text = 'Phi(°)';

            % Create PhiKnob
            app.PhiKnob = uiknob(app.ResultsTab, 'continuous');
            app.PhiKnob.Limits = [-30 30];
            app.PhiKnob.ValueChangedFcn = createCallbackFcn(app, @PhiKnobValueChanged, true);
            app.PhiKnob.Position = [731 298 60 60];

            % Create OptimaleConfigurationPanel
            app.OptimaleConfigurationPanel = uipanel(app.ResultsTab);
            app.OptimaleConfigurationPanel.TitlePosition = 'centertop';
            app.OptimaleConfigurationPanel.Title = 'Optimale Configuration';
            app.OptimaleConfigurationPanel.FontWeight = 'bold';
            app.OptimaleConfigurationPanel.FontSize = 16;
            app.OptimaleConfigurationPanel.Position = [650 31 195 147];

            % Create DLabel_4
            app.DLabel_4 = uilabel(app.OptimaleConfigurationPanel);
            app.DLabel_4.Position = [40 56 28 22];
            app.DLabel_4.Text = 'D =';

            % Create mLabel_13
            app.mLabel_13 = uilabel(app.OptimaleConfigurationPanel);
            app.mLabel_13.Position = [131 56 25 22];
            app.mLabel_13.Text = 'm';

            % Create X_D
            app.X_D = uieditfield(app.OptimaleConfigurationPanel, 'numeric');
            app.X_D.Position = [71 56 49 22];
            app.X_D.Value = 10;

            % Create mLabel_14
            app.mLabel_14 = uilabel(app.OptimaleConfigurationPanel);
            app.mLabel_14.Position = [131 23 25 22];
            app.mLabel_14.Text = 'm';

            % Create X_d
            app.X_d = uieditfield(app.OptimaleConfigurationPanel, 'numeric');
            app.X_d.Position = [71 23 49 22];
            app.X_d.Value = 1;

            % Create dLabel_4
            app.dLabel_4 = uilabel(app.OptimaleConfigurationPanel);
            app.dLabel_4.Position = [40 23 25 22];
            app.dLabel_4.Text = 'd =';

            % Create CurrentPositionLabel
            app.CurrentPositionLabel = uilabel(app.ResultsTab);
            app.CurrentPositionLabel.Position = [53 40 95 22];
            app.CurrentPositionLabel.Text = 'Current Position:';

            % Create XLabel
            app.XLabel = uilabel(app.ResultsTab);
            app.XLabel.Position = [170 13 25 22];
            app.XLabel.Text = 'X=';

            % Create currentX
            app.currentX = uieditfield(app.ResultsTab, 'numeric');
            app.currentX.Position = [193 12 58 22];

            % Create mLabel_15
            app.mLabel_15 = uilabel(app.ResultsTab);
            app.mLabel_15.Position = [253 12 25 22];
            app.mLabel_15.Text = 'm';

            % Create mLabel_16
            app.mLabel_16 = uilabel(app.ResultsTab);
            app.mLabel_16.Position = [372 12 25 22];
            app.mLabel_16.Text = 'm';

            % Create currentY
            app.currentY = uieditfield(app.ResultsTab, 'numeric');
            app.currentY.Position = [310 12 58 22];

            % Create YLabel
            app.YLabel = uilabel(app.ResultsTab);
            app.YLabel.Position = [287 12 25 22];
            app.YLabel.Text = 'Y=';

            % Create Label_8
            app.Label_8 = uilabel(app.ResultsTab);
            app.Label_8.Position = [492 18 25 22];
            app.Label_8.Text = '°';

            % Create currentPhi
            app.currentPhi = uieditfield(app.ResultsTab, 'numeric');
            app.currentPhi.Position = [433 12 58 22];

            % Create PhiLabel
            app.PhiLabel = uilabel(app.ResultsTab);
            app.PhiLabel.Position = [406 12 30 22];
            app.PhiLabel.Text = 'Phi=';

            % Create GoToButton
            app.GoToButton = uibutton(app.ResultsTab, 'push');
            app.GoToButton.ButtonPushedFcn = createCallbackFcn(app, @GoToButtonPushed, true);
            app.GoToButton.Interruptible = 'off';
            app.GoToButton.Position = [53 13 96 22];
            app.GoToButton.Text = 'Go To';

            % Create OptimisationresultsLabel
            app.OptimisationresultsLabel = uilabel(app.ResultsTab);
            app.OptimisationresultsLabel.FontSize = 24;
            app.OptimisationresultsLabel.FontWeight = 'bold';
            app.OptimisationresultsLabel.Position = [53 476 291 30];
            app.OptimisationresultsLabel.Text = 'Optimisation results';

            % Create ValidationTab
            app.ValidationTab = uitab(app.TabGroup);
            app.ValidationTab.Title = 'Validation';

            % Create ValidationAxis
            app.ValidationAxis = uiaxes(app.ValidationTab);
            title(app.ValidationAxis, '')
            xlabel(app.ValidationAxis, 'X')
            ylabel(app.ValidationAxis, 'Y')
            app.ValidationAxis.Position = [334 20 525 487];

            % Create AnalyseSolutionButton
            app.AnalyseSolutionButton = uibutton(app.ValidationTab, 'push');
            app.AnalyseSolutionButton.ButtonPushedFcn = createCallbackFcn(app, @AnalyseSolutionButtonPushed, true);
            app.AnalyseSolutionButton.BusyAction = 'cancel';
            app.AnalyseSolutionButton.Interruptible = 'off';
            app.AnalyseSolutionButton.Position = [67 217 201 70];
            app.AnalyseSolutionButton.Text = 'Analyse Solution';

            % Create ThistakesafewminutesPleasewaitLabel
            app.ThistakesafewminutesPleasewaitLabel = uilabel(app.ValidationTab);
            app.ThistakesafewminutesPleasewaitLabel.Position = [67 181 210 22];
            app.ThistakesafewminutesPleasewaitLabel.Text = 'This takes a few minutes. Please wait.';

            % Show the figure after all components are created
            app.RPRAnalysistoolByAMazouniHMartinezUIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = interface_exported

            % Create UIFigure and components
            createComponents(app)

            % Register the app with App Designer
            registerApp(app, app.RPRAnalysistoolByAMazouniHMartinezUIFigure)

            % Execute the startup function
            runStartupFcn(app, @startupFcn)

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.RPRAnalysistoolByAMazouniHMartinezUIFigure)
        end
    end
end