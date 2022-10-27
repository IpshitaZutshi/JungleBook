classdef Inspector_CheckRipples < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        PlaceField                    matlab.ui.Figure
        NextEventButton               matlab.ui.control.Button
        LastEventButton               matlab.ui.control.Button
        Label                         matlab.ui.control.Label
        savetofileButton              matlab.ui.control.Button
        UIAxes                        matlab.ui.control.UIAxes
        RemovecurrenteventButton      matlab.ui.control.Button
        RestorealldeletedeventButton  matlab.ui.control.Button
        GotoeventEditFieldLabel       matlab.ui.control.Label
        GotoeventEditField            matlab.ui.control.NumericEditField
        ExtendEditFieldLabel          matlab.ui.control.Label
        ExtendEditField               matlab.ui.control.NumericEditField
        LFPAmplitudePanel             matlab.ui.container.Panel
        increase                      matlab.ui.control.Button
        decrease                      matlab.ui.control.Button
    end

    
    properties (Access = private)
        ip; LFPs; Ripples; xmlinfo; Ripple_bak; basepath;
        channelxy; rmlist; iplist; ext; lfpamp;
    end
   
    methods (Access = private)
        
        function plotim(app)
            textLabel = sprintf('Ripple # %d (total %d); Duration: %.2f ms.', app.ip, size(app.Ripples.peaks,1),diff(app.Ripples.timestamps(app.ip,:))*1000);
            app.Label.Text = textLabel;
            app.GotoeventEditField.Value = app.ip;
            app.ExtendEditField.Value = app.ext;
            
            cla(app.UIAxes); extend = app.ext*diff(app.Ripples.timestamps(app.ip,:));
            epoch = app.Ripples.epochsidx(app.ip,1:2);
            plotrange = epoch + round([-extend,extend]*app.LFPs.samplingRate);
            ct = diff(plotrange)/app.LFPs.samplingRate;
            for i = 1:size(app.LFPs.data,2)
                ic = app.LFPs.channels == app.channelxy(i,3);
                idx = plotrange(1):plotrange(2);
                time = app.channelxy(i,1)*1.1*ct + app.LFPs.timestamps(idx);
                data = app.channelxy(i,2)*app.lfpamp + app.LFPs.data(idx,ic);
                plot(app.UIAxes,time,data,'k'); hold(app.UIAxes,'on');
            end
            for is = 1:length(app.xmlinfo.AnatGrps)
                time = is*1.1*ct + app.LFPs.timestamps(idx);
                if (app.ext > 0.01)
                    icx = round([extend*app.LFPs.samplingRate,length(time)-extend*app.LFPs.samplingRate]);
                    pch = patch(app.UIAxes,[time(icx(1)),time(icx(2)),time(icx(2)),time(icx(1))],[0,0,1.1*app.lfpamp*length(app.xmlinfo.AnatGrps(1).Channels),1.1*app.lfpamp*length(app.xmlinfo.AnatGrps(1).Channels)],'g');
                    pch.FaceAlpha = 0.2;  hold(app.UIAxes,'on'); pch.EdgeColor = 'none';
                end
                plot(app.UIAxes,(is*1.1*ct+app.Ripples.peaks(app.ip))*[1,1],[0,1.1*app.lfpamp*length(app.xmlinfo.AnatGrps(1).Channels)],'b');
            end
            xrange = app.LFPs.timestamps(idx(1))+[0.8*ct,(length(app.xmlinfo.AnatGrps)+1)*ct*1.1+0.1*ct];
            yrange = [0,1.05*app.lfpamp*length(app.xmlinfo.AnatGrps(1).Channels)];
            app.UIAxes.XLim = xrange;
            app.UIAxes.YLim = yrange;
            app.UIAxes.XTick = [];
            app.UIAxes.YTick = [];    
        end
    end
 

    % Callbacks that handle component events
    methods (Access = private)

        % Code that executes after component creation
        function startupFcn(app, LFPs, Ripples, basepath)
            app.LFPs = LFPs;
            app.Ripples = Ripples;
            app.Ripple_bak = Ripples;
            app.basepath = basepath;
            
            app.xmlinfo = LoadParameters(basepath);
            app.LFPs.data = zscore(double(LFPs.data),[],'all');
            app.channelxy = zeros(size(app.LFPs.data,2),3);
            icc = 0;
            for ishk = 1:length(app.xmlinfo.AnatGrps)
                for ic = 1:length(app.xmlinfo.AnatGrps(ishk).Channels)
                    icc = icc + 1;
                    app.channelxy(icc,:) = [ishk,length(app.xmlinfo.AnatGrps(ishk).Channels)-ic+1,app.xmlinfo.AnatGrps(ishk).Channels(ic)];
                end
            end
          
            app.ip = 1; app.ext = 0.5; app.lfpamp = 3;
            app.iplist = 1:size(app.Ripples.timestamps,1);
            app.rmlist = [];
            
            plotim(app)
        end

        % Button pushed function: NextEventButton
        function NextEventButtonPushed(app, event)
            app.ip = app.ip+1;
            if app.ip > length(app.iplist)
                app.ip = 1;
            end
            plotim(app);
        end

        % Button pushed function: LastEventButton
        function LastEventButtonPushed(app, event)
            app.ip = app.ip-1;
            if app.ip < 1
                app.ip = length(app.iplist);
            end
            plotim(app);
        end

        % Button pushed function: savetofileButton
        function savetofileButtonPushed(app, event)
            ripples = app.Ripples;
            save([app.basepath,'\ripple.mat'],'ripples');
        end

        % Button pushed function: RemovecurrenteventButton
        function RemovecurrenteventButtonPushed(app, event)
            app.rmlist = [app.rmlist;app.ip];
            app.Ripples.timestamps(app.ip,:) = [];
            app.Ripples.peaks(app.ip,:) = [];
            app.Ripples.epochsidx(app.ip,:) = [];
            plotim(app);
        end

        % Value changed function: GotoeventEditField
        function GotoeventEditFieldValueChanged(app, event)
            value = app.GotoeventEditField.Value;
            app.ip = value;
            plotim(app)
        end

        % Button pushed function: RestorealldeletedeventButton
        function RestorealldeletedeventButtonPushed(app, event)
            app.Ripples = app.Ripple_bak;
            app.ip= 1;
            plotim(app)
        end

        % Value changed function: ExtendEditField
        function ExtendEditFieldValueChanged(app, event)
            value = app.ExtendEditField.Value;
            app.ext = value;
            plotim(app);
        end

        % Button pushed function: increase
        function increaseButtonPushed(app, event)
            app.lfpamp = app.lfpamp * 0.8;
            if app.lfpamp < 0.01, app.lfpamp = 0.01; end
            plotim(app);
        end

        % Button pushed function: decrease
        function decreaseButtonPushed(app, event)
            app.lfpamp = app.lfpamp * 1.2;
            if app.lfpamp > 20, app.lfpamp = 20; end
            plotim(app);
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create PlaceField and hide until all components are created
            app.PlaceField = uifigure('Visible', 'off');
            app.PlaceField.Position = [100 100 1000 800];
            app.PlaceField.Name = 'UI Figure';

            % Create NextEventButton
            app.NextEventButton = uibutton(app.PlaceField, 'push');
            app.NextEventButton.ButtonPushedFcn = createCallbackFcn(app, @NextEventButtonPushed, true);
            app.NextEventButton.Position = [802 581 100 22];
            app.NextEventButton.Text = 'Next Event';

            % Create LastEventButton
            app.LastEventButton = uibutton(app.PlaceField, 'push');
            app.LastEventButton.ButtonPushedFcn = createCallbackFcn(app, @LastEventButtonPushed, true);
            app.LastEventButton.Position = [802 629 100 22];
            app.LastEventButton.Text = 'Last Event';

            % Create Label
            app.Label = uilabel(app.PlaceField);
            app.Label.HorizontalAlignment = 'center';
            app.Label.FontSize = 22;
            app.Label.FontWeight = 'bold';
            app.Label.FontColor = [0.502 0.502 0.502];
            app.Label.Position = [246 731 510 27];

            % Create savetofileButton
            app.savetofileButton = uibutton(app.PlaceField, 'push');
            app.savetofileButton.ButtonPushedFcn = createCallbackFcn(app, @savetofileButtonPushed, true);
            app.savetofileButton.Position = [815 96 73 40];
            app.savetofileButton.Text = 'save to file';

            % Create UIAxes
            app.UIAxes = uiaxes(app.PlaceField);
            title(app.UIAxes, 'Title')
            xlabel(app.UIAxes, '')
            ylabel(app.UIAxes, '')
            app.UIAxes.Position = [65 54 690 653];

            % Create RemovecurrenteventButton
            app.RemovecurrenteventButton = uibutton(app.PlaceField, 'push');
            app.RemovecurrenteventButton.ButtonPushedFcn = createCallbackFcn(app, @RemovecurrenteventButtonPushed, true);
            app.RemovecurrenteventButton.Position = [802 242 100 36];
            app.RemovecurrenteventButton.Text = {'Remove current'; ' event'};

            % Create RestorealldeletedeventButton
            app.RestorealldeletedeventButton = uibutton(app.PlaceField, 'push');
            app.RestorealldeletedeventButton.ButtonPushedFcn = createCallbackFcn(app, @RestorealldeletedeventButtonPushed, true);
            app.RestorealldeletedeventButton.Position = [802 190 100 36];
            app.RestorealldeletedeventButton.Text = {'Restore all'; 'deleted event'};

            % Create GotoeventEditFieldLabel
            app.GotoeventEditFieldLabel = uilabel(app.PlaceField);
            app.GotoeventEditFieldLabel.HorizontalAlignment = 'right';
            app.GotoeventEditFieldLabel.Position = [811 535 74 22];
            app.GotoeventEditFieldLabel.Text = 'Goto event #';

            % Create GotoeventEditField
            app.GotoeventEditField = uieditfield(app.PlaceField, 'numeric');
            app.GotoeventEditField.ValueDisplayFormat = '%d';
            app.GotoeventEditField.ValueChangedFcn = createCallbackFcn(app, @GotoeventEditFieldValueChanged, true);
            app.GotoeventEditField.HorizontalAlignment = 'center';
            app.GotoeventEditField.Position = [802 514 100 22];

            % Create ExtendEditFieldLabel
            app.ExtendEditFieldLabel = uilabel(app.PlaceField);
            app.ExtendEditFieldLabel.HorizontalAlignment = 'right';
            app.ExtendEditFieldLabel.Position = [826 443 46 22];
            app.ExtendEditFieldLabel.Text = 'Extend ';

            % Create ExtendEditField
            app.ExtendEditField = uieditfield(app.PlaceField, 'numeric');
            app.ExtendEditField.ValueChangedFcn = createCallbackFcn(app, @ExtendEditFieldValueChanged, true);
            app.ExtendEditField.HorizontalAlignment = 'center';
            app.ExtendEditField.Position = [811 422 77 22];

            % Create LFPAmplitudePanel
            app.LFPAmplitudePanel = uipanel(app.PlaceField);
            app.LFPAmplitudePanel.TitlePosition = 'centertop';
            app.LFPAmplitudePanel.Title = 'LFP Amplitude';
            app.LFPAmplitudePanel.Position = [792 319 120 85];

            % Create increase
            app.increase = uibutton(app.LFPAmplitudePanel, 'push');
            app.increase.ButtonPushedFcn = createCallbackFcn(app, @increaseButtonPushed, true);
            app.increase.FontSize = 24;
            app.increase.Position = [6 16 45 37];
            app.increase.Text = '+';

            % Create decrease
            app.decrease = uibutton(app.LFPAmplitudePanel, 'push');
            app.decrease.ButtonPushedFcn = createCallbackFcn(app, @decreaseButtonPushed, true);
            app.decrease.FontSize = 24;
            app.decrease.Position = [65 16 45 37];
            app.decrease.Text = '-';

            % Show the figure after all components are created
            app.PlaceField.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = Inspector_CheckRipples(varargin)

            % Create UIFigure and components
            createComponents(app)

            % Register the app with App Designer
            registerApp(app, app.PlaceField)

            % Execute the startup function
            runStartupFcn(app, @(app)startupFcn(app, varargin{:}))

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.PlaceField)
        end
    end
end