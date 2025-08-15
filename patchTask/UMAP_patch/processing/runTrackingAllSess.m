%
% USAGE
%   runs tracking behavior for all folders
%
% INPUTS 
%    start in top directory of mouse
%
%    =========================================================================


direc = pwd;

folders = dir(direc);  % list all files/folders in direc

for i = 1:length(folders)
    if folders(i).isdir && ~ismember(folders(i).name, {'.', '..'}) % only folders, ignore '.' and '..'
        folderName = folders(i).name;
        if contains(folderName, 'sess')
            currFolder = fullfile(direc, folderName);
            fprintf('Processing folder: %s\n', currFolder);
            cd(currFolder);
            basepath = pwd;

            if ~isempty (dir(fullfile(basepath, '*session.mat')))
                % sleep score
                badChannels = [42 48 56:59 61 70:72]; %N17
                SleepScoreMaster(pwd,'stickytrigger',true,'rejectChannels',badChannels); % try to sleep score

                SleepScoreMaster(pwd,'stickytrigger',true,'rejectChannels',badChannels); % try to sleep score


                %getPatchTracking()
                % get tracking
                % if isempty (dir(fullfile(basepath, '*Tracking.Behavior.mat')))
                %     getPatchTracking() 
                % end
                % 
                % get just ripples
                % if isempty (dir(fullfile(basepath, '*ripples.events.mat')))
                %     if ~isempty (dir(fullfile(basepath, '*.lfp')))
                %         % Extract sharp wave ripples
                %         pyrCh = 121; % n17
                %         noiseCh = 111;
                %         [ripples] = bz_FindRipples(pwd,pyrCh,'noise',noiseCh,'savemat',true,'durations',[30 100],'passband',[130 200]);
                %     end
                % end

                %  % get lfp and ripples
                %  if isempty (dir(fullfile(basepath, '*.lfp')))
                %      %% 2. Extract LFP
                %     [sessionInfo] = bz_getSessionInfo(pwd, 'noPrompts', true); 
                %     sessionInfo.rates.lfp = 1250;  
                %     save(strcat(sessionInfo.session.name,'.sessionInfo.mat'),'sessionInfo');
                %     if isempty(dir('*.lfp'))
                %         try 
                %             bz_LFPfromDat(pwd,'outFs',1250); % generating lfp
                %         catch
                %             disp('Problems with bz_LFPfromDat, resampling...');
                %             ResampleBinary(strcat(sessionInfo.session.name,'.dat'),strcat(sessionInfo.session.name,'.lfp'),...
                %                 sessionInfo.nChannels,1,sessionInfo.rates.wideband/sessionInfo.rates.lfp);
                %         end
                %     end
                %     % Extract sharp wave ripples
                %     pyrCh = 121; % n17
                %     noiseCh = 111;
                %     [ripples] = bz_FindRipples(pwd,pyrCh,'noise',noiseCh,'savemat',true,'durations',[30 100],'passband',[130 200]);
                % end
            end

            cd(direc);
        end
    end
end


