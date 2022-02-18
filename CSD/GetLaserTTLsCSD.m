function [ON_ts, OFF_ts, RecStart, RecEnd] = GetLaserTTLsCSD(folder_name)

%This function finds TTL ON and OFF events of light stimulation
%Normally, ON: (0x0084), OFF : (0x0080)
%By I.Z. Aug 2015


% Load events file
file_name = strcat(folder_name,filesep, 'Events.nev');
[Events.timeStamps Events.eventString] = Nlx2MatEV(file_name,[1 0 0 0 1],0,1,[]);

% find unique events
X = unique(Events.eventString);
tmp = strfind(X, 'TTL Input');
InputStrs = {};
cnt = 1;
for hih = 1:length(tmp), 
    if tmp{hih} == 1,
       InputStrs{cnt} = X(hih);
       cnt = cnt+1;
    end
end
No_inputStrings = 0;
if length(InputStrs) < 2,
    ON_ts = [];
    OFF_ts = [];
    No_inputStrings = 1;
    return
end

RecStart = Events.timeStamps(1);
Rec_End = Events.timeStamps(size(Events.timeStamps));    
RecEnd = Rec_End(2);
% possible ON strings
T_ON1 = strcmp(X, 'TTL Input on AcqSystem1_0 board 0 port 1 value (0x0084).'); if sum(T_ON1) == 1, ON_string = 'TTL Input on AcqSystem1_0 board 0 port 1 value (0x0084).'; end
T_ON2 = strcmp(X, 'TTL Input on AcqSystem1_0 board 0 port 1 value (0x00F8).'); if sum(T_ON2) == 1, ON_string = 'TTL Input on AcqSystem1_0 board 0 port 1 value (0x00F8).'; end
T_ON3 = strcmp(X, 'TTL Input on AcqSystem1_0 board 0 port 1 value (0x0000).'); if sum(T_ON3) == 1, ON_string = 'TTL Input on AcqSystem1_0 board 0 port 1 value (0x0000).'; end

% possible OFF strings
T_OFF1 = strcmp(X, 'TTL Input on AcqSystem1_0 board 0 port 1 value (0x0080).'); if sum(T_OFF1) == 1, OFF_string = 'TTL Input on AcqSystem1_0 board 0 port 1 value (0x0080).'; end
T_OFF2 = strcmp(X, 'TTL Input on AcqSystem1_0 board 0 port 1 value (0x00FC).'); if sum(T_OFF2) == 1, OFF_string = 'TTL Input on AcqSystem1_0 board 0 port 1 value (0x00FC).'; end
T_OFF3 = strcmp(X, 'TTL Input on AcqSystem1_0 board 0 port 1 value (0x0004).'); if sum(T_OFF3) == 1, OFF_string = 'TTL Input on AcqSystem1_0 board 0 port 1 value (0x0004).'; end

if sum(sum([T_ON1 T_ON2 T_ON3 T_OFF1 T_OFF2 T_OFF3])) > 2, error('Multiple ON or OFF strings detected in root.event'), end
if sum(sum([T_ON1 T_ON2 T_ON3 T_OFF1 T_OFF2 T_OFF3])) < 2, display('No stimulation detected.  May need to include additional ON/OFF strings in GetStartStopStim'),No_inputStrings = 1; return,  end


% get logical for ON events
ON_eventIND = strcmp(Events.eventString, ON_string);
ON_IND = find(ON_eventIND);

% get logical for OFF events
OFF_eventIND = strcmp(Events.eventString, OFF_string);
OFF_IND = find(OFF_eventIND);

% get ON and OFF timestamps
ON_ts = Events.timeStamps(ON_IND);
OFF_ts = Events.timeStamps(OFF_IND);

end