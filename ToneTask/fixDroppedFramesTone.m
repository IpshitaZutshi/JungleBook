figure
subplot(2,7,[1 7])
plot(tracking.timestamps(1:length(tracking.position.y)),tracking.position.y)

startReward = behavTrials.timestamps(:,1);
endReward  = behavTrials.timestamps(:,2);
for ii = 1:length(startReward)
	[~,idxStart(ii)] = min(abs(tracking.timestamps-startReward(ii)));
end

for ii = 1:length(startReward)
	[~,idxEnd(ii)] = min(abs(tracking.timestamps-endReward(ii)));
end
hold on;plot(tracking.timestamps(idxStart),tracking.position.y(idxStart),'r.')
hold on;plot(tracking.timestamps(idxEnd),tracking.position.y(idxEnd),'b.')

lickLoc = behavTrials.lickLoc;
for ii = 1:6
    subplot(2,7,7+ii)
    posLick = tracking.position.y(idxEnd(lickLoc==(ii-1)));
    histogram(posLick,0:1:120)
    title(num2str(median(posLick)))
end
subplot(2,7,14)
posLick = tracking.position.y(idxStart);
histogram(posLick,0:1:120)
title(num2str(median(posLick)))

% Start of every trial, find the FIRST index where it reaches 8, versus the
% actual index. Plot the distribution of this difference across trials - at
% some point, when frames are dropped, you will see this difference
% increase. 
% 30 second 

% If the difference > a certain amount, try to correct it by REMOVING those
% timestamps. Repeat for following trials, after fixing the current trial
% % to see if maybe frames were only dropped once and can fix the error. 
% 
% figure
% plot(tracking.position.y)
% idxtoMatch =72232;
% idxProbTrial  = idxStart(90);
% idxDropped = idxProbTrial-idxtoMatch-1;
% % % % Move the position values back by idxDropped units and fill those with
% % % % nans
% InsertArr =[];
% InsertArr(1:idxDropped,1) = nan;
% x = [];y = [];vx = [];vy=[]; v=[];
% x = [tracking.position.x(1:idxtoMatch);InsertArr;tracking.position.x(idxtoMatch+1:end)];
% y = [tracking.position.y(1:idxtoMatch);InsertArr;tracking.position.y(idxtoMatch+1:end)];
% vx = [tracking.position.vx(1:idxtoMatch);InsertArr;tracking.position.vx(idxtoMatch+1:end)];
% vy = [tracking.position.vy(1:idxtoMatch);InsertArr;tracking.position.vx(idxtoMatch+1:end)];
% v = [tracking.position.v(1:idxtoMatch);InsertArr;tracking.position.v(idxtoMatch+1:end)];
% % 
% tracking.position.x = x;%(1:length(tracking.timestamps));
% tracking.position.y = y;%(1:length(tracking.timestamps));
% tracking.position.vx = vx;%(1:length(tracking.timestamps));
% tracking.position.vy = vy;%(1:length(tracking.timestamps));
% tracking.position.v = v;%(1:length(tracking.timestamps));

%tracking.timestamps = tracking.timestamps(1:length(tracking.position.x));
% tracking.position.x = x(1:length(tracking.timestamps));
% tracking.position.y = y(1:length(tracking.timestamps));
% tracking.position.vx = vx(1:length(tracking.timestamps));
% tracking.position.vy = vy(1:length(tracking.timestamps));
% tracking.position.v = v(1:length(tracking.timestamps));