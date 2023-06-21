function compileMiceControlRipples

expPath = 'Z:\Homes\zutshi01\Recordings\CA1_silencing\';
pathToSessionsAll = {'IZ50\IZ50_230501_sess1','IZ50\IZ50_230502_sess2','IZ50\IZ50_230503_sess3','IZ50\IZ50_230504_sess4'...
    'IZ51\IZ51_230501_sess1','IZ51\IZ51_230502_sess2','IZ51\IZ51_230503_sess3','IZ51\IZ51_230504_sess4'};
%{'IZ25\IZ25_0um_201012_sess11','IZ25\IZ25_0um_201013_sess12'};    
    %,'IZ12\IZ12_288um_200122_sess3','IZ12\IZ12_288um_200127_sess4',...
   %  'IZ13\IZ13_216um_200304_sess3','IZ13\IZ13_360um_200306_sess4','IZ16\IZ16_144um_200616_sess1','IZ25\IZ25_0um_201012_sess11','IZ25\IZ25_0um_201013_sess12'};

figure

for kk = 1%:2
    for i=1:size(pathToSessionsAll,2)
        %cd(strcat(expPath,pathToSessionsAll{i}))
        cd(strcat(expPath,pathToSessionsAll{i}))
        file = dir(['*pulses.events.mat']);
        load(file(1).name)

        file = dir(['*sessioninfo.mat']);
        load(file(1).name)

        file = dir(['*ripples.events.mat']);
        load(file(1).name)    

        pulTr = (pulses.stimComb==kk);
        events = pulses.intsPeriods(:,pulTr)';

        ripple_pre = InIntervals(ripples.peaks,events-10);
        ripple_prestim = InIntervals(ripples.peaks,events-5);
        ripple_post = InIntervals(ripples.peaks,events);
        ripple_poststim = InIntervals(ripples.peaks,events+10);

        rate(i,:)= [sum(ripple_pre) sum(ripple_prestim) sum(ripple_post) sum(ripple_poststim)]./(numel(events)*5);

    end    
    subplot(1,2,kk)
    data{1} = rate(:,1);
    data{2} = rate(:,2);
    data{3} = rate(:,3);
    data{4} = rate(:,4);

    stats = groupStats(data,[],'repeatedMeasures',true,'plotType','boxLinesSEM','inAxis',true);
    ylim([0 0.3])
end
end