function checkPGAMfits

sess= {'IZ39\Final\IZ39_220622_sess8','IZ39\Final\IZ39_220624_sess10','IZ39\Final\IZ39_220629_sess12',...
    'IZ39\Final\IZ39_220702_sess14','IZ39\Final\IZ39_220714_sess18',...
    'IZ39\Final\IZ39_220705_sess16','IZ39\Final\IZ39_220707_sess17',...   
    'IZ40\Final\IZ40_220705_sess15','IZ40\Final\IZ40_220707_sess16',...
    'IZ40\Final\IZ40_220708_sess17','IZ40\Final\IZ40_220714_sess18',...
    'IZ43\Final\IZ43_220826_sess2','IZ43\Final\IZ43_220828_sess4',...
    'IZ43\Final\IZ43_220830_sess6','IZ43\Final\IZ43_220901_sess8',...
    'IZ43\Final\IZ43_220911_sess9','IZ43\Final\IZ43_220913_sess11','IZ43\Final\IZ43_220919_sess14',...
    'IZ43\Final\IZ43_220915_sess13','IZ43\Final\IZ43_220920_sess15',...    
    'IZ44\Final\IZ44_220827_sess4', 'IZ44\Final\IZ44_220828_sess5',...
    'IZ44\Final\IZ44_220829_sess6','IZ44\Final\IZ44_220830_sess7',...
    'IZ44\Final\IZ44_220912_sess10','IZ44\Final\IZ44_220913_sess11','IZ44\Final\IZ44_220919_sess14',...
    'IZ44\Final\IZ44_220915_sess13','IZ44\Final\IZ44_220920_sess15',...
    }; 

expPath = 'Z:\Homes\zutshi01\Recordings\Auditory_Task\';
pgamPATH = 'C:\Data\PGAMAnalysis\results\new_vars2';
%'Z:\Buzsakilabspace\LabShare\AthinaApostolelli\PGAM\results\new_vars2';
saveloc = 'Z:\Homes\zutshi01\Recordings\Auditory_Task\Compiled\PGAM fits_new\';

numcol = 10*2;
numrows = 24;

for ii = 1:length(sess)
    
    fig2  = figure;
    set(fig2,'Renderer','painters')
    set(fig2,'Color','w')
    set(fig2,'Position',[1 41 1920 970]);

    cd(strcat(expPath,sess{ii}))    
    file = dir(['*cell_metrics.cellinfo.mat']);
    load(file(1).name);        
    
    %Find 15 random cells 
    idx = find(cell_metrics.firingRate>1 & cell_metrics.firingRate<5);
    
    randNum = randsample(idx,10);
    
    for rr = 1:length(randNum)
        try
            cellNum =randNum(rr);
            file = dir(['*.spikeData.cellinfo.mat']);
            load(file.name);
            file = dir(['*.Tracking.Behavior.mat']);
            load(file(1).name);
            file = dir(['*TrialBehavior.Behavior.mat']);
            load(file.name);
            file = dir(['*.rateMapsAvg.cellinfo.mat']);
            load(file.name);
            plotExampleCell(1,((rr-1)*2)+1, spikeData, tracking, behavTrials, firingMaps,fig2,numrows,numcol,cellNum)

            pgamloc = strcat(pgamPATH,sess{ii}(11:end),'\spatial_neuron_',num2str(cellNum),'_fit_k-fold.csv');
            results = readtable(pgamloc);  

            subplot(numrows,numcol,[15*numcol+((rr-1)*2)+1,15*numcol+((rr-1)*2)+2])
            %First ylin
            pval = table2array(results(3,10));
            rates = table2array(results(3,13));
            modelrates = readRates(rates);
            rates = table2array(results(3,14));
            rawrates = readRates(rates);
            rates = table2array(results(3,12));
            x = readRates(rates);
            %x = x(2:end);
            plot(rawrates,'Color',[0.5 0.5 0.5])
            hold on
            if pval <0.05 && pval >0
                col = 'r';
            elseif pval ==0
                col = 'b';
            else
                col = 'k';
            end
            plot(modelrates,'Color',col)
            title('Ylin')
            mut_info = table2array(results(3,11));
            title(num2str(mut_info))


            subplot(numrows,numcol,[15*numcol+((rr-1)*2)+1+2*numcol,15*numcol+((rr-1)*2)+2+2*numcol])

            pval = table2array(results(1,10));
            rates = table2array(results(1,13));
            modelrates = readRates(rates);
            rates = table2array(results(1,14));
            rawrates = readRates(rates);
            rates = table2array(results(1,12));
            x = readRates(rates);
            x = x(2:end);
            plot(rawrates,'Color',[0.5 0.5 0.5])
            hold on
            if pval <0.05 && pval >0
                col = 'r';
            elseif pval ==0
                col = 'b';
            else
                col = 'k';
            end
            plot(modelrates,'Color',col)
            mut_info = table2array(results(1,11));
            title(num2str(mut_info))


            subplot(numrows,numcol,[15*numcol+((rr-1)*2)+1+4*numcol,15*numcol+((rr-1)*2)+2+4*numcol])
            %First ylin

            pval = table2array(results(4,10));
            rates = table2array(results(4,13));
            modelrates = readRates(rates);
            rates = table2array(results(4,14));
            rawrates = readRates(rates);
            rates = table2array(results(4,12));
            x = readRates(rates);
            %x = x(2:end);
            plot(rawrates,'Color',[0.5 0.5 0.5])
            hold on
            if pval <0.05 && pval >0
                col = 'r';
            elseif pval ==0
                col = 'b';
            else
                col = 'k';
            end
            plot(modelrates,'Color',col)
            mut_info = table2array(results(4,11));
            title(num2str(mut_info))


            subplot(numrows,numcol,[15*numcol+((rr-1)*2)+1+6*numcol,15*numcol+((rr-1)*2)+2+6*numcol])
            %First ylin

            pval = table2array(results(10,10));
            rates = table2array(results(10,13));
            modelrates = readRates(rates);
            rates = table2array(results(10,14));
            rawrates = readRates(rates);
            rates = table2array(results(10,12));
            x = readRates(rates);
            %x = x(2:end);
            plot(rawrates,'Color',[0.5 0.5 0.5])
            hold on
            if pval <0.05 && pval >0
                col = 'r';
            elseif pval ==0
                col = 'b';
            else
                col = 'k';
            end
            plot(modelrates,'Color',col)
            title('Licks')
            mut_info = table2array(results(10,11));
            title(num2str(mut_info))     
            
            subplot(numrows,numcol,[15*numcol+((rr-1)*2)+1+8*numcol,15*numcol+((rr-1)*2)+2+8*numcol])
            title(strcat('Cell num ',num2str(cellNum)))
            axis off
        catch 
        end
    end
    
    saveas(fig2,strcat(saveloc,sess{ii}(12:end),'.png'));
%    saveas(gcf,strcat(expPath,'Compiled\returnpopulationMaps.eps'),'epsc');
    saveas(fig2,strcat(saveloc,sess{ii}(12:end),'.fig'));
    close all
    

end
end

function rawrates = readRates(rates)

rates = rates{1};
stringData = strrep(rates, '[', '');
stringData = strrep(stringData, ']', '');
stringData = strrep(stringData, '↵', '');
rawrates = regexp(stringData, '\s+', 'split');        
rawrates = str2double(rawrates);
end