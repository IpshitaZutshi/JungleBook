function Summary = compileMice2DPlaceFields

sess= {'IZ39\Final\IZ39_220622_sess8','IZ39\Final\IZ39_220624_sess10','IZ39\Final\IZ39_220629_sess12',...
    'IZ39\Final\IZ39_220702_sess14','IZ39\Final\IZ39_220714_sess18',...
    'IZ39\Final\IZ39_220705_sess16','IZ39\Final\IZ39_220707_sess17',...   23
    'IZ40\Final\IZ40_220705_sess15','IZ40\Final\IZ40_220707_sess16',...
    'IZ40\Final\IZ40_220708_sess17','IZ40\Final\IZ40_220714_sess18',...27
    'IZ43\Final\IZ43_220826_sess2','IZ43\Final\IZ43_220828_sess4',...
    'IZ43\Final\IZ43_220830_sess6','IZ43\Final\IZ43_220901_sess8',...
    'IZ43\Final\IZ43_220911_sess9','IZ43\Final\IZ43_220913_sess11','IZ43\Final\IZ43_220919_sess14',...
    'IZ43\Final\IZ43_220915_sess13','IZ43\Final\IZ43_220920_sess15',...    36
    'IZ44\Final\IZ44_220827_sess4', 'IZ44\Final\IZ44_220828_sess5',...
    'IZ44\Final\IZ44_220829_sess6','IZ44\Final\IZ44_220830_sess7',...
    'IZ44\Final\IZ44_220912_sess10','IZ44\Final\IZ44_220913_sess11','IZ44\Final\IZ44_220919_sess14',...
    'IZ44\Final\IZ44_220915_sess13','IZ44\Final\IZ44_220920_sess15',... 45
    'IZ47\Final\IZ47_230707_sess24',...
    'IZ47\Final\IZ47_230710_sess25','IZ47\Final\IZ47_230712_sess27',...49
    'IZ48\Final\IZ48_230628_sess17','IZ48\Final\IZ48_230703_sess21',...
    'IZ48\Final\IZ48_230705_sess22','IZ48\Final\IZ48_230714_sess28'}; 

expPath = 'Z:\Homes\zutshi01\Recordings\Auditory_Task\';

Summary.linCorr = [];
Summary.toneNoToneCorr = [];
Summary.linlinEndCorr = [];
Summary.tonelinEndCorr = [];

for ii = 1:length(sess)
    %% Load files
    cd(strcat(expPath,sess{ii}))    
    file = dir(['*.rateMapsAvg2D.cellinfo.mat']);
    Maps2D = load(file(1).name);
    
    file = dir(['*.rateMapsAvg.cellinfo.mat']);
    load(file(1).name);   
    
    file = dir(['*cell_metrics.cellinfo.mat']);
    load(file(1).name);
    
    for kk=1:length(cell_metrics.UID)
        if strcmp(cell_metrics.putativeCellType{kk},'Pyramidal Cell') == 1 
            cellType = 1;
        else
            cellType = 0;
        end

        Field_Info = detectFields(firingMaps.forward.rateMaps{kk}{1});
        if isempty(Field_Info)
            linField = 0;
        else 
            linField = 1;
        end

        Field_Info = detectFields(firingMaps.forward.rateMaps{kk}{7});
        if isempty(Field_Info)
            toneField = 0;
        else 
            toneField = 1;
        end

        if length(firingMaps.forward.rateMaps{kk})==28
            Field_Info = detectFields(firingMaps.forward.rateMaps{kk}{28});
            if isempty(Field_Info)
                linEndField = 0;
            else 
                linEndField = 1;
            end
        else
            linEndField = 0;
        end    

        if cellType == 1 && linField==1
            corrtemp = corrcoef(Maps2D.firingMaps.forward.rateMaps{kk}{2}(:),Maps2D.firingMaps.forward.rateMaps{kk}{3}(:),'rows','pairwise');
            Summary.linCorr = [Summary.linCorr;corrtemp(1,2)];
        end
        
        if cellType == 1 && (linField==1 || toneField==1)
            corrtemp = corrcoef(Maps2D.firingMaps.forward.rateMaps{kk}{1}(:),Maps2D.firingMaps.forward.rateMaps{kk}{5}(:),'rows','pairwise');
            Summary.toneNoToneCorr = [Summary.toneNoToneCorr;corrtemp(1,2)];
        end
        
        if cellType == 1 && (linEndField==1 || toneField==1)
            corrtemp = corrcoef(Maps2D.firingMaps.forward.rateMaps{kk}{4}(:),Maps2D.firingMaps.forward.rateMaps{kk}{5}(:),'rows','pairwise');
            Summary.tonelinEndCorr = [Summary.tonelinEndCorr;corrtemp(1,2)];
        end
        
        if cellType == 1 && (linField==1 && linEndField==1)
            corrtemp = corrcoef(Maps2D.firingMaps.forward.rateMaps{kk}{1}(:),Maps2D.firingMaps.forward.rateMaps{kk}{4}(:),'rows','pairwise');
            Summary.linlinEndCorr = [Summary.linlinEndCorr;corrtemp(1,2)];
        end        
    end
end
end