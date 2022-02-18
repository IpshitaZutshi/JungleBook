function [StimFreq, PhaseData, cell_count, Field] = getPhasePrecessionPlots2(root, cell_count, StimFreq, PhaseData, Field,plotfig)

    minFieldSize = 6;
    maxFieldSize = 40;
    % For each cell...
    for Cell = 1:size(root.b_myvar.Cells,2)

        %%ONLY get the same cells as those used in the place cell analysis
        SmoothedFiringRate = root.b_myvar.Cells(Cell).SmoothedBinnedRate;
        if nanmean(SmoothedFiringRate)== 0 || nanmean(SmoothedFiringRate)> 5
            continue;
        end
        
        SmoothedFiringRate = root.b_myvar.Cells(Cell).Lap.OFF_PooledSmoothedBinnedFiringRate;%SmoothedBinnedRate;%Lap.ON_PooledSmoothedBinnedFiringRate;
        [peakValues, peakLocations] = findpeaks(SmoothedFiringRate, 'minpeakheight',5, 'minpeakdistance', 20);
        Field_Info = [];
        for j = 1:length(peakLocations)
            FieldPeak = peakLocations(j);
            % FieldPeak must be 5 Hz or more
            if peakValues(j) < 5, continue, end
            LookForward = [FieldPeak+1:116,1:FieldPeak-1];
            LookBack = fliplr(LookForward);
            PercentPeakRate_Forward = SmoothedFiringRate(LookForward)./peakValues(j);
            PercentPeakRate_Back = SmoothedFiringRate(LookBack)./peakValues(j);
            tempInd1 = find(PercentPeakRate_Forward < .2);
            if isempty(tempInd1), continue, end
            FieldEnd = LookForward(tempInd1(1)); % this is the first bin forward of the animal that has a FR less than 20% of the peak rate
            tempInd2 = find(PercentPeakRate_Back < .2);
            FieldStart = LookBack(tempInd2(1)); % this is the first bin forward of the animal that has a FR less than 20% of the peak rate
            % Field must be at least 10 bins and less than 40 bins in the length (more than 20 cm and less than 80cm)
            if FieldEnd>FieldStart && FieldEnd-FieldStart < minFieldSize, continue, end
            if FieldEnd<FieldStart && FieldEnd+(116-FieldStart) < minFieldSize, continue, end
            if FieldEnd>FieldStart && FieldEnd-FieldStart > maxFieldSize, continue, end
            if FieldEnd<FieldStart && FieldEnd+(116-FieldStart) > maxFieldSize, continue, end
            if FieldEnd<FieldStart
                FieldEnd = 116+FieldEnd;
            end
            Field_Info = [Field_Info;FieldStart, FieldEnd, FieldPeak];
        end
        if isempty(Field_Info)
             continue;
        end
        
        minFieldSizePP = 3;      
        minSpikesPPPooled = 3;
        
        
      %%Initialize
      bSp_ON = [];
      bSp_OFF = [];
      pSp_ON = [];
      pSp_OFF = [];

      %Build array for spike phase and bins
      TT = root.b_myvar.Cells(Cell).Tetrode;
      C = root.b_myvar.Cells(Cell).Cell;
      Allbins = root.b_myvar.Cells(Cell).AllSpikeBins;
      Allphases = root.b_myvar.Cells(Cell).AllSpikePhases;
      tSp = root.spike(TT,C).ts;
      StimONOFF = root.b_myvar.StimPeriod_i_vid;
      for t = 1:length(tSp)
          [~, ind] = min(abs(root.ts - tSp(t)));
          if StimONOFF(ind) == 1 %%Stim was ON;
              bSp_ON = [bSp_ON; Allbins(t)];
              pSp_ON = [pSp_ON; Allphases(t)];
          elseif StimONOFF(ind) == 0 %%Stim was OFF;
              bSp_OFF = [bSp_OFF; Allbins(t)];
              pSp_OFF = [pSp_OFF; Allphases(t)];
          end
      end

            %%Get phase precession data
            %for ij = 1:root.b_myvar.Cells(Cell).nfields
            for ij = 1:size(Field_Info,1)
            flag = 0;
                
              Field = Field+1;

              LOGS = ismember(bSp_OFF, Field_Info(ij,1):1:Field_Info(ij,2));
              AllPosition = bSp_OFF(LOGS);
              AllPhase = pSp_OFF(LOGS);
              PositionField = Field_Info(ij,1):1:Field_Info(ij,2);
                    
             
              %if length(AllPosition)>minSpikesPPPooled && range(AllPosition)>minFieldSizePP,
              %if length(AllPosition)>minSpikesPPPooled && range(PositionField)>minFieldSizePP
                    [aopt_pooledB phi0_pooledB rho_pooled R_pooled p_pooledB lindatB phidatB] = Thetapacing.PhasePrecessionStats(AllPosition, AllPhase, PositionField);  
               %     else aopt_pooledB = NaN; phi0_pooled = NaN; rho_pooled = NaN; R_pooled = NaN; p_pooledB = NaN; lindatB = NaN; phidatB = NaN;       
              %end    
              
              if aopt_pooledB<0 && p_pooledB<0.05 && plotfig==1
                  figure
                  subplot(2,2,1)
                  Functions.plotPhasePrecession(lindatB,phidatB,aopt_pooledB,phi0_pooledB,1)
                  %plot([lindatB; lindatB], [rad2deg(phidatB); rad2deg(phidatB)+360],'k.')
                  %xlim([0 1])
                  title(strcat('Slope: ',num2str(aopt_pooledB),' pVal: ',num2str(p_pooledB)));
                  subplot(2,2,2)
                  imagesc(root.b_myvar.Cells(Cell).Lap.OFF_PooledSmoothedBinnedFiringRate')
                  %imagesc(root.b_myvar.Cells(Cell).SmoothedBinnedRate);
                  colorbar
                  title(strcat('Slope: ',num2str(aopt_pooledB),' pVal: ',num2str(p_pooledB),' FieldS-E ', num2str(Field_Info(ij,1)),' ',num2str(Field_Info(ij,2))));                  
                  %title(strcat('cellnum', num2str(root.b_myvar.Cells(Cell).Tetrode),' ',num2str(root.b_myvar.Cells(Cell).Cell)));
                  flag = 1;
               
              end

              PhaseData.nSpikes.Baseline{Field} = numel(AllPosition);
              PhaseData.PooledSlopes.Baseline{Field} = aopt_pooledB; %pooled slope
              PhaseData.PooledpVals.Baseline{Field} = p_pooledB; %pooled p
              PhaseData.PooledCors.Baseline{Field} = rho_pooled; %pooled corr
              PhaseData.Pooledphi.Baseline{Field} = phi0_pooledB; %pooled y intercept

              LOGS = ismember(bSp_ON, Field_Info(ij,1):1:Field_Info(ij,2));
              AllPosition = bSp_ON(LOGS);
              AllPhase = pSp_ON(LOGS);
              PositionField = Field_Info(ij,1):1:Field_Info(ij,2);
              
              %if length(AllPosition)>minSpikesPPPooled && range(AllPosition)>minFieldSizePP,
             % if length(AllPosition)>minSpikesPPPooled && range(PositionField)>minFieldSizePP
                    [aopt_pooled phi0_pooled rho_pooled R_pooled p_pooled lindat phidat] = Thetapacing.PhasePrecessionStats(AllPosition, AllPhase, PositionField);  
                    %else aopt_pooled = NaN; phi0_pooled = NaN; rho_pooled = NaN; R_pooled = NaN; p_pooled = NaN; lindat = NaN; phidat = NaN;    
              %end    
              freq = mode(root.b_myvar.StimFreq_lfp(root.b_myvar.StimFreq_lfp>0));
                            
              if ((aopt_pooled<0 && p_pooled<0.05) || flag ==1) && plotfig ==1
                  if flag==0
                        figure
                        subplot(2,2,1)
                        Functions.plotPhasePrecession(lindatB,phidatB,aopt_pooledB,phi0_pooledB,1)
                        %plot([lindatB; lindatB], [rad2deg(phidatB); rad2deg(phidatB)+360],'k.')
                        xlim([0 1])
                        title(strcat('Slope: ',num2str(aopt_pooledB),' pVal: ',num2str(p_pooledB)));
                        subplot(2,2,2)
                        imagesc(root.b_myvar.Cells(Cell).SmoothedBinnedRate);
                        colorbar
                        title(strcat('Freq', num2str(freq), 'Slope: ',num2str(aopt_pooledB),' pVal: ',num2str(p_pooledB),' FieldS-E ', num2str(Field_Info(ij,1)),' ',num2str(Field_Info(ij,2))));             
                  end
                                          
                  subplot(2,2,3)
                  Functions.plotPhasePrecession(lindat,phidat,aopt_pooled,phi0_pooled,1)
                  %plot([lindat; lindat], [rad2deg(phidat); rad2deg(phidat)+360],'k.')
                  xlim([0 1])
                  title(strcat('Freq', num2str(freq),'Slope: ',num2str(aopt_pooled),' pVal: ',num2str(p_pooled),' FieldS-E ', num2str(Field_Info(ij,1)),' ',num2str(Field_Info(ij,2))));
                  subplot(2,2,4)
                  imagesc(root.b_myvar.Cells(Cell).Lap.ON_PooledSmoothedBinnedFiringRate')
                  colorbar
                  title(root.name);

                    set(gcf,'renderer','painters')
      
                    saveas(gcf,['C:\Users\Ipshita\Documents\PhD documents\Theta project compiled data\Revision analysis\Revision plots\ExamplesRev\',num2str(Field),'.jpg'],'jpg');
                    saveas(gcf,['C:\Users\Ipshita\Documents\PhD documents\Theta project compiled data\Revision analysis\Revision plots\ExamplesRev\',num2str(Field),'.eps'],'epsc');
                    close all
              end

              PhaseData.nSpikes.Stim{Field} = numel(AllPosition);
              PhaseData.PooledSlopes.Stim{Field} = aopt_pooled; %pooled slope
              PhaseData.PooledpVals.Stim{Field} = p_pooled; %pooled p
              PhaseData.PooledCors.Stim{Field} = rho_pooled; %pooled corr
              PhaseData.Pooledphi.Stim{Field} = phi0_pooled; %pooled y intercept
              
              StimFreq(Field) = mode(root.b_myvar.StimFreq_lfp(root.b_myvar.StimFreq_lfp>0));

            end
       
        %%Increment cell count
        cell_count = cell_count+1;

     end
    
end
    










