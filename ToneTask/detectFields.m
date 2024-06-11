function Field_Info = detectFields(SmoothedFiringRate,varargin)

    p = inputParser;
    addParameter(p,'maxRate',5,@isnumeric);
    addParameter(p,'minFieldSize',4,@isnumeric);
    addParameter(p,'maxFieldSize',35,@isnumeric);
    addParameter(p,'percentRate',0.2,@isnumeric);

    parse(p,varargin{:});
    maxRate = p.Results.maxRate;
    minFieldSize = p.Results.minFieldSize;
    maxFieldSize = p.Results.maxFieldSize;
    percentRate = p.Results.percentRate;

    % Pad on each end with zeros for edge effects. If SmoothedFiringRate
    % has nans at the end, pad around the nans. 
    % Find the index of the last non-NaN element
    last_nonNaN_index = find(~isnan(SmoothedFiringRate), 1, 'last');
    SmoothedFiringRate = SmoothedFiringRate(1:last_nonNaN_index);
    SmoothedFiringRate = [0 0 SmoothedFiringRate 0 0];
    if length(SmoothedFiringRate) < 12
        minpeakdistance = 0;
    else
        minpeakdistance = 10;
    end
    [peakValues, peakLocations] = findpeaks(SmoothedFiringRate, 'minpeakheight',maxRate, 'minpeakdistance', minpeakdistance);
    Field_Info = [];
    for j = 1:length(peakLocations)
        FieldPeak = peakLocations(j);
        % FieldPeak must be 5 Hz or more
        if peakValues(j) <= maxRate, continue, end
        LookForward = FieldPeak+1:length(SmoothedFiringRate);
        LookBack = 1:FieldPeak-1;
        PercentPeakRate_Forward = SmoothedFiringRate(LookForward)./peakValues(j);
        PercentPeakRate_Back = SmoothedFiringRate(LookBack)./peakValues(j);
        tempInd1 = find(PercentPeakRate_Forward < percentRate);
        if isempty(tempInd1), continue, end
        FieldEnd = FieldPeak+tempInd1(1); % this is the first bin forward of the animal that has a FR less than 10% of the peak rate
        tempInd2 = find(PercentPeakRate_Back < percentRate);
        if isempty(tempInd2)
            FieldStart = 1;
        else
            FieldStart = tempInd2(end); % this is the first bin forward of the animal that has a FR less than 10% of the peak rate
        end
        % Field must be more than 10 cm and less than 80cm
        if FieldEnd>FieldStart && FieldEnd-FieldStart < minFieldSize, continue, end
        if FieldEnd>FieldStart && FieldEnd-FieldStart > maxFieldSize, continue, end        
        if FieldEnd>50
            FieldEnd = 50;
        else
            FieldEnd = FieldEnd-2;
        end
        if FieldStart<3
            FieldStart = 1;
        else
            FieldStart = FieldStart-2;
        end
        Field_Info = [Field_Info; peakValues(j), FieldStart, FieldEnd, FieldPeak-2];
    end

end