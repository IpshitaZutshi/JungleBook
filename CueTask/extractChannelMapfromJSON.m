%% Code to read the output .json files from the ElectrodeLocationsImportSaveOpenEphys.py 
% code using code from https://github.com/SpikeInterface/probeinterface/ to import
% Channel map information into Matlab

clf
clear

for ns=1
    fname = 'Z:\Homes\zutshi01\Recordings\Cue task\T17\xml files\CA1_cortex.json';
    fid = fopen(fname);
    raw = fread(fid,inf);
    str = char(raw');
    fclose(fid);
    val = jsondecode(str);

    xcoords1=val.probes.contact_positions(:,1);
    ycoords1=val.probes.contact_positions(:,2);

  
        xcoords=[xcoords1];
        ycoords=[ycoords1];

    plot(xcoords,...
        ycoords,'.')
    hold on
    text(xcoords,...
        ycoords,...
        num2str((1:length(xcoords))'))
    xcoords=xcoords';
    ycoords=ycoords';
    save(['E:\DataForNatProtocol\',FileNA{ns},'\',FileNA{ns},'_ChannelMap.mat'],'xcoords','ycoords')
end
