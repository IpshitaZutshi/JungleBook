function session = removeChannelSess(session,ch)

session.nChannels = session.nChannels-1;
session.channels = session.channels(1:(end-1));

for i=1:length(session.AnatGrps)
    d = find(session.AnatGrps(i).Channels==ch);

    if ~isempty(d)
        session.AnatGrps(i) = [];
        session.ElecGp(i) = [];
    end
end

session.nElecGps = session.nElecGps-1;

end