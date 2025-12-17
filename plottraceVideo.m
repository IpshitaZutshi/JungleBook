function M = plottraceVideo_fullwidth(lfp)

% --- params ---
Fs = round(1/median(diff(lfp.timestamps)));
win_len_s   = 0.5;
win_len_samp = round(win_len_s*Fs);
samples_per_frame = 1;
target_fps  = 20;
line_width  = 1.5;

% --- data prep (your indexing as before) ---
[~,idx] = min(abs(lfp.timestamps-10826));
data = double(lfp.data(idx:idx+12500, :));
data = data(:, [1:7 9:64]);
nChan = size(data,2);
nSamp = size(data,1);
arr   = (0:nChan-1)*1000;
data  = 2.5*data - arr;

% colormap
% cm = cbrewer('seq','PuBuGn',90);
% cm(cm<0) = 0;

nColors = 64;  % number of levels (like in cbrewer)

navy   = [0.10 0.10 0.35];  % deep navy
teal   = [0.00 0.60 0.60];  % teal
gold   = [0.85 0.65 0.10];  % warm gold
purple = [0.45 0.20 0.50];  % muted purple

%anchor_pos   = [1, round(0.45*nColors), round(0.50*nColors), nColors];
anchor_pos   = [1, round(0.50*nColors), nColors];
anchor_colors = [navy; teal; purple];

% Interpolate smoothly across the three anchor colors
cm = interp1(anchor_pos, anchor_colors, 1:nColors, 'linear');

% --- figure/axes (full HD) ---
fig = figure('Color','w','Position',[100 100 1920 1080], ...
             'MenuBar','none','ToolBar','none','Renderer','opengl');
ax  = axes('Parent',fig,'Position',[0 0 1 1]); axis(ax,'off'); hold(ax,'on');

% resampling to 1920 px width
Xfixed = linspace(0, win_len_samp-1, 1920);
% use gridded interpolant for speed (1D, per channel)
Xsrc = (0:win_len_samp-1)';
F = cell(1,nChan);
firstWin = data(1:win_len_samp,:);
for j=1:nChan, F{j} = griddedInterpolant(Xsrc, firstWin(:,j), 'linear'); end

% initialize lines
h = gobjects(nChan,1);
Y0 = zeros(1920,nChan);
for j=1:nChan
    Y0(:,j) = F{j}(Xfixed');  % evaluate first window
    h(j) = plot(ax, Xfixed, Y0(:,j), 'LineWidth', line_width, 'Color', cm(j,:));
end
xlim(ax,[Xfixed(1) Xfixed(end)]);
ylim(ax,[min(data(:))-200, max(data(:))+200]);

% --- video writer (streaming: no frame array in memory) ---
V = VideoWriter('C:\Users\ipshi\NYU Langone Health Dropbox\Ipshita Zutshi\Website\trace_fullHD_stream','MPEG-4');
V.FrameRate = target_fps; V.Quality = 100; open(V);

start_max = nSamp - win_len_samp;
k = 1;

while k <= start_max
    idx_win = k:(k+win_len_samp-1);
    % update interpolants’ data (no realloc of lines)
    for j=1:nChan
        F{j}.Values = data(idx_win,j);      % replace window values
        set(h(j),'YData', F{j}(Xfixed'));   % evaluate & set
    end

    drawnow limitrate nocallbacks;
    frame = getframe(fig);      % capture
    writeVideo(V, frame);       % write immediately
    % release the large frame struct
    clear frame

    k = k + samples_per_frame;
end

close(V);
close(fig);
end