function plotTimeSeries(data, time, tcp, ecp, title_text)
%
%   Overview
%       Plot time series and true and estimated change points
%   
%   Input
%       rr:             1D vector of RR intervals
%       tcp:            1D vector of true change points (indices)
%       ecp:            1D vector of estimated change points (indices)
%       title_text:     string; becomes title of plot      
%   
%   Authors
%       Erik Reinertsen <er@gatech.edu>
%
%   Dependencies
%       https://github.com/erikrtn/dataviz
%   
%   Copyright (C) 2017 Erik Reinertsen <er@gatech.edu>
%   All rights reserved.
%   
%   This software may be modified and distributed under the terms
%   of the BSD license.  See the LICENSE file in this repo for details.

figure('Position', [10 10 1200 175]);

% Load google colors
google_colors = loadGoogleColors();

% Calculate time for true and estimated changepoints
tcp = time(tcp);
ecp = time(ecp);

% Plot RR intervals
plot(time, data, '.-', 'color', google_colors.blue, ...
    'markersize', 1.5, 'linewidth', 1.5);

% Plot actual change points over subplot in red
tcp_len = length(tcp);
lh = line(repmat(tcp', 2, 1), [repmat(-2,1,tcp_len); repmat(2,1,tcp_len)], ...
          'linewidth', 1.5);
set(lh, 'color', google_colors.red);

% Plot estimated change points over subplot in green
if ~isempty(ecp)
    ecp_len = length(ecp);
    lh = line(repmat(ecp', 2, 1), [repmat(-2,1,ecp_len); repmat(2,1,ecp_len)], ...
              'linewidth', 1.5);
    set(lh, 'color', google_colors.green);
end

% Settings for the figure
title(title_text);
xlabel('Time from start (sec)');

if isempty(ecp)
    legend({'Samples', 'Changepoints'});
else
    h = zeros(3, 1);
    hold on;
    h(1) = plot(NaN, NaN, '-', 'Color', google_colors.blue, 'linewidth', 2);
    h(2) = plot(NaN, NaN, '-', 'Color', google_colors.red, 'linewidth', 2);
    h(3) = plot(NaN, NaN, '-', 'Color', google_colors.green, 'linewidth', 2);
    legend(h, 'Samples', 'True CPs', 'Estimated CPs');
end

% Adjust figure limits
xlim([0 max(time)]);
ylim([min(data) max(data)]);
xticks(0:200:length(time));

% Set font size
set(gca, 'fontsize', 18);

% Set background color to white
set(gcf, 'Color', 'w');

% Expand axes to fill figure
tightfig;

end % end function