function [out] = convert_rr_to_act(rr, snr)

%   [out] = convert2actgen(rr)
%   
%   OVERVIEW:
%   	1. Converts RR time series data into activity data
%
%   INPUT
%       rr      - Input data timeseries
%
%   OUTPUT
%       act     - Activity data timeseries
%
%   DEPENDENCIES
%       https://github.com/cliffordlab/rrgen
%
%   AUTHORS
%       Ayse Selin Cakmak <acakmak3@gatech.edu>
%
%   COPYRIGHT (C) 2016 AUTHOR(S)
%
%   LICENSE
%       This file is covered by the LICENSE file
%       in the parent directory of this GitHub repo.

act = zeros(length(rr),1);

% Step 1: Mirroring rr interval data with respect to rr's mean(making its
% max new min)
m = mean(rr);
for i=1:length(rr)
    dist = abs(m - rr(i));
    if rr(i) > m
        act(i) = m - dist;
    elseif rr(i) < m
        act(i) = m + dist;
    else 
        act(i) = rr(i);
    end
end

% Step 2: Subtract 2 (Not sure)
act = act-2;

% Step 3: Inverse log
act = exp(act);

%Step 4: Add Laplacian noise
out = SNR_Set(act', snr)';

end

