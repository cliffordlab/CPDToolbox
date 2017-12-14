function ecp = segWin(data,lambda,winSize)
%
% OVERVIEW
%   Finds change points using original Bayesian Online Change Point Detection.
%   Uses windowing to evaluate long signals
%
% INPUT
%   data
%   lambda:     Parameter for BOCPD method
%   winSize:    Window size to be analyzed
%
% REFERENCES
%
% Copyright (C) 2017 Ayse Cakmak <acakmak3@gatech.edu>
% All rights reserved.
%   
% This software may be modified and distributed under the terms
% of the BSD license.  See the LICENSE file in this repo for details.
%
ecp = [];
winStart = 1;
iWin = 0;

while winStart < length(data)
    winEnd = winStart+winSize;
    if winEnd > length(data)
        winEnd = length(data);
    end
    
    dataWindow = data(winStart:winEnd);
    [~, maxes] = bocd(dataWindow, lambda);
    [~, ecpNew] = findpeaks(maxes(:,1));
    
    ecpNew = ecpNew + iWin * winSize; % Adjust location of ECPs 
    ecp = [ecp; ecpNew]; % Store new ECP
    winStart = winEnd;
    
    iWin = iWin+1;
    fprintf('\nWindow #%d completed\n', iWin);
end

end