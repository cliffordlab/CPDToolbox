function ecp = binarySegmentation(data, ctype, penalty, minThreshold)
% 
% Overview
%   Implementation of Binary Segmentation algorithm. Changepoints are the
%   last points of each segment.
%     
% Input
%   data
%   ctype:  == 1: Looking for mean changes in Gaussian data
%           == 2: Looking for variance changes in Gaussian data
%   penalty: Penalty term, can be 'BIC', 'AIC' or Hannan-Quinn ('HQ')
%   minThreshold: Limiting the number of returned significant changes by 
%   applying the additional penalty to each prospective changepoint
%
% Output
%   ecp: Estimated changepoints
%
% Reference(s)4
%   [1] Vostrikova, L. I. "Detection of the disorder in multidimensional 
%   random-processes." Doklady Akademii Nauk SSSR 259.2 (1981): 270-274.
%   [2] Chen, Jie, and Arjun K. Gupta. Parametric statistical change point 
%   analysis: with applications to genetics, medicine, and finance. Springer
%   Science & Business Media, 2011.
%
% Authors
%   Ayse Cakmak <acakmak3@gatech.edu>
% 
% Copyright (C) 2018 Authors, all rights reserved.
%
% This software may be modified and distributed under the terms
% of the BSD license. See the LICENSE file in this repo for details.
%
if (ctype == 1 || ctype == 2); paramNum = 1; end

if strcmp(penalty, 'BIC')
    beta = paramNum * log(length(data));
elseif strcmp(penalty, 'AIC')
    beta = 2 * paramNum;
elseif strcmp(penalty, 'HQ') % Hannan-Quinn
    beta = 2 * paramNum * log(log(length(data)));
end

% Make a queue to hold intervals to be segmented
queue{1} = [1 length(data)];

ecp = [];
iSeg = 1;
iQueue = 1;
while iQueue <= length(queue)
    
    % Re-initialize
    clear t;

    % Move a sliding pointer from left to right of the interval
    sigStart = queue{iQueue}(1);
    sigEnd = queue{iQueue}(2);
    iCost = 1;
    costSplit = zeros(sigEnd - sigStart + 1, 1);
    for i = sigStart:sigEnd
        
        % At each position measure total cost of each resulting segment
        if ctype == 1 % Changes in mean
            v = var(data);
            costSplit(iCost) = cost1(data(sigStart:i), v) + cost1(data((i+1):sigEnd), v);       
        elseif ctype == 2 % Changes in variance
            m = mean(data);
            costSplit(iCost) = cost2(data(sigStart:i), m) + cost2(data((i+1):sigEnd), m);       
        end
        
        iCost = iCost + 1;
    end

    % Find min of costSplit and split point 
    [costMin, iCost] = min(costSplit);
    iSig = iCost + sigStart - 1;

    % Find cost without splitting
    if ctype == 1
        costAll = cost1(data(sigStart:sigEnd), v);
    elseif ctype == 2
        costAll = cost2(data(sigStart:sigEnd), m);
    end
    
    % If cost of splitting is lower, split and put subsegments in queue
    if costAll >= costMin + beta + minThreshold
        if (iSig > sigStart+1) && (iSig < sigEnd-1)
            
            % Add changepoint
            ecp(iSeg) = iSig;
            iSeg = iSeg + 1;
            % Queue sub-segments
            queue{length(queue)+1} = [sigStart iSig];
            queue{length(queue)+1} = [iSig sigEnd];
        end
    end

    iQueue = iQueue + 1;
end

% Organize changepoints
ecp = sort(ecp);
ecp = ecp(:);

% If end of signal is not found as changepoint, add it 
if isempty(find((ecp - length(data)) == 0, 1))
    ecp = [ecp; length(data)];
end

end

% Cost function 1: Mean changes in Gaussian
function c = cost1(segData, v)
    c = 0;
    m = mean(segData);
    for i = 1:length(segData)
        c = c + (segData(i) - m).^2;
    end
    
    c = c / v;
end

% Cost function 2: Variance changes in Gaussian
function c = cost2(segData, m)
   tempSum = 0;
   n = length(segData);
   for i = 1:length(segData)
       tempSum = tempSum + (segData(i) - m).^2;
   end
   
   c = n * log(tempSum / n) + n;
end
