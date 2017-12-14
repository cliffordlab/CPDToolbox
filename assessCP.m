function results = assessCP(tcp, ecp, tolerance)
%   
%   OVERVIEW
%   	Quantifies how well change points are detected
%
%   INPUT
%       tcp,ecp         - true change point, estimated change point
%       tolerance       - if estimated change point is within this distance
%                         of a true change point, then the estimated c.p. is
%                         considered a true positive
%       * tcp, ecp and tolerance must be in same units, for example all of
%       them in seconds or samples
%
%   OUTPUT
%       results.tpr     - true positive rate (field of structure)
%       results.fpr     - false positive rate (field of structure)
%
%   DEPENDENCIES
%
%   AUTHORS
%       Ayse Selin Cakmak <acakmak3@gatech.edu>
%
%   COPYRIGHT (C) 2016 AUTHOR(S)
%
%   LICENSE
%       This file is covered by the LICENSE file
%       in the parent directory of this GitHub repo.
%

% Initialize field values in results structure
results.tpr = 0;
results.fpr = 0;
results.fnr = 0;

% Initialize true and false positive counters
tp = 0; % hits
fp = 0; % misses
fn = 0;
idx = [];

% Loop through each estimated change point
for i = 1:length(ecp)
    
    % Isolate index of i'th estimated change point
	ii = ecp(i);
    
    % Check if i'th estimated change point is within
    % 'tolerance' of any true change points
    ii_idx_detected = find(abs(tcp - ii) <= tolerance, 1);
    
    % If our estimated change point is not among true,
    % increment false positive count, or misses
    if isempty(ii_idx_detected)
        fp = fp + 1;
        idx(i) = 0;
    % If we find a change point, the above variable is not empty;
    % increment true positive count, or hits
    else
        idx(i) =  ii_idx_detected;
        if i>1
            if idx(i)~=idx(i-1)
                tp = tp+1;
            end
        else
        tp = tp + 1;
        end
    end
end

idx = idx';

for i=1:length(tcp)
    k = find(idx==i);
    if isempty(k)
        fn = fn + 1;
    end 
end

% Count number of tcp's
results.num_tcps = length(tcp);

% Count number of ecp's
results.num_ecps = length(ecp);

if ~isempty(ecp)
    % Calculate true positive rate, or true positives / all positives
    % == sensitivity
    results.tpr = tp / (tp + fn);

    % Calculate positive predictive value
    results.ppv = tp / (tp + fp);

    % Calculate false negatives rate, or number of false negatives / all positives
    results.fnr = fn / (tp + fn);

    % Number of false positives & true postives
    results.fp = fp;
    results.tp = tp;
    
    % F1 score
    results.f1 = 2 * tp / (2 * tp + fp + fn);
    
else
    results.tpr = 0;
    results.ppv = 0;
    results.fnr = 1;
    results.fp = 0;
    results.tp = 0;
    results.f1 = 2 * tp / (2 * tp + fp + fn);
end

end