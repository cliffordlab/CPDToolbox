function [meanTPR, meanFNR, meanPPV, meanF1] = accuracyAnalysis(ecpAllTime, t_nnAll, t_SleepCPAll)
% 
% Overview
%   Calculates accuracy of all methods by comparing estimated changepoints
%   with sleep state change points
%     
% Input
%   ecpAllTime: Time of estimated changepoints for all methods during sleep
%   t_nnAll: Time of cleaned RR (NN) intervals
%   t_SleepCPAll: Sleep state changepoints which will be used as true
%   positives
%   pathToRepo: Path to CPD Toolbox
%
% Output
%   meanTPR: Mean True Positive Rate
%   meanFNR: Mean False Negative Rate
%   menPPV: Mean Positive Predictive Value
%
% Authors
%   Ayse Cakmak <acakmak3@gatech.edu>
% 
% Copyright (C) 2017 Authors, all rights reserved.
%
% This software may be modified and distributed under the terms
% of the BSD license. See the LICENSE file in this repo for details.

% Test with 5 different tolerance values
tolerance = [3 6 9 12 15];
resultsTPRAll = cell(length(tolerance),5);
resultsFNRAll = cell(length(tolerance),5);
resultsPPVAll = cell(length(tolerance),5);

for iTolerance = 1:length(tolerance)
    for iMethod = 1:8
        resultsTPR = zeros(length(size(ecpAllTime,1)),1);
        resultsFNR = zeros(length(size(ecpAllTime,1)),1);
        resultsPPV = zeros(length(size(ecpAllTime,1)),1);
        
        for iSubject = 1:size(ecpAllTime,1)
            ecpTime = ecpAllTime{iSubject,iMethod};
            nnTime = t_nnAll{iSubject};
            
            % Correct Changepoints that are too close together
            ecpCorrected = ecpTimeCorrector(ecpTime);

            % Sleep state change times are true changepoints
            tcp = t_SleepCPAll{iSubject};
            tcp = tcp(tcp < nnTime(end));

            results = assessCP(tcp, ecpCorrected, tolerance(iTolerance));
            resultsTPR(iSubject) = results.tpr; 
            resultsFNR(iSubject) = results.fnr; 
            resultsPPV(iSubject) = results.ppv; 
        end
        resultsTPRAll{iTolerance,iMethod} = resultsTPR;
        resultsFNRAll{iTolerance,iMethod} = resultsFNR;
        resultsPPVAll{iTolerance,iMethod} = resultsPPV;
    end
end

meanTPR = round(cellfun(@mean, resultsTPRAll),2);
meanFNR = round(cellfun(@mean, resultsFNRAll),2);
meanPPV = round(cellfun(@mean, resultsPPVAll),2);
meanF1 = 2 .* meanTPR .* meanPPV ./ (meanTPR + meanPPV);

function ecpCorrected = ecpTimeCorrector(ecp)
j = 1;
i = 1;
while i < length(ecp)
    if ecp(i+1) - ecp(i) < 2
        ecpCorrected(j) = (ecp(i) + ecp(i+1)) / 2;
        i = i + 2;
    else
        ecpCorrected(j) = ecp(i);
        i = i + 1;
    end
    j = j + 1;
end

end

end