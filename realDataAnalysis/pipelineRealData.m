function [accuracy, AUC, meanTPR, meanFNR, meanPPV] = ...
                        pipelineRealData(pathToRepo, pathToECG, pathToHyp)
% 
% Overview
%   Pipeline of real data analysis in 'A comparison of changepoint detection algorithms'
%
% Input
%   pathToECG: Path to folder of ECG files ex.:'~/repos/CPDToolbox/testData/physionetCapslpdb/ecg/';
%   pathToRepo: Path to CPD Toolbox ex.: '~/repos/CPDToolbox'
%   pathToHyp:  Path to hypnogram files ex.: '~/repos/CPDToolbox/testData/physionetCapslpdb/hypnogram/'
%
% Output
%   Outputs from analysis 1 (accuracy): meanTPR, meanFNR, meanPPV
%   Outputs from analysis 2 (classification): accuracy, AUC
%   
% Dependencies
%   https://github.com/cliffordlab/HRVToolbox
%   https://www.mathworks.com/matlabcentral/fileexchange/31900-edfread
%
% Authors
%   Ayse Cakmak <acakmak3@gatech.edu>
% 
% Copyright (C) 2017 Authors, all rights reserved.
%
% This software may be modified and distributed under the terms
% of the BSD license. See the LICENSE file in this repo for details.

addpath(genpath(pathToRepo));
close all;  % Close all figures
clc;        % Clear command window
cd(pathToRepo);

%% Step 1: Read ECG data 
% pathToECG should contain all the ECG files in text format (n16 does not have ECG)
[ECGAll, ECGTimeAll, subjectID] = readECG(pathToECG);

% This matrix has parameters necessary for analysis
load([pathToRepo filesep 'realDataAnalysis' filesep 'parameters.mat']);

%% Step 2: Read hypnograms
[t_SleepCPAll, hypAll] = extractHypnogram(pathToHyp);

%% Step 3: Extract & clean RR Intervals
[nnAll, t_nnAll] = extractRR(ECGAll, ECGTimeAll, parameters);

%% Step 4: Segment RR Intervals
ecpAll = segmentRR(nnAll, subjectID, pathToRepo);

%% Step 5: Get results of bcp R code
if exist([pathToRepo filesep 'realDataAnalysis' filesep 'bcpResultsRDA' filesep '1']) ~= 0 
    ecpBCP = cell(1);
    for iFile = 1:length(subjectID)
        fileID = fopen([pathToRepo filesep 'realDataAnalysis' filesep 'bcpResultsRDA' filesep num2str(iFile)]);
        C = cell2mat(textscan(fileID,'%f'));
        ecpBCP{iFile,1} = C';
        fclose(fileID);
    end
    % Combine with other method ECPs and organize according to order in paper
    ecpAll1 = ecpAll(:,1:5);
    ecpAll2 =  ecpAll(:,6:7);
    clear ecpAll;
    ecpAll = [ecpAll1 ecpBCP ecpAll2];
    clear ecpAll1; clear ecpAll2; clear C;
else
    fprintf('Run bcpCp in RStudio please')
end

%% Step 6: Remove subjects that yielded too low number of NN intervals during sleep

% Adjust everything according to the offset
ecpTimeSleep = cell(1); t_nnAllSleep = cell(1); nnAllSleep = cell(1);
t_SleepCPAllSleep = t_SleepCPAll;
for iSubject = 1:length(nnAll)
    nnTime = t_nnAll{iSubject};
    nn = nnAll{iSubject};
    sleepLength = length(hypAll{iSubject}) * 30;
    
    offset = parameters(iSubject).timeDiff;
    nnIdx = find(nnTime >= offset & nnTime < sleepLength + offset);
    t_nnAllSleep{iSubject,1} = nnTime(nnIdx) - offset;
    nnAllSleep{iSubject,1} = nn(nnIdx);
    
    for iMethod = 1:8
       ecp = ecpAll{iSubject,iMethod};  
       ecpTime = nnTime(ecp);
       ecpTimeSleep{iSubject,iMethod} = ecpTime(ecpTime >= offset & ...
              ecpTime < (sleepLength + offset)) - offset;
    end
end

% Remove records
idxRemoveSleep = [];
for iSubject = 1:length(nnAllSleep)
    if length(nnAllSleep{iSubject}) < 10000
        idxRemoveSleep = [idxRemoveSleep; iSubject];
    end
end

fprintf('Subject %s was removed from sleep analysis. \n', subjectID{idxRemoveSleep});

nnAllSleep(idxRemoveSleep) = [];
t_nnAllSleep(idxRemoveSleep) = [];
subjectID(idxRemoveSleep) = [];
ecpTimeSleep(idxRemoveSleep,:) = [];
t_SleepCPAllSleep(idxRemoveSleep,:) = [];
ecpAll(idxRemoveSleep,:) = [];
t_nnAll(idxRemoveSleep) = [];

%% Step 7: Accuracy Analysis
[meanTPR, meanFNR, meanPPV, meanF1] = accuracyAnalysis(ecpTimeSleep, t_nnAllSleep, t_SleepCPAllSleep);

%% Step 8: Classification
[accuracy, AUC, AP, TPR_classify, PPV_classify] = ...
                                   classificationAnalysis(ecpAll, t_nnAll);

end