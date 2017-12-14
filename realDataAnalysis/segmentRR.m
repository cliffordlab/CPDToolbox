function ecpAll = segmentRR(nnAll, subjectID, pathToRepo)
% 
% Overview
%   Segments cleaned RR time series using 7 methods. Writes RR time series
%   into text files so that these files can be used with bcp method in
%   RStudio.
%     
% Input
%   nnAll: Cleaned RR time series/NN time series
%   subjectID: ID of all patients
%
% Output
%   ecpAll: Estmated changepoints of all methods
%
% Authors
%   Ayse Cakmak <acakmak3@gatech.edu>
% 
% Copyright (C) 2017 Authors, all rights reserved.
%
% This software may be modified and distributed under the terms
% of the BSD license. See the LICENSE file in this repo for details.

ecpAll = cell(length(subjectID),7);

%% RMDM
disp('RMDM');
for iSubject = 1:length(subjectID)
    disp(subjectID{iSubject});
    NN = nnAll{iSubject};
    NN(NN < 0) = []; % Interpolation leads to negative NN sometimes
    
    P_0 = 0.95;
    l0 = 7;
    ecpAll{iSubject,1} = (rmdm(NN, P_0, l0))';
end

%% BiS
disp('BiS');
for iSubject = 1:length(subjectID)
    disp(subjectID{iSubject});
    NN = nnAll{iSubject};
    NN(NN < 0) = []; % Interpolation leads to negative NN sometimes

    ctype = int64(1); 
    n = int64(length(NN));
    p = 1;
    beta = 2 * p * log(log(length(NN)));
    minss = int64(2);

    [ecp, ~, ~] = nag_tsa_cp_binary(ctype, zscore(NN), 'n', n, 'beta', beta, ...
        'minss', minss);
    
    ecpAll{iSubject,2} = double(ecp);
end

%% PELT1
disp('PELT1');
for iSubject = 1:length(subjectID)
    disp(subjectID{iSubject});
    NN = nnAll{iSubject};
    NN(NN < 0) = []; % Interpolation leads to negative NN sometimes

    ctype = int64(1);
    n = int64(length(NN));
    p = 1;
    beta = p*log(length(NN)); 
    minss = int64(2);
    
    [ecp, ~, ~] = nag_tsa_cp_pelt(ctype, zscore(NN), 'n', n, ...
                           'beta', beta, 'minss', minss);
    ecpAll{iSubject,3} = double(ecp);
end

%% findchangepts.m / PELT2
disp('PELT2');  
for iSubject = 1:length(subjectID)
    disp(subjectID{iSubject});
    NN = nnAll{iSubject};
    NN(NN < 0) = []; % Interpolation leads to negative NN sometimes

    ecpAll{iSubject,4} = findchangepts(NN, 'Statistic', 'rms', 'MinThreshold', 0.1);
end

% Find change points with other segmentation methods
%% BBlocks
disp('BBlocks');    
for iSubject = 1:length(subjectID)
    disp(subjectID{iSubject});
    NN = nnAll{iSubject};
    NN(NN < 0) = []; % Interpolation leads to negative NN sometimes
    
    ncp_prior = 4;

    ecpAll{iSubject,5} = bblocks(NN', ncp_prior);
end

%% BCP method is on RStudio
fprintf('Starting to bcp implementation in R for real data analysis \n');

% Make folder to store results
folder = [pathToRepo filesep 'realDataAnalysis' filesep 'bcpResultsRDA'];
if exist(folder, 'dir') == 0; mkdir(folder); end
        
% Write data on text files so we can go to RStudio and use them
% Make directory to store these files
disp('BCP');
if exist([pathToRepo filesep 'realDataAnalysis' filesep 'RRIntervals'], 'dir') == 0 
    mkdir([pathToRepo filesep 'realDataAnalysis' filesep 'RRIntervals']); end
for iSubject = 1:length(subjectID)
    fileID = fopen([pathToRepo filesep 'realDataAnalysis' filesep 'RRIntervals' filesep , ...
        subjectID{iSubject},'.txt'],'w');
    fprintf(fileID,'%f\n',nnAll{iSubject});
    fclose(fileID);
end

system(['Rscript realDataAnalysis' filesep 'bcpCP.R'])

%% BOCD Original
disp('BOCD');    
for iSubject = 1:length(subjectID)
    disp(subjectID{iSubject});
    NN = nnAll{iSubject};
    NN(NN < 0) = []; % Interpolation leads to negative NN sometimes
    
    lambda = 1840;
    winSize = 30000;
    
    ecpAll{iSubject,6} = segWin(NN,lambda,winSize);
end

%% BOCD Modified
disp('mBOCD');    
for iSubject = 1:length(subjectID)
    disp(subjectID{iSubject});
    NN = nnAll{iSubject};

    NN(NN < 0) = []; % Interpolation leads to negative NN sometimes

    lambda = 80;
    ecpAll{iSubject,7} = msegWin(NN, lambda, 30000);
end

end