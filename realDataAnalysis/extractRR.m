function [nnAll, t_nnAll] = extractRR(ECGAll, ECGTimeAll, parameters)
% 
% Overview
%   Extracts and cleans RR interval time series from ECG records
%     
% Input
%   ECGAll: ECG records for all subjects
%   ECGTimeAll: Time of ECG samples (sec)
%   parameters: Matrix that contains information about getting peaks of ECG
%   signal to 1mV and sampling frequency
%
% Output
%   nnAll: Cleaned RR intervals (NN intervals)
%   t_nnAll: Time of cleaned RR intervals
%
% Dependencies
%   https://github.com/cliffordlab/HRVToolbox
%
% Authors
%   Ayse Cakmak <acakmak3@gatech.edu>
% 
% Copyright (C) 2017 Authors, all rights reserved.
%
% This software may be modified and distributed under the terms
% of the BSD license. See the LICENSE file in this repo for details.

nnAll = cell(1);
t_nnAll = cell(1);

% Put 'CPDAnalysis' as project name in InitializeHRVparams.m
s = InitializeHRVparams('CPDAnalysis');
s.debug = 0;
s.gen_figs = 0;
for iSubject = 1:length(ECGAll)
    
    fprintf('Looking at ECG file of subject %d\n', iSubject);
    
    s.Fs = parameters(iSubject).fs;

    ECG = ECGAll{iSubject};
    ECGTime = ECGTimeAll{iSubject};
    
    ECG = ECG - mean(ECG);
    
    if parameters(iSubject).unit_uV == 1
        ECG = ECG / 1000;
    end
    
    % Make peaks 1mV
    ECG = ECG * parameters(iSubject).scaling;
    
    % QRS Dection 1
    disp('QRS detection 1');
    qrs = run_qrsdet_by_seg(ECG, s); 
    
    % QRS Detection 2
    disp('QRS detection 2');
    wqrs = wqrsm(ECG, s.Fs);
    
    % QRS SQI
    disp('Looking at quality of extracted RR intervals..');
    [sqijw, StartIdxSQIwindows] = bsqi(qrs, wqrs, s);
    
    rr = diff(qrs./s.Fs);
    t_rr = ECGTime(qrs(1:end-1));
    
    % Start RR interval preprocessing
    disp('RR interval preprocessing..');
    [nn, t_nn, ~] = RRIntervalPreprocess(rr, t_rr, [], s);

    nnAll{iSubject,1} = nn;
    t_nnAll{iSubject,1} = t_nn;
end

end       
