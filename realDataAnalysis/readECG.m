function [ECGAll, ECGTimeAll, subjectID] = readECG(pathToECG)
% 
% Overview
%  Before running this code, save all ECG files as text in pathToECG.
%   
% Input
%   pathToECG: Path to folder of ECG files
%
% Output
%   ECGAll: Cell array of ECG files
%   ECGTimeAll: Seconds ellapsed since the start of recording
%   ECGFrequency: Cell array of ECG sampling frequencies
%   subjectID: Cell array of subject IDs
%
% Authors
%   Ayse Cakmak <acakmak3@gatech.ecu>
% 
% Copyright (C) 2017 Authors, all rights reserved.
%
% This software may be modified and distributed under the terms
% of the BSD license. See the LICENSE file in this repo for details.

% Find names of ECG files in the input path
files = dir(pathToECG);
files = files(~cellfun('isempty', {files.date})); 
files(strncmp({files.name}, '.', 1)) = [];
idxECG = find(endsWith({files.name}, {'csv'}, 'ignorecase', true))';

% Read data
ECGAll = cell(1);
ECGTimeAll = cell(1);
subjectID = cell(1);
for iSubject = 1:length(idxECG)
    
    fprintf('Looking at ECG file of subject %d\n', iSubject);

    ECGFilename = files(idxECG(iSubject)).name;
    fullPathToECG = [pathToECG filesep ECGFilename];
    
    fid = fopen(fullPathToECG);
    C = textscan(fid,'%f %f','Delimiter',',','TreatAsEmpty',{'-'});    
    fclose(fid);
    ecg = C{1,2};
    ecgTime = C{1,1};
    
    % Remove nan entries
    ecg(isnan(C{1,2})) = [];
    ecgTime(isnan(C{1,2})) = [];
    
    ECGAll{iSubject,1} = ecg;
    ECGTimeAll{iSubject,1} = ecgTime;     
    subjectID{iSubject,1} = ECGFilename(1:end-4);
    
    clear C;
end

end