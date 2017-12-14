function [t_SleepCPAll, hypAll] = extractHypnogram(pathToHyp)
% 
% Overview
%  Download all hypnogram text files from following link:
%  https://physionet.org/pn6/capslpdb/
%

% Load hypnogram files 
files = dir(pathToHyp);
files = files(~cellfun('isempty', {files.date})); 
files(strncmp({files.name}, '.', 1)) = [];
idxHyp = find(endsWith({files.name}, {'txt'}, 'ignorecase', true))';

hypAll =  cell(length(idxHyp),1);
hypFilename =  cell(length(idxHyp),1);
for iSubject = 1:length(idxHyp)
    hypFilename{iSubject} = files(idxHyp(iSubject)).name;
    fullPathToHyp = [pathToHyp, hypFilename{iSubject}];
    hyp = ScoringReader(fullPathToHyp);
    hyp = hyp(:,1);

    % Combine NREM sleep stages
    hypMerged = hyp;
    hypMerged(hypMerged == 1) = 7; 
    hypMerged(hypMerged == 2) = 7; 
    hypMerged(hypMerged == 3) = 7; 
    hypMerged(hypMerged == 4) = 7; 

    hypAll{iSubject} = hypMerged;
end

% Extract sleep state changepoints
t_SleepCPAll = cell(length(hypAll),1);
sleepCPAll = cell(length(hypAll),1);
for iSubject = 1:length(hypAll)    
    hyp = hypAll{iSubject,1};
    
    j = 1;
    cp = [];
    for i = 1:length(hyp)-1
        if hyp(i)~=hyp(i+1)
            cp(j) = i;
            j = j+1;
        end
    end
    sleepCPAll{iSubject} = cp(:);
    t_SleepCPAll{iSubject} = 30 * sleepCPAll{iSubject}; % Each epoch is 30 sec long
end

end