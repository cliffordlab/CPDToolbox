function [accuracy, AUC, AP, TPR, PPV] = classificationAnalysis(ecpAll, t_nnAll)
% 
% Overview
%   Classificarion of health controls and RBD patients using features 
%   derived from distribution of time lengths between changepoints
%     
% Input
%   ecpAll: Estimated changepoints for all changeoint detection methods
%
% Output
%   accuracy: Percent classifcation accuracy
%   AUC: Derived AUC for k-NN
%   TPR: True positive rate
%   PPV: Positive Predictive Value
%
% Authors
%   Ayse Cakmak <acakmak3@gatech.edu>
% 
% Copyright (C) 2017 Authors, all rights reserved.
%
% This software may be modified and distributed under the terms
% of the BSD license. See the LICENSE file in this repo for details.

numSubjects = size(ecpAll,1);

% Find segment lengths
segAll = cell(numSubjects,8);
for i = 1:8
    for iSubject = 1:numSubjects
        ecpSubject = ecpAll{iSubject,i};

        tnnSubject = t_nnAll{iSubject};

        seg = zeros(length(ecpSubject)-1,1);
        for iEcp = 1:length(ecpSubject)-1
            seg(iEcp) = tnnSubject(ecpSubject(iEcp + 1)) - tnnSubject(ecpSubject(iEcp)); % Found segment length in seconds
        end
        
        segAll{iSubject,i} = seg;
    end
end

% Fit Pareto distirbution to all & get the parameters
kAll = zeros(numSubjects,8); sigmaAll = zeros(numSubjects,8);
for i = 1:8
    for iSubject = 1:numSubjects
        seg = segAll{iSubject,i};
        seg(seg == 0) = [];
        parmhat = gpfit(seg);
        kAll(iSubject,i) = parmhat(1);
        sigmaAll(iSubject,i) = parmhat(2);
    end    
end

for iMethod = 1:8
    % Linear regression to seperate two groups (control/rem behaviour disorder)
    k = kAll(:,iMethod);
    sigma = sigmaAll(:,iMethod);

    % First plot data to see if k-NN is a good idea
    control_k = k(1:15); % First 15 is from control class
    control_sigma = sigma(1:15);
    X1 = [control_k, control_sigma];

    rbd_k = k(16:end);
    rbd_sigma = sigma(16:end);
    X2 = [rbd_k, rbd_sigma];

    figure(iMethod)
    scatter(control_k,control_sigma,'filled','go');
    hold on;
    scatter(rbd_k,rbd_sigma,'filled','ro');

    X = [X1; X2];
    Y = cell(numSubjects,1);
    Y(1:15) = {'control'};
    Y(16:end) = {'rdb'};
    
    % Optimize k-NN
    rng(1)
    Mdl = fitcknn(X,Y,'OptimizeHyperparameters','auto',...
        'HyperparameterOptimizationOptions',...
        struct('AcquisitionFunctionName','expected-improvement-plus'));
    optimalNeighbor{iMethod,1} = Mdl.NumNeighbors;
    optimalDistance{iMethod,1} = Mdl.Distance;

    % Cross validate model & evaluate
    kNNMdl = fitcknn(X,Y,'Distance',optimalDistance{iMethod,1}, ...
        'NumNeighbors',optimalNeighbor{iMethod,1});
    rng(1); % For reproducibility
    CVKNNMdl = crossval(kNNMdl);
    classError(iMethod,1) = kfoldLoss(CVKNNMdl);
    predictedLabel = kfoldPredict(CVKNNMdl);
    trueLabel = CVKNNMdl.Y;
    accuracy(iMethod,1) = sum(strcmp(predictedLabel, trueLabel))/length(trueLabel);  
    
    TP = sum(strcmp(predictedLabel(1:15),trueLabel(1:15)));
    FP = 15 - TP;
    TN = sum(strcmp(predictedLabel(16:end),trueLabel(16:end)));
    FN = 22 - TN;
    
    % Calculate metrics
    TPR(iMethod,1) = TP / (TP + FN);
    PPV(iMethod,1) = TP / (TP + FP);
    
    % Find distance of each point to all other points
    distanceMetric = optimalDistance{iMethod,1};
    neighborNum = optimalNeighbor{iMethod,1};
    distances = pdist2(X, X, distanceMetric);
    
    % Find score for each point = # of points from positive class in the nearest neighbors / k
    score = zeros(size(X,1),1);
    for iPoint = 1:size(X,1)
        distPoint = distances(:,iPoint);
        
        [~,closest] = sort(distPoint);
        closest = closest(2:neighborNum + 1);
        
        labelsOfClosest = trueLabel(closest);
        labelsOfClosest = strcmp(labelsOfClosest,'control');
        score(iPoint) = sum(labelsOfClosest)/neighborNum;
    end
    
    % Use derived scores to get AUC
    [~, ~, ~, AUC(iMethod,1)] = perfcurve(trueLabel,score,'control');
    [~, ~, ~, ~, AP(iMethod,1)] = prc_stats(15/37, mean(score(16:end)), ...
                  std(score(16:end)), mean(score(1:15)), std(score(1:15)));
end

TPR = round(TPR,2);
PPV = round(PPV,2);
accuracy = round(accuracy,2);
AUC = round(AUC,2);
AP = round(AP,2);

close all;

end


