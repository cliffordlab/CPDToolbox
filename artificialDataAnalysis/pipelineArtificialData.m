%function pipelineArtificialData(pathToRepo)
%   
%   Overview:
%   	Loads artificially generated data and applies change point
%   	detection methods
%
%   Input:
%       pathToRepo: Path of CPDToolbox
%       ex : '~/repos/CPDToolbox'
%
%   Dependencies
%       https://github.com/cliffordlab/rrgen
%       https://github.com/erikrtn/dataviz
%
%   Authors
%       Erik Reinertsen <er@gatech.edu>
%       Ayse Selin Cakmak <acakmak3@gatech.edu>
%
%   Copyright (C) 2017 Authors
%   All rights reserved.
%   
%   This software may be modified and distributed under the terms
%   of the BSD license.  See the LICENSE file in this repo for details.

addpath(genpath(pathToRepo));
close all;  % Close all figures
clc;        % Clear command window
cd(pathToRepo);

% Settings
settings.tolerance = 5;     % tolerance in seconds for determining if change point is detected;
                            % if estimated change point (ecp) is within true change point (tcp),
                            % the change point is detected.
settings.act = 1;           % 1: work with actigraphy data 
settings.figs = 0;          % 1: plot tcp & ecp on time series
settings.iterations = 1000;

% Choose parameter set according to each method
if settings.act == 0 % Work on RR interval time series
    lambda = 80;        % Modified BOCD
    lambda2 = 1840;     % Original BOCD
    ncpPrior = 4;       % Bayesian Blocks
else
    lambda = 170;       % Modified BOCD
    lambda2 = 4490;     % Original BOCD
    ncpPrior = 3.5;     % Bayesian Blocks    
end

fprintf('Starting change point detection pipeline...\n\n');

%% Load data and apply changepoint detection methods
for iterCounter = 1:settings.iterations

    % Initialize constants for rrgen
    
    % Create random seed
    seed = rand(1)*1e3;
    
    % Generate artificial data
    if settings.act ~= 1
        % 24 hour RR time series
        vector_length = 86400;

        % Create RR time series with 'rrgen_sys' exe (compiled from c file)
        [data, tcp, sleepStart, sleepEnd] = rrgenV3_wrapper(seed, vector_length, 0, 0, pathToRepo,0);
        
        % Calculate sleep length
        sleepLength(iterCounter) = (sleepEnd - sleepStart) / 3600;

        % Extract sleep portion of data 
        data = data(sleepStart:sleepEnd);
                
        % Adjust changepoints for sleep portion
        tcp_idx = find((tcp > sleepStart) & (tcp < sleepEnd));
        tcp = tcp(tcp_idx) - sleepStart;
        
        time = cumsum(data); % Time of RR time series is found by cumulative summation
    else
        vector_length = 24000;
        [data, tcp, sleepStart, sleepEnd] = rrgenV3_wrapper(seed, vector_length, 0, 0, pathToRepo,1);
        time = 0:1:24000; 
    end
    
    % Converting rrgen data to actigraphy data
    if settings.act == 1
        % Make folders to store data if it does not exist already
        folder = [pathToRepo filesep 'artificialDataAnalysis' filesep 'generatedData' filesep 'act'];
        if exist(folder, 'dir') == 0; mkdir(folder); end
        
        folder = [pathToRepo filesep 'artificialDataAnalysis' filesep 'generatedData' filesep 'act_tcp'];
        if exist(folder, 'dir') == 0; mkdir(folder); end

        snr = 30; % choose SNR level
        data = convert_rr_to_act(data, snr);
        fprintf('Test %1.0f/%1.0f: generating actigraphy samples (length %1.0f)...\n\n', ...
        iterCounter, settings.iterations, length(data));
    
        % Write files so they can be tested with BCP R code later
        filename = strcat([pathToRepo filesep 'artificialDataAnalysis' filesep 'generatedData' filesep 'act' filesep 'act'],num2str(iterCounter),'.txt');
        fileID = fopen(filename,'w');
        fprintf(fileID,'%f \n',data);
        fclose(fileID);

        filename = strcat([pathToRepo filesep 'artificialDataAnalysis' filesep 'generatedData' filesep 'act_tcp' filesep 'tcp'],num2str(iterCounter),'.txt');
        fileID = fopen(filename,'w');
        fprintf(fileID,'%d \n',tcp);
        fclose(fileID);
        
    else
        % Make folders to store data if it does not exist already
        folder = [pathToRepo filesep 'artificialDataAnalysis' filesep 'generatedData' filesep 'rr'];
        if exist(folder, 'dir') == 0; mkdir(folder); end
        
        folder = [pathToRepo filesep 'artificialDataAnalysis' filesep 'generatedData' filesep 'rr_tcp'];
        if exist(folder, 'dir') == 0; mkdir(folder); end
        fprintf('Test %1.0f/%1.0f: generating RR samples (length %1.0f)...\n\n', ...
        iterCounter, settings.iterations, length(data));
        
        % Write files so they can be tested with BCP R code later
        filename = strcat([pathToRepo filesep 'artificialDataAnalysis' filesep 'generatedData' filesep 'rr' filesep 'rr'],num2str(iterCounter),'.txt');
        fileID = fopen(filename,'w');
        fprintf(fileID,'%f \n',data);
        fclose(fileID);

        filename = strcat([pathToRepo filesep 'artificialDataAnalysis' filesep 'generatedData' filesep 'rr_tcp' filesep 'tcp'],num2str(iterCounter),'.txt');
        fileID = fopen(filename,'w');
        fprintf(fileID,'%d \n',tcp);
        fclose(fileID);
    end
        
    % Plot data before analysis
    if settings.figs == 1
        plotTimeSeries(data, time, tcp, [], 'Time Series');
        pause(1);
        close;
    end
    
    %% Modified Bayesian Online Changepoint Detection
        
    ecp = msegWin(data, lambda, 30000);
    
    % Plot true and estimated change points
    if settings.figs == 1    
        plotTimeSeries(data, time, tcp, ecp, 'Modified BOCD');
        pause(1);
        close;
    end
    
    % Assess performance of change point detection
    % and save results to an array of structs
    results_bocd(iterCounter) = assessCP(time(tcp), time(ecp), settings.tolerance);
    
    %% BOCD Original Code
    
    ecp = segWin(data,lambda2, 30000);
    
    % Plot true and estimated change points
    if settings.figs == 1    
        plotTimeSeries(data, time, tcp, ecp, 'Original BOCD');
        pause(1);
        close;
    end
    
    results_bocdo(iterCounter) = assessCP(time(tcp), time(ecp), settings.tolerance);
    
    %% Bayesian Blocks
    
    ecp = bblocks(data, ncpPrior);
     
    % Plot true and estimated change points
    if settings.figs == 1
        plotTimeSeries(data, time, tcp, ecp, 'Bayesian Blocks');
        pause(1);
        close;
    end
    
    % Assess performance of change point detection
    results_bblocks(iterCounter) = assessCP(time(tcp), time(ecp), settings.tolerance);
    
    %% Pruned Exact Linear Time / PELT 1 
    ctype = int64(1); 
    n = int64(length(data));
    p = 1;
    beta = p*log(length(data));
    minss = int64(2);
    
    [ecp, ~, ~] = nag_tsa_cp_pelt(ctype, zscore(data), 'n', n, ...
                           'beta', beta, 'minss', minss);
    ecp = double(ecp);
                       
    % Plot true and estimated change points
    if settings.figs == 1    
        plotTimeSeries(data, time, tcp, ecp, 'PELT 1');
        pause(1);
        close;
    end
    
    % Assess performance of change point detection
    results_pelt1(iterCounter) = assessCP(time(tcp), time(ecp), settings.tolerance);
  
    %% Matlab function: findchangepts.m / PELT2

    ecp = findchangepts(data, 'Statistic', 'rms', 'MinThreshold', 0.1);
    
    % Plot true and estimated change points
    if settings.figs == 1    
        plotTimeSeries(data, time, tcp, ecp, 'PELT 2');
        pause(1);
        close;
    end
    
    % Assess performance of change point detection
    results_pelt2(iterCounter) = assessCP(time(tcp), time(ecp), settings.tolerance);
          
    %% Binary Segmentation
        
    ctype = int64(1); 
    n = int64(length(data));
    p = 1;
    beta = 2 * p * log(log(length(data)));
    minss = int64(2);

    [ecp, ~, ~] = nag_tsa_cp_binary(ctype, zscore(data), 'n', n, 'beta', beta, ...
        'minss', minss);
    
    ecp = double(ecp);

    % Plot true and estimated change points
    if settings.figs == 1    
        plotTimeSeries(data, time, tcp, ecp, 'Binary Segmentation');
        pause(1);
        close;
    end
    
    % Assess performance of change point detection
    results_bis(iterCounter) = assessCP(time(tcp), time(ecp), settings.tolerance);   
    
    %% Recursive Mean Difference Maximization
    
    P_0 = 0.95;
    l0 = 7; 
    ecp = (rmdm(data, P_0, l0))';

    % Plot true and estimated change points
    if settings.figs == 1    
        plotTimeSeries(data, time, tcp, ecp, 'RMDM');
        pause(1);
        close;
    end
    
    % Assess performance of change point detection   
    results_rmdm(iterCounter) = assessCP(time(tcp), time(ecp), settings.tolerance);    

    %% Close all figures
    close all;
    
end 
%%
fprintf('Starting to bcp implementation in R \n');

% Make folder to store results
folder = [pathToRepo filesep 'artificialDataAnalysis' filesep 'bcpResultsADA'];
if exist(folder, 'dir') == 0; mkdir(folder); end

if settings.act == 0 % Work on RR interval time series
    system(['Rscript CPDMethods' filesep 'bcp' filesep 'bcp_testRRGEN.R'])
else
    system(['Rscript CPDMethods' filesep 'bcp' filesep 'bcp_testACTGEN.R'])
end

% Load bcp results ~ Generate these results in R using codes in BCP folder
    if settings.act == 0 % Working on RR interval time series
        fileID = fopen([pathToRepo filesep 'artificialDataAnalysis' filesep 'bcpResultsADA' filesep 'rrgenTestResults']);
        C = cell2mat(textscan(fileID,'%f %f %f %f'));
        fclose(fileID);   
    else % Working on actigraphy time series
        fileID = fopen([pathToRepo filesep 'artificialDataAnalysis' filesep 'bcpResultsADA' filesep 'actgenTestResults']);
        C = cell2mat(textscan(fileID,'%f %f %f %f'));
        fclose(fileID);          
    end

results_bcp.tpr = C(:,1);
results_bcp.fnr = C(:,2);
results_bcp.ppv = C(:,3);
results_bcp.fp = C(:,4);

%% Plot change point detection performance across algorithms
figure('position', [100 100 1200 800]);

matrix_tpr = [[results_rmdm.tpr]', ...
              [results_bis.tpr]', ...
              [results_pelt1.tpr]', ...              
              [results_pelt2.tpr]', ...
              [results_bblocks.tpr]', ...
              results_bcp.tpr, ...
              [results_bocdo.tpr]', ...
              [results_bocd.tpr]'];

subplot(2,2,1);
bp1 = bplot(matrix_tpr, 'outliers');
xlabel('\bf (a)');
xticks(1:8);
xlim([0 9])
xticklabels({'RMDM', 'BiS', 'PELT1', 'PELT2', 'BBlocks', 'BCP', 'BOCD', 'mBOCD'});
xtickangle(30)
ylabel('TPR');
ylim([0 1]);
legend(bp1,'location','best');
legend('boxoff');
set(gca, 'fontsize', 16);

matrix_fnr = [[results_rmdm.fnr]', ...
              [results_bis.fnr]', ...
              [results_pelt1.fnr]', ...              
              [results_pelt2.fnr]', ...
              [results_bblocks.fnr]', ...
              results_bcp.fnr, ...              
              [results_bocdo.fnr]', ...
              [results_bocd.fnr]'];
          
subplot(2,2,2);
bp2 = bplot(matrix_fnr, 'outliers');
xlabel('\bf (b)');
xticks(1:8);
xticklabels({'RMDM', 'BiS', 'PELT1', 'PELT2', 'BBlocks', 'BCP', 'BOCD', 'mBOCD'});
xtickangle(30)
xlim([0 9])
ylabel('FNR');
ylim([0 1]);
legend(bp2,'location','best');
legend('boxoff');
set(gca, 'fontsize', 16);

matrix_ppv = [[results_rmdm.ppv]', ...
              [results_bis.ppv]', ...
              [results_pelt1.ppv]', ...              
              [results_pelt2.ppv]', ...
              [results_bblocks.ppv]', ...
              results_bcp.ppv, ...              
              [results_bocdo.ppv]', ...
              [results_bocd.ppv]'];
          
subplot(2,2,3);
bp3 = bplot(matrix_ppv, 'outliers');
xlabel('\bf (c)');
xticks(1:8);
xticklabels({'RMDM', 'BiS', 'PELT1', 'PELT2', 'BBlocks', 'BCP', 'BOCD', 'mBOCD'});
xtickangle(30)
xlim([0 9])
ylabel('PPV');
ylim([0 1]);
legend(bp3,'location','best');
legend('boxoff');
set(gca, 'fontsize', 16);

matrix_fp = [[results_rmdm.fp]', ...
              [results_bis.fp]', ...
              [results_pelt1.fp]', ...              
              [results_pelt2.fp]', ...
              [results_bblocks.fp]', ...
              results_bcp.fp, ...             
              [results_bocdo.fp]', ...
              [results_bocd.fp]'];

subplot(2,2,4);
bp4 = bplot(matrix_fp, 'outliers');
xlabel('\bf (d)');
xticks(1:8);
xticklabels({'RMDM', 'BiS', 'PELT1', 'PELT2', 'BBlocks', 'BCP', 'BOCD', 'mBOCD'});
xtickangle(30)
xlim([0 9])
ylabel('False Positive Count');
legend(bp4,'location','best');
legend('boxoff');
set(gca, 'fontsize', 16);

tightfig;

%end