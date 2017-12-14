% Simple demonstration of different methods for estimating a
% precision-recall curve (PRC).
% 
% Generates decision values from a binormal distribution and plots true vs.
% empirical vs. model-based quantities derived from the confusion matrix of
% the simulated classifier.
% 
% Usage:
%     prc_demo
% 
% Literature:
%     K.H. Brodersen, C.S. Ong, K.E. Stephan, J.M. Buhmann (2010). The
%     binormal assumption on precision-recall curves. In: Proceedings of
%     the 20th International Conference on Pattern Recognition (ICPR).

% Kay H. Brodersen & Cheng Soon Ong, ETH Zurich, Switzerland
% $Id: prc_demo.m 5529 2010-04-22 21:10:32Z bkay $
% -------------------------------------------------------------------------

% Specify distributions of decision values
n = 100;                 % number of predictions
alpha = 0.4;             % fraction of positive decision values
muN = -1; sigmaN = 2;    % mean and std of negative class
muP = 1; sigmaP = 2;     % mean and std of positive class

% Generate classifier output
[targs, dvs] = prc_generate_dvs(n, alpha, muN, sigmaN, muP, sigmaP);


% -------------------------------------------------------------------------
% Compute empirical curves
[TPR_emp, FPR_emp, PPV_emp] = prc_stats_empirical(targs, dvs);

% Compute true curves
[TPR_true, FPR_true, PPV_true] = prc_stats(alpha, muN, sigmaN, muP, sigmaP);

% Compute smooth curves (binormal model)
[TPR_bin, FPR_bin, PPV_bin] = prc_stats_binormal(targs, dvs, false);

% Compute smooth curves (alpha-binormal model)
[TPR_abin, FPR_abin, PPV_abin] = prc_stats_binormal(targs, dvs, true);


% -------------------------------------------------------------------------
% Plot results
cols = [200 45 43; 37 64 180; 0 176 80; 0 0 0]/255;

% Plot decision-value histograms
% (red: negative class; green: positive class)
figure; hold on;
binWidth = 0.5;
pdfRange = [muN-3*sigmaN:0.01:muP+3*sigmaP];
prc_conthist(dvs(targs==-1), binWidth, cols(1,:), 2);
prc_conthist(dvs(targs==+1), binWidth, cols(3,:), 2);
plot(pdfRange, normpdf(pdfRange, muN, sigmaN)*binWidth*sum(targs==-1), 'linewidth', 2, 'color', cols(1,:));
plot(pdfRange, normpdf(pdfRange, muP, sigmaP)*binWidth*sum(targs==+1), 'linewidth', 2, 'color', cols(3,:));
xlabel('decision value'); ylabel('density/frequency'); title('Decision values');
set(gca, 'box', 'on');
axis tight;

% Plot ROC curves
figure; hold on;
plot(FPR_true, TPR_true, '-', 'color', cols(4,:), 'linewidth', 2);
plot(FPR_emp, TPR_emp, '-o', 'color', cols(1,:), 'linewidth', 2);
plot(FPR_bin, TPR_bin, '-', 'color', cols(2,:), 'linewidth', 2);
plot(FPR_abin, TPR_abin, '--', 'color', cols(3,:), 'linewidth', 2);
plot([0 1], [0 1], '-', 'color', [0.7 0.7 0.7], 'linewidth', 2);
axis([0 1 0 1]);
xlabel('FPR'); ylabel('TPR'); title('ROC curves');
set(gca, 'box', 'on');

% Plot PR curves
figure; hold on;
plot(TPR_true, PPV_true, '-', 'color', cols(4,:), 'linewidth', 2);
plot(TPR_emp, PPV_emp, '-o', 'color', cols(1,:), 'linewidth', 2);
plot(TPR_bin, PPV_bin, '-', 'color', cols(2,:), 'linewidth', 2);
plot(TPR_abin, PPV_abin, '-', 'color', cols(3,:), 'linewidth', 2);
axis([0 1 0 1]);
xlabel('TPR (recall)'); ylabel('PPV (precision)'); title('PR curves');
set(gca, 'box', 'on');

