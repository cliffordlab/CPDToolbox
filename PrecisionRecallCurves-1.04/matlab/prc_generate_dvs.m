% Generates an imaginary output of a classifier, given a parameterization
% of the distributions underlying its decision values.
%
% Usage:
%     [targs, dvs] = prc_generate_dvs(n, alpha, muN, sigmaN, muP, sigmaP)
%
% Arguments:
%     n: number of examples
%     alpha: fraction of examples from the positive class (range: 0..1)
%         0   = total class imbalance in favour of the negative class
%         0.5 = no class imbalance
%         1   = total class imbalance in favour of the positive class
%     muN, sigmaN: mean and std.dev. of decision values assigned to
%         excamples from negative class 
%     muP, sigmaP: mean and std.dev. of decision values assigned to
%         examples from positive class 
% 
% Return values:
%     targs: true class labels (targets)
%     dvs: decision values

% Kay H. Brodersen & Cheng Soon Ong, ETH Zurich, Switzerland
% $Id: prc_generate_dvs.m 6393 2010-06-14 15:24:46Z bkay $
% -------------------------------------------------------------------------
function [targs, dvs] = prc_generate_dvs(n, alpha, muN, sigmaN, muP, sigmaP)
    
    % Check input
    assert(sigmaN>0 && sigmaP>0, 'variances must be strictly positive');
    assert(0<alpha && alpha<1, 'alpha must be strictly between 0 and 1');
    
    % Generate negative class
    n1 = round((1-alpha)*n);
    if n1<2, warning('at least 2 values from either class recommended'), end
    targs = -ones(1,n1);
    dvs = normrnd(muN, sigmaN, 1, n1);
    
    % Generate positive class
    n2 = n-n1;
    if n2<2, warning('at least 2 values from either class recommended'), end
    targs = [targs, ones(1,n2)];
    dvs = [dvs, normrnd(muP, sigmaP, 1, n2)];
    
    % Check return value
    assert(size(targs,1)==1);
    assert(size(dvs,1)==1);
    assert(size(targs,2)==size(dvs,2));
    assert(size(targs,2)==n);
    
end
