function [prob_run_lengths, maxes] = bocd(rr, lambda)
%   
%   OVERVIEW
%      Bayesian online changepoint detection:
%      Demonstration of online detection of a change in 1D Gaussian parameters,
%      updated for incremental data updates. Instead of proposed elimination
%      of minor probabilities, just limit the maximum run length and restart if exceeded.
%   
%   REFERENCE
%       Adams RP, MacKay DJC. Bayesian Online Changepoint Detection.
%       arXiv. 2007:7. doi:arXiv:0710.3742v1
%       Original code: http://hips.seas.harvard.edu/content/bayesian-online-changepoint-detection
%       
%   INPUT
%       rr      - Input data timeseries
%       lambda  - Timescale or expected run length
%
%   OUTPUT
%       prob_run_lengths - Probability of run lengths
%       maxes - 
%
%   DEPENDENCIES
%
%   AUTHORS
%   Ryan P. Adams (original code)
%	Erik Reinertsen <er@gatech.edu> (modularized function)
%   
%   COPYRIGHT (C) 2017 AUTHORS (see above)
%       This code (and all others in this repository)
%       are under the GNU General Public License v3.
%       See LICENSE in this repo for details.

% Turn input 'rr' into a timeseries object
t = cumsum(rr);
ts = timeseries(rr, t);

% Interpolate and resample data
ts.DataInfo.Interpolation = tsdata.interpolation('zoh');
tsrs = resample(ts, ts.Time);

x = tsrs.Data;
x(isnan(x)) = nanmean(x); % Replace NaNs if any exist
T = length(x);

% Calculate initial parameters of prior assuming normal gamma distribution:
% https://www.cs.ubc.ca/~murphyk/Teaching/CS340-Fall07/reading/NG.pdf
x      = zscore(x); % z-score the data
sigma0 = std(x);	% standard deviation; should be 1 after zscore
var0   = sigma0^2;	% variance; should be 1 after zscore
tau    = 1 / var0;	% precision is inverse of variance; should be 1 after zscore
mu0    = 0;         % normal prior on mu; should be 0 after zscore
kappa0 = 1;         % model is insensitive to kappa0
alpha0 = 1;         % if zscore, should be 1
beta0  = 1e-4;         % if zscore, should be 1

muT    = mu0;
kappaT = kappa0;
alphaT = alpha0;
betaT  = beta0;

% Specify the hazard function.
% This is a handle to a function that takes one argument - the number of
% time increments since the last changepoint - and returns the probability
% of changepoint.  Generally you might want to have your hazard function
% take parameters, so using an anonymous function is helpful.  We're going
% to just use the simple constant-rate hazard function that gives
% geomtrically-drawn intervals between changepoints.
% We'll specify the rate via a mean.
hazard_func  = @(r) constant_hazard(r, lambda);

% Matrix that holds probabilities of run lengths at every step,
% at t=1 the run length is zero with 100% probability.
prob_run_lengths = zeros([T+1 T]);
prob_run_lengths(1,1) = 1;

% Keep track of the maximums.
maxes  = zeros([T+1]);
% Loop over the data like we're seeing it all for the first time.
for t=1:T
    
    % Estimate posterior predictive
    % In the case of m = 1 (new observations), this follows the Student's T-distribution
	% t_nu = studentpdf(x, mu, var, nu)
    variance = betaT .* (kappaT+1) ./ (alphaT .* kappaT);
    nu = 2 * alphaT; % degrees of freedom
    predprobs = studentpdf(x(t), muT, variance, nu);
    
    % Evaluate the hazard function for this interval.
    H = hazard_func([1:t]');
    
    % Evaluate the growth probabilities; shift the probabilities down and
    % to the right, scaled by the hazard function and the predictive probabilities.
    prob_run_lengths(2:t+1,t+1) = prob_run_lengths(1:t,t) .* predprobs .* (1-H);
    
    % Evaluate the probability that there *was* a changepoint and we're
    % accumulating the mass back down at r = 0.
    prob_run_lengths(1,t+1) = sum( prob_run_lengths(1:t,t) .* predprobs .* H );
    
    % Renormalize the run length probabilities for improved numerical stability.
    prob_run_lengths(:,t+1) = prob_run_lengths(:,t+1) ./ sum(prob_run_lengths(:,t+1));

    % Update the parameter sets for each possible run length.
    muT0    = [ mu0    ; (kappaT .* muT + x(t)) ./ (kappaT+1) ];
    kappaT0 = [ kappa0 ; kappaT + 1 ];
    alphaT0 = [ alpha0 ; alphaT + 0.5 ];
    betaT0  = [ beta0  ; betaT + (kappaT .* (x(t)-muT).^2) ./ (2*(kappaT+1)) ];
    
    muT     = muT0;
    kappaT  = kappaT0;
    alphaT  = alphaT0;
    betaT   = betaT0;
    
    % Store the maximum to plot later.
    maxes(t) = find(prob_run_lengths(:,t) == max(prob_run_lengths(:,t)));
end

end % end bocp.m

%% Other functions %%

% This form is taken from Kevin Murphy's lecture notes:
% https://www.cs.ubc.ca/~murphyk/Teaching/CS340-Fall07/reading/NG.pdf
function p = studentpdf(x, mu, var, nu)
    % original below; mathematically equivalent to Murphy's notes
    c = exp(gammaln(nu/2 + 0.5) - gammaln(nu/2)) .* (nu .* pi .* var).^(-0.5);
    
    % Erik's rewrite: identical to expression in Murphy's notes, but
    % does not work well for larger values of nu because values get too large for Matlab
%     c = gamma(nu/2 + 0.5) ./ gamma(nu/2) ./ sqrt(nu .* pi .* var);

    p = c .* (1 + (1./(nu.*var)).*(x-mu).^2).^(-(nu+1)/2);
end

% Calculates hazard function
function p = constant_hazard(r, lambda)
    p = 1/lambda * ones(size(r));
end
