function run_length = mbocd(rr, lambda, prior_type)
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
%
%       Murphy K. Conjugate Bayesian analysis of the Gaussian distribution.
%       https://www.cs.ubc.ca/~murphyk/Papers/bayesGauss.pdf
%       
%   INPUT
%       rr      - Input data timeseries
%       lambda  - Timescale or expected run length
%       prior_type - Type of prior: should be "chi" for inverse chi-squared
%       distribution, "gamma" for inverse gamma distribution
%
%   OUTPUT
%       prob_run_lengths - Probability of run lengths
%       run_length
%
%   DEPENDENCIES
%
%   AUTHORS
%       Ryan P. Adams (original code)
%       Ayse Cakmak <acakmak3@gatech.edu> (modified code)
%   
%   COPYRIGHT (C) 2017 AUTHORS (see above)
%       This code (and all others in this repository)
%       are under the GNU General Public License v3.
%       See LICENSE in this repo for details.

% Turn input 'rr' into a timeseries object
%t = cumsum(rr);
%ts = timeseries(rr, t);

% Interpolate and resample data
%ts.DataInfo.Interpolation = tsdata.interpolation('zoh');
%tsrs = resample(ts, ts.Time);

%x = tsrs.Data;

x = rr;
x(isnan(x)) = nanmean(x); % Replace NaNs if any exist
T = length(x);

% Calculate initial parameters of prior:
% https://www.cs.ubc.ca/~murphyk/Teaching/CS340-Fall07/reading/NG.pdf
x      = zscore(x); % z-score the data
sigma0 = std(x);	% standard deviation; should be 1 after zscore
var0   = sigma0^2;	% variance; should be 1 after zscore
tau    = 1 / var0;	% precision is inverse of variance; should be 1 after zscore
mu0    = 0;         % normal prior on mu; should be 0 after zscore
kappa0 = 1;         % model is insensitive to kappa0

%%%%%% Modified by Ayse Cakmak
if strcmp(prior_type,'gamma')
    % Uninformative Priors
    alpha0 = 1e-4;         
    beta0  = 1e-4;         

    muT    = mu0;
    kappaT = kappa0;
    alphaT = alpha0;
    betaT  = beta0;
elseif strcmp(prior_type,'chi')
    v0 = 6;  
    muT    = mu0;
    kappaT = kappa0;
    vT = v0;
    varT = var0;
end
%%%%%% End of mofied part

% Matrix that holds probabilities of run lengths at every step,
% at t=1 the run length is zero with 100% probability.
prob_run_lengths = zeros([T+1 T]);
prob_run_lengths(1,1) = 1;

% Keep track of the maximums.
run_length  = zeros([T+1]);

% Loop over the data like we're seeing it all for the first time.
for t=1:T
    
    %%%%%% Modified by Ayse Cakmak
    % Estimate posterior predictive    
    % Parameterize student PDF by variance
    if strcmp(prior_type,'gamma')
        precision = (alphaT .* kappaT) ./ (betaT .* (kappaT + 1)); % page 10, near eqn 110
        variance = 1 ./ precision;
        nu = 2 * alphaT; % degrees of freedom    
    elseif strcmp(prior_type,'chi')
        variance = varT;
        nu = vT;
    end 
    %%%%%% End of mofied part

    % Update according to calculated values
    predprobs = studentpdf(x(t), muT, variance, nu);
    
    % Evaluate the hazard function for this interval;
    % hazard function is constant
    H = 1/lambda;
    
    % Evaluate the growth probabilities; shift the probabilities down and
    % to the right, scaled by the hazard function and the predictive probabilities.
    prob_run_lengths(2:t+1, t+1) = prob_run_lengths(1:t, t) .* predprobs .* (1-H);
    
    % Evaluate the probability that there *was* a changepoint and we're
    % accumulating the mass back down at r = 0.
    prob_run_lengths(1, t+1) = sum(prob_run_lengths(1:t, t) .* predprobs .* H);
    
    % Renormalize the run length probabilities for improved numerical stability.
    prob_run_lengths(:, t+1) = prob_run_lengths(:, t+1) ./ sum(prob_run_lengths(:, t+1));

    % Update the parameter sets for each possible run length.
    if strcmp(prior_type,'gamma')
        muT0    = [ mu0    ; (kappaT .* muT + x(t)) ./ (kappaT+1) ];
        kappaT0 = [ kappa0 ; kappaT + 1 ];
        alphaT0 = [ alpha0 ; alphaT + 0.5 ];
        betaT0  = [ beta0  ; betaT + (kappaT .* (x(t) - muT).^2) ./ (2*(kappaT + 1)) ];

        muT     = muT0;
        kappaT  = kappaT0;
        alphaT  = alphaT0;
        betaT   = betaT0;
        
    %%%%%% Modified by Ayse Cakmak    
    elseif strcmp(prior_type,'chi')
        muT0    = [ mu0    ; (kappaT .* muT + x(t)) ./ (kappaT+1) ];
        vT0 = [ v0 ; vT + 1 ];
        varT0 = [var0; vT.*varT + (kappaT .* (x(t)-muT).^2) ./ (kappaT+1) ]./vT0;
        kappaT0 = [ kappa0 ; kappaT + 1 ];
        
        muT     = muT0;
        kappaT  = kappaT0;
        vT  = vT0;
        varT   = varT0;
    end
    
    % Run length calculation - New Idea
    if t == 1
        run_length(t) = 1;
        rl_dummy = 2;
    else
        % Calculate probabilities of continuing the run or starting a new
        % run
        new_run = prob_run_lengths(1,t);
        cont_run = sum(prob_run_lengths(1:t,t) .* predprobs .* (1-H));

        %Compare and fix run length according to this comparison
        if cont_run > new_run
            run_length(t) = rl_dummy;
            rl_dummy = rl_dummy+1;
        else
            run_length(t) = 1;
            rl_dummy = 2;
        end
    end
    %%%%%% End of mofied part
end

end

%% Other functions %%

% This form is taken from Kevin Murphy's lecture notes:
% https://www.cs.ubc.ca/~murphyk/Teaching/CS340-Fall07/reading/NG.pdf
function p = studentpdf(x, mu, var, nu)
    % Original expression does not work well for larger values of nu
    % c = gamma(nu/2 + 0.5) ./ gamma(nu/2) ./ sqrt(nu .* pi .* var);
    
    % Equivalent expression (exp of log gamma) which works better in Matlab
    c = exp(gammaln(nu/2 + 0.5) - gammaln(nu/2)) .* (nu .* pi .* var).^(-0.5);
    
    p = c .* (1 + (1./(nu.*var)).*(x-mu).^2).^(-(nu+1)/2);
end

% Calculates hazard function
function p = constant_hazard(r, lambda)
    p = 1/lambda * ones(size(r));
end