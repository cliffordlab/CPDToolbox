function seg = rmdm_v1(sig, P_0)

% Description:
%   Recursive mean difference maximization by Pedro Bernaola-Galvan.
%
% Bernaola-GalvÃ¡n, P., Ivanov, P., Nunes Amaral, L., & Stanley, H. (2001).
% Scale Invariance in the Nonstationarity of Human Heart Rate.
% Physical Review Letters, 87(16), 168105.
%
% Arguments:
%   ts_in - Input data timeseries
%   P_O - significance threshold
%
% Results (all optional):
%   ts_out - Timeseries of data segments
%
%   Dependencies
%       https://github.com/cliffordlab/change_point_detection
%
%   Reference(s)
%
%   Authors
%       Maxim Osipov
%   
%   Copyright (C) 2017 Authors
%   All rights reserved.
%   
%   This software may be modified and distributed under the terms
%   of the BSD license.  See the LICENSE file in this repo for details.


%sig = ts_in.Data;
seg = [];
seg_idx = 1;

% queue to hold intervals to be segmented
q{1} = [1 length(sig)];
q_idx = 1;

while q_idx <= length(q)
    % re-initialize
    clear t;

    % move a sliding pointer from left to right of the interval
    sig_start = q{q_idx}(1);
    sig_end = q{q_idx}(2);
    t_idx = 1;
    for i=sig_start:sig_end
        % at each position measure significance of difference between means
        t(t_idx) = func_t(sig(sig_start:i), sig((i+1):sig_end));
        t_idx = t_idx + 1;
    end

    % find maximum of t (t_max) and split point (t_idx, sig_idx in complete
    % signal)
    [t_max, t_idx] = max(t);
    sig_idx = t_idx + sig_start - 1;

    % find significance of the maximum
    p = func_p(P_0, sig_end-sig_start+1);

    % if significance exceeds limit, repeat the above for both sub-signals
    if p <= (1-P_0)
        % add here checking for the neighbouring segments
        if (sig_idx > sig_start+1) && (sig_idx < sig_end-1)
            % add change point
            seg(seg_idx) = sig_idx;
            seg_idx = seg_idx + 1;
            % queue sub-segments
            q{length(q)+1} = [sig_start sig_idx];
            q{length(q)+1} = [sig_idx sig_end];
        end
    end

    %  if inner segment
    %   check if significance between left and lefter and right and righter
    %   is also significant
    q_idx = q_idx + 1;
end

seg = sort(seg);
%ts_out = timeseries(ts_in.Time(seg(1:end)), ts_in.Time([1 seg(1:end-1)]));

end

% compute the significance level
function p = func_p(tau, N)
    nu = N-2;
    gamma = 4.19*log(N) - 11.54;
    delta = 0.4;
    p = (1 - betainc(nu/(nu+tau^2), delta*nu, delta))^gamma;
end

% compute the difference statistic
function t = func_t(sig_left, sig_right)
    N_left = length(sig_left);
    N_right = length(sig_right);
    mu_left = mean(sig_left);
    mu_right = mean(sig_right);
    s_left = std(sig_left);
    s_right = std(sig_right);

    s_d = sqrt( (s_left^2 + s_right^2) / (N_left + N_right - 2) )*...
            sqrt(1/N_left + 1/N_right);
    t = abs((mu_left - mu_right)/s_d);
end

