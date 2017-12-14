function ecp = rmdm(data, P_0, l0)
%
% Overview:
%   Recursive mean difference maximization by Pedro Bernaola-Galvan.
%
% Input:
%   data - Input data timeseries
%   P_O - significance threshold
%   l0 - min sample number between ECPs
%
%   Dependencies
%       https://github.com/cliffordlab/change_point_detection
%
%   Reference(s)
%       Bernaola-GalvÃ¡n, P., Ivanov, P., Nunes Amaral, L., & Stanley, H. (2001).
%       Scale Invariance in the Nonstationarity of Human Heart Rate.
%       Physical Review Letters, 87(16), 168105.
%
%   Authors
%       Maxim Osipov (Original code)
%       Ayse Selin Cakmak (Modified code)
%   
%   Copyright (C) 2017 Authors
%   All rights reserved.
%   
%   This software may be modified and distributed under the terms
%   of the BSD license.  See the LICENSE file in this repo for details.

%sig = ts_in.Data;
ecp = [1; length(data)];
seg_idx = 3;

% queue to hold intervals to be segmented
clear q;
q{1} = [1 length(data)];
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
        t(t_idx) = func_t(data(sig_start:i), data((i+1):sig_end));
        t_idx = t_idx + 1;
    end

    % find maximum of t (t_max) and split point (t_idx, sig_idx in complete
    % signal)
    [t_max, t_idx] = max(t);
    sig_idx = t_idx + sig_start - 1;

    % find significance of the maximum
    p = func_p(t_max, sig_end-sig_start+1);

    %%%%%% Modified by Ayse Cakmak    
    accept = 0;  
    % if significance exceeds limit, repeat the above for both sub-signals
    if p > (P_0)
        if (sig_idx > sig_start+l0) && (sig_idx < sig_end-l0)
           
            % Check if there is cp on the right
            if length(ecp(ecp>sig_idx))>1
                % Check the significance of new seg and the one on the
                % right
                new = sig_idx:sig_end;
                e = find(ecp>sig_idx);
                e = ecp(e(2));
                right = sig_end+1:e;
                t1 = func_t(data(new), data(right));
                p1 = func_p(t1,length(new)+length(right));
                if p1 > (P_0)
                    accept = accept+1;
                end
            else
                accept = accept+1;
            end
            
            if length(ecp(ecp<sig_idx))>1
                % Check the significance of new seg and the one on the
                % left
                new = sig_start:sig_idx;
                e = find(ecp<sig_idx);
                e = ecp(e(end-1));
                left = e:sig_start-1;
                t2 = func_t(data(left), data(new));
                p2 = func_p(t2,length(new)+length(left));
                if p2 > (P_0)
                    accept = accept+1;
                end                    
            else
                accept = accept+1;
            end
            
            if accept==2
                % add change point
                ecp(seg_idx) = sig_idx;
                ecp = sort(ecp);
                seg_idx = seg_idx + 1;
                % queue sub-segments
                q{length(q)+1} = [sig_start sig_idx];
                q{length(q)+1} = [sig_idx sig_end];
            end
        end
    end
    %%%%%% End of mofied part

    %  if inner segment
    %   check if significance between left and lefter and right and righter
    %   is also significant
    q_idx = q_idx + 1;
end

ecp = sort(ecp);

end

% compute the significance level
function g = func_p(tau, N)
    nu = N-2;
    gamma = 4.19*log(N) - 11.54;
    delta = 0.4;
    p = (1 - betainc(nu/(nu+tau^2), delta*nu, delta));
    g = p^gamma;
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

