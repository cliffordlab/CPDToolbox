 function [rr, tcp, sleepStart, sleepEnd] = rrgenV3_wrapper(seed, vector_length, ...
                                        prob_ectopy, prob_noise, pathToRepo, act)
%   
%   OVERVIEW:
%       Matlab wrapper for calling compiled rrgen executable
%
%       First, run:
%           gcc -Wall rrgenV3.c -lm -o rrgenV3
%   
%   REFERENCE
%
%   INPUT
%       seed          - number to initialize random number generator
%       vector_length - length of RR output
%       prob_ectopy   - probability of a beat being ectopic
%       prob_noise    - probability of a beat being noisy
%       act           - option to work with actigraphy data
%
%   OUTPUT
%       rr - RR interval time series
%       tcp - True change points in time series
%       sleepStart - Index of sleep start in artificial RR interval code
%       sleepEnd - Index of sleep end in artificial RR interval code
%
%   REPO
%       https://github.com/cliffordlab/rrgen
%  
%   REFERENCES
%   Please refer to C code this function calls for more details:
%   https://github.com/cliffordlab/CPDToolbox/blob/master/artificialDataAnalysis/rrGenerator/rrgenV3.c
%

% Initialize constants for rrgen
if nargin < 4
    prob_noise = 0.0048;
end

if nargin < 3
    prob_ectopy = 0.0003;
end

if nargin < 2
    seed = now; % Create seed from current unix timestamp
end

if nargin < 1
    vector_length = 500;
end

if vector_length < 200
    vector_length = 200;
end

generate = 1;
while generate == 1
    command = [pathToRepo filesep 'artificialDataAnalysis' filesep 'rrGenerator' filesep 'rrgenV3' ' ' num2str(seed) ' ' num2str(vector_length) ' ' num2str(prob_ectopy) ' ' num2str(prob_noise)];

    % Call rrgen2 via system and save output
    [~, systemout] = system(command);

    % Read the system output into cells
    systemout_cell = textscan(systemout,'%f %f %f','Delimiter',',');

    % Convert cells of doubles into matrix
    rrgen_mat = cell2mat(systemout_cell);

    % Isolate columns of matrix into output vectors
    trpeaks = rrgen_mat(rrgen_mat(:,3)==7,1);
    rr = rrgen_mat(rrgen_mat(:,3)~=7,1);

    % Output true change points depending on noise probability
    tcp_idx = find(rrgen_mat(rrgen_mat(:,3)~=7,2)==1);

    % Find sleep start and end indexes
    sleep_idx = find(rrgen_mat(rrgen_mat(:,3)~=7,3)==0);

    if isempty(sleep_idx)
        disp('Regenerating RR intervals...');
    else 
        generate = 0; % Artificial data generation complete;
    end
end

% Adjust for noise if necessary
tcp = zeros(length(tcp_idx),1);
if prob_noise==0
    tcp = tcp_idx;
else
    tcp_time = trpeaks(tcp_idx);
    t = cumsum(rr);
    for i=1:length(tcp)
        [~,I] = min(abs(t-tcp_time(i)));
        tcp(i) = I(1);
    end
end

sleep = zeros(length(sleep_idx),1);
if prob_noise==0
    sleep = sleep_idx;
else
    sleep_time = trpeaks(sleep_idx);
    t = cumsum(rr);
    for i=1:length(sleep)
        [~,I] = min(abs(t-sleep_time(i)));
        sleep(i) = I(1);
    end
end

if act ~= 1
    sleepStart = sleep(1);
    sleepEnd = sleep(end-1);
else
    sleepStart = [];
    sleepEnd = [];
end

end 