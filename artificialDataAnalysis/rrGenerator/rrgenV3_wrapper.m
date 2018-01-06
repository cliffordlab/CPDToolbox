 function [rr, tcp] = rrgenV3_wrapper(seed, vector_length, prob_ectopy, prob_noise, pathToRepo)
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
%
%   OUTPUT
%       rr - RR interval time series
%       tcp - True change points in time series
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

end 