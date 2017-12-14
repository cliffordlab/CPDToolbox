function ecp = bblocks(rr, ncp_prior)
%   
%   OVERVIEW
%      Bayesian Blocks Changepoint Detection:
%      - Finds optimal segmentation of "Point Measurement" data in the 
%      observation interval
%     
%   REFERENCE
%   Scargle, J.D., Norris, J.P., Jackson, B. and Chiang, J., 2013. Studies in 
%   astronomical time series analysis. vi. bayesian block representations. The 
%   Astrophysical Journal, 764(2), p.167.
%       
%   INPUT
%       data_in:
%       data_in.cell_data - Nx2 array of measurement value and time
%       data_in.fp_rate - false positives rate (optional)
%       data_in.ncp_prior - log gamma for the number of blocks prior (optional)
%       data_in.do_iter - number of iterations to estimate number of blocks (optional)
%
%   OUTPUT
%       data_out:
%       data_out.change_points 
%       data_out.num_vec
%       data_out.rate_vec
%       data_out.best
%       data_out.last
%       data_out.ncp_prior
%
%   AUTHORS
%   Jeffrey D. Scargle (original code)
%	Ayse Selin Cakmak <acakmak3@gatech.edu> (modified function)
%   
%   COPYRIGHT (C) 2017 AUTHORS (see above)
%       This code (and all others in this repository)
%       are under the GNU General Public License v3.
%       See LICENSE in this repo for details.

% Convert input data into strange structure format
cell_data = zscore(rr);
cell_data(:, 2) = ones(length(rr), 1);
data_in.cell_data = cell_data;

[ num_points, dummy ] = size( cell_data );
tt = 1: num_points; % nominal evenly spaced time points
nn_vec = [];

if isfield( data_in, 'fp_rate')
    fp_rate = data_in.fp_rate;
else
    fp_rate = .05; % Default value
end

%{
if isfield( data_in, 'ncp_prior')
    ncp_prior = data_in.ncp_prior;
else
   % ncp_prior = 4 - log( fp_rate / ( 0.0136 * num_points .^ (0.478 ) ) );
   
   % This formula comes from the end of section 3 of the referenced
   % article.
   ncp_prior = 1.32 + 0.577*log10( num_points );
end
%}

if isfield( data_in, 'do_iter')
   do_iter = data_in.do_iter;
else
    do_iter = 0; % Default: do not iterate on ncp_prior
end

%===========================================
change_points = [];
count_vec = [];

iter_count = 0;
iter_max = 10;
while 1
    
    best = [];
    last = [];
    cpu_0 = cputime;
    
    for R = 1:num_points 

% Calculate a and b from formulas (31) and (32)
        sum_x_1 = cumsum( cell_data( R:-1:1, 1 ) )'; % sum( x    / sig^2 )
        sum_x_0 = cumsum( cell_data( R:-1:1, 2 ) )'; % sum[ 1    / sig^2 )
        
% Calculate the fitness function by using formula (41)
        fit_vec = ( ( sum_x_1(R:-1:1) ) .^ 2 ) ./ ( 2* sum_x_0(R:-1:1));

        [ best(R), last(R)] = max( [ 0 best ] + fit_vec - ncp_prior );     
    end

    %----------------------------------------------------------------------
    % Now find changepoints by iteratively peeling off the last block
    %----------------------------------------------------------------------
    
    index = last( num_points );
    change_points = [];

    while index > 1
        change_points = [ index change_points ];
        index = last( index - 1 );
    end
    
    %----------------------------------------------------------------------
    %            Iterate if desired - Optional
    %----------------------------------------------------------------------
    if do_iter == 0
        break
    else
        iter_count = iter_count + 1;
        num_cp = length( change_points );
        if num_cp < 1
            num_cp = 1;
        end

        if exist('cpt_old' )

            if num_cp == length( cpt_old ) % compare with previous iteration
                err_this = sum( abs( change_points - cpt_old ) );
            else
                err_this = Inf;
            end

            if err_this == 0
                fprintf(1,'Converged at %3.0f\n', iter_count )
                break
            end

            if iter_count > iter_max
                fprintf(1,'Did not converge at %3.0f\n', iter_count )
                break
            end

        end

        fp_rate = 1 - ( 1 - fp_rate ) .^ ( 1 / num_cp );
        ncp_prior = 1.32 + 0.577*log10( num_points );
        cpt_old    = change_points;
        
    end
    
end

num_changepoints = length( change_points );
num_blocks = num_changepoints + 1;

 rate_vec = zeros( num_blocks, 1 );
  num_vec = zeros( num_blocks, 1 );
   dt_vec = zeros( num_blocks, 1 );
 tt_1_vec = zeros( num_blocks, 1 );
 tt_2_vec = zeros( num_blocks, 1 );

cpt_use = [ 1 change_points ];

for id_block = 1: num_blocks
    
    ii_1 = cpt_use( id_block ); % start
    if id_block < num_blocks
        ii_2 = cpt_use( id_block + 1 ) - 1;
    else
        ii_2 = num_points;
    end
    
    xx_this = cell_data( ii_1:ii_2, 1);
    wt_this = cell_data( ii_1:ii_2, 2); 
    rate_vec( id_block ) = sum( wt_this .* xx_this ) / sum( wt_this );
    
end

% Assign output values
% data_out.change_points = change_points;
% data_out.num_vec       = num_vec;
% data_out.rate_vec      = rate_vec;
% data_out.best = best;
% data_out.last = last;
% data_out.ncp_prior     = ncp_prior;
% data_out.nn = nn_vec;

% Return the only thing we care about right now
ecp = change_points';
    
end

