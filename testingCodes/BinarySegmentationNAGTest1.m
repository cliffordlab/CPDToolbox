%% NAG Binary Segmentation

fprintf('g13nd example results\n\n');

% Input series
y = [ 0.00; 0.78;-0.02; 0.17; 0.04;-1.23; 0.24; 1.70; 0.77; 0.06;
      0.67; 0.94; 1.99; 2.64; 2.26; 3.72; 3.14; 2.28; 3.78; 0.83;
      2.80; 1.66; 1.93; 2.71; 2.97; 3.04; 2.29; 3.71; 1.69; 2.76;
      1.96; 3.17; 1.04; 1.50; 1.12; 1.11; 1.00; 1.84; 1.78; 2.39;
      1.85; 0.62; 2.16; 0.78; 1.70; 0.63; 1.79; 1.21; 2.20;-1.34;
      0.04;-0.14; 2.78; 1.83; 0.98; 0.19; 0.57;-1.41; 2.05; 1.17;
      0.44; 2.32; 0.67; 0.73; 1.17;-0.34; 2.95; 1.08; 2.16; 2.27;
     -0.14;-0.24; 0.27; 1.71;-0.04;-1.03;-0.12;-0.67; 1.15;-1.10;
     -1.37; 0.59; 0.44; 0.63;-0.06;-0.62; 0.39;-2.63;-1.63;-0.42;
     -0.73; 0.85; 0.26; 0.48;-0.26;-1.77;-1.53;-1.39; 1.68; 0.43];

% Type of change point(s) being looked for
% (change in mean, assuming a Normal distribution)
ctype = int64(1);

% Standard deviation to use for Normal distribution
param = 1;

[ecp_nag ,sparam, ifail] = g13nd(ctype, y, 'param', param);

% Print the results
fprintf('  -- Change Points --         --- Distribution ---\n');
fprintf('  Number     Position              Parameters\n');
fprintf(' ==================================================\n');
for i = 1:numel(ecp_nag)
  fprintf('%5d%13d%16.2f%16.2f\n', i, ecp_nag(i), sparam(1:2,i));
end

% Plot the results
fig1 = figure;

% Plot the original series
plot(y,'Color','red');

% Mark the change points, drop the last one as it is always
% at the end of the series
xpos = transpose(double(ecp_nag(1:end-1))*ones(1,2));
ypos = diag(ylim)*ones(2,numel(ecp_nag)-1);
line(xpos,ypos,'Color','black');

% Plot the estimated mean in each segment
xpos = transpose(cat(2,cat(1,1,ecp_nag(1:end-1)),ecp_nag));
ypos = ones(2,1)*sparam(1,:);
line(xpos,ypos,'Color','green');

% Add labels and titles
title({'{\bf g13nd Example Plot}', ...
       'Simulated time series and corresponding changes in mean'});
xlabel('{\bf Time}');
ylabel('{\bf Value}');

%% My Binary Segmentation

ecp = binarySegmentation(y, 1, 'BIC');




