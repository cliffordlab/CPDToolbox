

%% Generate artificial time series
r1 = normrnd(3,10,[1,40]);
r2 = normrnd(3,2,[1,80]);
r3 = normrnd(3,40,[1,80]);
r4 = normrnd(3,10,[1,100]);
y = [r1(:); r2(:); r3(:); r4(:)];
tcp = [40; 120; 200; 300];

%% NAG Binary Segmentation
 
% Type of change point(s) being looked for
% (change in mean, assuming a Normal distribution)
ctype = int64(2);
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
line(xpos,ypos,'Color','black','linewidth',2);

% Plot the estimated mean in each segment
xpos = transpose(cat(2,cat(1,1,ecp_nag(1:end-1)),ecp_nag));
ypos = ones(2,1)*sparam(1,:);
line(xpos,ypos,'Color','green','linewidth',2);

% Add labels and titles
title({'{\bf g13nd Example Plot}', ...
       'Simulated time series and corresponding changes in mean'});
xlabel('{\bf Time}');
ylabel('{\bf Value}');

%% My Binary Segmentation

ecp = binarySegmentation(y, 2, 'BIC');




