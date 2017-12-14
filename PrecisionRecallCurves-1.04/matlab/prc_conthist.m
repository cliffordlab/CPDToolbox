% Plots a contour histogram.
% 
% Usage:
%     prc_conthist(x, binWidth, color, linewidth)
% 
% Arguments:
%     x: data
%     binWidth: width of each histogram bin
%     color: line color of contour
%     linewidth: width of contour

% Kay H. Brodersen, ETH Zurich, Switzerland
% $Id: prc_conthist.m 5516 2010-04-22 17:17:24Z bkay $
% -------------------------------------------------------------------------
function prc_conthist(x, binWidth, color, linewidth)
    
    edges = floor(min(x)/binWidth)*binWidth:binWidth:ceil(max(x)/binWidth)*binWidth;
    v = sort(repmat(edges,1,2));
    n = histc(x, edges);
    n = [0, reshape(repmat(n(1,1:end-1),2,1),1,(length(n)-1)*2), 0];
    plot(v,n,'-','color',color,'linewidth',linewidth);
    try; plotfill(v,zeros(size(v)),n, 'color', color, 'transparency', 0.2); end
    
end
