% Given adjacency matrix A, work out numbers nn of pairs of nodes at each
% non-zero, non-maximum distance aa. 
%
%
% Associated with 
%
% "Correlation dimension in empirical networks" 
% by 
% Jack Murdoch Moore, Haiying Wang, Michael Small, Gang Yan, Huijie Yang, 
% and Changgui Gu
% 
function [ss, nn] = count_distances(A)

G = graph(A);

dd = G.distances; dd = dd((dd > 0) & (dd < Inf)); dd = dd(:); minDist = min(dd); maxDist = max(dd);

ss = minDist:maxDist;
if isempty(ss)
    nn = zeros(size(ss));
else
    nn = histcounts(dd, (min(ss) - 0.5):(max(ss) + 0.5));
end

end