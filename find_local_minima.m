% Function to find local minima, used by est_corr_dim_3, est_corr_dim_4,
% and show_fit_func
% 
% 
function [mm, locMinInd] = find_local_minima(xx)
xx = xx(:)';%Make sure xx is a row vector
nonNANInd = find(~isnan(xx));
xx2 = xx(nonNANInd);
if (numel(xx2) == 0)
    mm = []; locMinInd = [];
elseif (numel(xx2) == 1)
    locMinInd = 1;
    locMinInd = nonNANInd(locMinInd);
    mm = xx(locMinInd);
else
    locMinInd = find(islocalmin(xx2, 'FlatSelection', 'all'));
    if xx2(1) <= xx2(2); locMinInd = [1, locMinInd]; end
    if xx2(end) <= xx2(end - 1); locMinInd = [locMinInd, numel(locMinInd)]; end
    locMinInd = unique(locMinInd);
    locMinInd = nonNANInd(locMinInd);
    mm = xx(locMinInd);
end
end