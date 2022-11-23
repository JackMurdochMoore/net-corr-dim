% EST_CORR_DIM_4: estimate correlation dimension using four different
% methods. 
%
% Uses the model
% C(s) = s^D,
% c(s) = s^D - (s - 1)^D
% where D is correlation dimension and C(s) (c(s)) is correlation integral (correlation) at distance s.
% 
% REQUIRED INPUTS:
% ss - Vector of network distances
% nn - Number of node pairs at each network distance
% 
% OPTIONAL INPUTS:
% DLims - Upper and lower bound for correlation dimension, entered as the
%         two-vector [DLower, DUpper] 
% 
% OUTPUTS:
% DVec - Vector of estimated correlation dimensions
% sMaxVec - Vector of estimated upper cutoffs sMax
% D2Mat - D2Mat(i, j) records network correlation dimension D estimated
%         using method i with upper cutoff sMax chosen using method j
% codeCell - Code for method of estimation
% descriptionCell - Description of method of estimation
% DMat - DMat(i, j) records estimate of network correlation dimension D
%        made using method i using upper cutoff number j
% objFuncMat - Value at optimum of the function which is minimised (either
%              locally or globally) to choose the upper cutoff sMax 
% 
% Example use:
% D = 2; N0 = 1000; k = 10; p = 0;
% A = small_world_manhattan(N0, k, D, p); N = size(A, 1);
% [ss, nn] = count_distances(A);
% DLims = [-Inf, Inf]; [DVec, sMaxVec, D2Mat, codeCell, descriptionCell, DMat, objFuncMat] = est_corr_dim_4(ss, nn, DLims);
% 
% 
% Associated with 
% "Correlation dimension in empirical networks" 
% by 
% Jack Murdoch Moore, Haiying Wang, Michael Small, Gang Yan, Huijie Yang, 
% and Changgui Gu. 
%
function [DVec, sMaxVec, D2Mat, codeCell, descriptionCell, DMat, objFuncMat] = est_corr_dim_4(ss, nn, varargin)
if (numel(varargin) == 0)
    DLims = [-Inf, Inf];
else
    DLims = varargin{1};
end
DLower = DLims(1); DUpper = DLims(2);

codeCell = {'CE', 'KS', 'RC'};
descriptionCell = {...
    'Minimum negative log-likelihood per observation using model C(s) ∝ s^D',...
    'Minimum KS distance and maximum likelihood using model C(s) ∝ s^D',...
    'Linear regression to log-log plot of correlation sum (i.e., using model C(s) ∝ s^D)',...
    };
numDTypes = numel(codeCell);

num_s = numel(ss);

DMat = NaN(numDTypes, num_s);
objFuncMat = NaN(numDTypes, num_s);

options = optimoptions('fminunc', 'Display', 'off');
for ii_s = 2:num_s
    ss1 = ss(1:ii_s);
    nn1 = nn(1:ii_s);
    NN1 = cumsum(nn1);
    if (numel(ss1) >= 2)
        fun = @(g) -log_like_4(g, ss1, nn1);
        if sum(nn1) > 0
            [DML, ~] = fminunc(fun, 2, options);
        else
            DML = NaN;
        end
        DLimImposed = 0;
        if (DML < DLower); DML = DLower; DLimImposed = 1; end
        if (DML > DUpper); DML = DUpper; DLimImposed = 1; end
        negLogLike = fun(DML);
        scaledNegLogLike = negLogLike - log(1 + max(ss1) - min(ss1));
        if ~isequal(size(scaledNegLogLike), [1, 1]); scaledNegLogLike = NaN; end
        objFuncMat(1, ii_s) = scaledNegLogLike;%Objective function for CE
        
        DMat(1:2, ii_s) = DML;%Dimension estimated from CE and KS
        if (~DLimImposed) && (ii_s <= 2)
            KS1 = 0;
        else
            n1CumFit = @(a) a.^DML*(sum(nn1)/(max(ss1).^DML));
            nn1CumFit = n1CumFit(ss1);
            nn1Cum = cumsum(nn1);
            KS1 = max(abs(nn1Cum - nn1CumFit), [], 'includenan')/sum(nn1);
        end
        objFuncMat(2, ii_s) = KS1;%Objective function for KS
        
        linFit = polyfit(log(ss1), log(NN1), 1);
        DRC = linFit(1);
        if (DRC < DLower); DRC = DLower; end
        if (DRC > DUpper); DRC = DUpper; end
        DMat(3, ii_s) = DRC;
        
        objFuncMat(3, ii_s) = -nn(ii_s);
    end
end

DVec = NaN(numDTypes, 1); sMaxVec = NaN(numDTypes, 1);
D2Mat = NaN(numDTypes, numDTypes);%D2Mat(i, j) records D estimated using method i with sMax chosen using method j.
for iiDType = 1:numDTypes
    objFunc = objFuncMat(iiDType, :);
    try
    if (iiDType ~= 2)%Final global minimum:
        [~, ii_s] = min(fliplr(objFunc)); ii_s = numel(objFunc) - (ii_s - 1);
    else
        %Final local minimum:
        [~, locMinInd] = find_local_minima(objFunc);
        ii_s = locMinInd(end);
    end
    DVec(iiDType) = DMat(iiDType, ii_s); sMaxVec(iiDType) = ss(ii_s);
    D2Mat(:, iiDType) = DMat(:, ii_s);
    catch
    end
end
end