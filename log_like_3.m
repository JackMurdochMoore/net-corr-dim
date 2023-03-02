% Log-likelihood per observation for model
% c(s) = s^(d-1),
% where c(s) is correlation at distance s.
% 
% d is exponent (d + 1 is correlation dimension)
% ss is vector of distances from s = 1 to s = s_max
% nn is vector of counts of distances
%
% 
% Associated with 
% "Correlation dimension in empirical networks" 
% by 
% Jack Murdoch Moore, Haiying Wang, Michael Small, Gang Yan, Huijie Yang, 
% and Changgui Gu. 
%
function logLikePerObs = log_like_3(d, ss, nn)
unnormedPMat = ss.^d;
PMat = unnormedPMat./sum(unnormedPMat);
logLikePerObs = log(PMat)*(nn')/sum(nn);
end
