% Log-likelihood per observation for model
% C(s) = s^D,
% c(s) = s^D - (s - s)^D
% where C(s) (c(s)) is correlation integral (correlation) at distance s.
% 
% D is exponent/correlation dimension
% ss is vector of distances from s = 1 to s = s_max
% nn is vector of counts of distances
%
function logLikePerObs = log_like_4(D, ss, nn)
unnormedPMat = ss.^D - (ss - 1).^D;
PMat = unnormedPMat./sum(unnormedPMat);
logLikePerObs = log(PMat)*(nn')/sum(nn);
end