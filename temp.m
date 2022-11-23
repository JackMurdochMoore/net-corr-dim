clear;
N = 1000; k = 2; D = 1; p = 0.01; sigma = 0; kappa = 1; omega = 0; rewireFlag = 1; lowerAndUpperQuantile = [-eps, 1 + eps]; A = small_world_manhattan(N, k, D, p, sigma, kappa, omega, rewireFlag, lowerAndUpperQuantile); G = graph(A); figure; plot(G, 'Layout', 'Force3'); allDeg = G.degree; minDeg = min(allDeg); maxDeg = max(allDeg);
figure; histogram(allDeg, minDeg:maxDeg); xlabel('Degree'); ylabel('Count');
allDist = G.distances; allDist = allDist(:); minDist = min(allDist); maxDist = max(allDist);
figure; nn = histcounts(allDist, (minDist - 0.5):(maxDist + 0.5));
figure; loglog(minDist:maxDist, nn); xlabel('Distance'); ylabel('Count');
figure; plot(minDist:maxDist, nn); xlabel('Distance'); ylabel('Count');