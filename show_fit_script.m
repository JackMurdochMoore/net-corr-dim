close all;
clear;

rng(0);

sigma = 0; kappa = 1; omega = 0; rewireFlag = 1; lowerAndUpperQuantile = [-eps, 1 + eps];

% D00 = 1; N00 = 10000; k00 = 2*D00; p00 = 0;
% D00 = 1; N00 = 10000; k00 = 2*D00; p00 = 0.001;
% D00 = 1; N00 = 10000; k00 = 2*D00; p00 = 0.01;
% D00 = 1; N00 = 1000; k00 = 10; p00 = 0;
% D00 = 1; N00 = 5000; k00 = 10; p00 = 0;
% D00 = 1; N00 = 10000; k00 = 10; p00 = 0;

% D00 = 2; N00 = 1000; k00 = 2*D00; p00 = 0.00;
% D00 = 2; N00 = 10000; k00 = 2*D00; p00 = 0.00;
% D00 = 2; N00 = 10000; k00 = 2*D00 + 2; p00 = 0.00;
% D00 = 2; N00 = 1000; k00 = 10; p00 = 0.00;
% D00 = 2; N00 = 10000; k00 = 10; p00 = 0.00;
% D00 = 2; N00 = 10000; k00 = 10; p00 = 0.00;
% D00 = 2; N00 = 10000; k00 = 10; p00 = 0.01;

% D00 = 3; N00 = 1000; k00 = 2*D00; p00 = 0.00;
% D00 = 3; N00 = 10000; k00 = 2*D00; p00 = 0.00;
% D00 = 3; N00 = 1000; k00 = 10; p00 = 0.00;
% D00 = 3; N00 = 10000; k00 = 10; p00 = 0.01;
% D00 = 3; N00 = 10000; k00 = 10; p00 = 0.01;

D00 = 3;
k00 = 2*D00;
N00 = 10000;
p00 = 0.01;

N0 = N00; D = D00; k = k00; p = p00;

% networkFlag = -3;%Largest connected component of small world

% networkFlag = 0;%Small world

% networkFlag = 1;% TV show% N = 3,892 
% networkFlag = 2;% Power grid% N = 4,941 
% networkFlag = 3;% Politician% N = 5,908 
% networkFlag = 12;% Computer science PhD % N = 1,025 
% networkFlag = 13;% Erdos collaboration network % N = 4,991 
% networkFlag = 16;% School friendship network # 27 % N = 1,152 
% networkFlag = 17;% School friendship network # 67 % N = 439 
% networkFlag = 18;% Mouse visual cortex % N = 987 
% networkFlag = 19;% Yeast protein interactions % N = 1,458 
% networkFlag = 20;% WormNet DM-LC% N = 483 
% networkFlag = 21;% WormNet DM-HT% N = 2,831 
% networkFlag = 22;% Human disease % N = 516 
% networkFlag = 23;% WormNet CE-LC% N = 993 
% networkFlag = 24;% WormNet CE-HT% N = 2,194 
% networkFlag = 25;% EPA% N = 4,772 
% networkFlag = 26;% road-minnesota.mtx% N = 2,640 
% networkFlag = 27;% road-euroroad.edges% N = 1,039 
% networkFlag = 28;% s208_st.txt% N = 122 
% networkFlag = 29;% s420_st.txt% N = 252 
networkFlag = 30;% s838_st.txt% N = 512 

if (networkFlag == 0)
    A = small_world_manhattan(N0, k, D, p, sigma, kappa, omega, rewireFlag, lowerAndUpperQuantile);
    nameStr = ['small-world-manhattan', '_p-', num2str(p), '_D-', num2str(D)];
    N = size(A, 1); k = sum(A(:))/N;
    if (p > 0)
        titleStr = ['\begin{tabular}{c}Small world from ', num2str(D), '-dim. lattice\\$\langle k \rangle=', num2str(k), ', p=', num2str(p), ', N=', num2str(N), '$\end{tabular}'];
    else
        titleStr = ['\begin{tabular}{c}', num2str(D), '-dimensional lattice\\', '$k=', num2str(k), ', N=', num2str(N), '$\end{tabular}'];
    end
elseif (networkFlag > 0)
    [A, nameStr, titleStr, loadStr] = load_network(networkFlag);
    N = size(A, 1); k = sum(A(:))/N;
    titleStr = ['\begin{tabular}{c}', titleStr, '\\$\langle k \rangle=', num2str(k, 3), ', N=', num2str(N), '$\end{tabular}'];
elseif (networkFlag == -3)
    A = small_world_manhattan_lcc(N0, k, D, p, sigma, kappa, omega, rewireFlag, lowerAndUpperQuantile);
    nameStr = ['small-world-manhattan-lcc', '_p-', num2str(p), '_D-', num2str(D)];
    N = size(A, 1); k = sum(A(:))/N;
    if (p > 0)
        titleStr = ['\begin{tabular}{c}Small world from ', num2str(D), '-dim. lattice\\$\langle k \rangle=', num2str(k), ', p=', num2str(p), ', N=', num2str(N), '$\end{tabular}'];
    else
        titleStr = ['\begin{tabular}{c}', num2str(D), '-dimensional lattice\\', '$k=', num2str(k), ', N=', num2str(N), '$\end{tabular}'];
    end
end
saveStr = [nameStr, '_N-', num2str(N), '_k-', num2str(k, 3)];

[aa, nn] = count_distances(A);

saveFigFolder = 'figures';%Where to save figures

saveFigFlag = 0;%Whether to save figures

disp(nameStr);

% close all;
[d1, d3, d1P, d3P, d1RP, d3RP] = show_fit_func(aa, nn, saveStr, saveFigFlag, saveFigFolder, titleStr);

% return;

numTriangles = trace(A^3)/6; B = A^2; numConnectedTriples = (sum(B(:)) - trace(B))/2; phi = numTriangles/numConnectedTriples;%Global clustering coefficient/transitivity
G = graph(A);
distMat0 = G.distances; distMat = distMat0(~eye(N)); dMaxFinite = max(distMat(isfinite(distMat)));
E = mean(distMat.^(-1)); L = mean(distMat); %'E' is network efficiency. 'NEffective' is a new measure called the effective network size. The effective network size NEffective is the mean number of communications a node can receive per time step
[~, C] = clusteringcoef(G);%Mean local clustering coefficient

disp(['N', char(9), 'k', char(9), '<s>', char(9), 'E', char(9), 'diam', char(9), 'phi', ':'])
disp([num2str(N), char(9), num2str(k), char(9), num2str(L), char(9), num2str(E), char(9), num2str(dMaxFinite), char(9), num2str(phi)]);