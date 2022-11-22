% Estimate correlation dimension of synthetic networks using different
% methods and save results in folder
% results-est-dim-3
% 
% compare_corr_dim_est_3 uses the model
% c(s) = s^(D-1),
% where D is correlation dimension and c(s) is correlation at distance s.
%
% Example use:
% D00 = 3; N00 = 10000; k00 = 8; numRep = 100; compare_corr_dim_est_3(D00, N00, k00, numRep); 
% 
function compare_corr_dim_est_3(D00, N00, k00, numRep)

saveDataFolder = 'results-est-dim-3';

rng(0);

maxFailures = numRep;

sigma = 0; kappa = 1; omega = 0; rewireFlag = 1; lowerAndUpperQuantile = [-eps, 1 + eps];

DLims = [-Inf, Inf];
[~, ~, ~, codeCell, descriptionCell, ~, ~] = est_corr_dim_3(1, 2, DLims);
numMeth = numel(codeCell);

p00 = 0;

% versus dimension:

N0 = N00; D = NaN; k = k00; p = p00;
DD = 1:5; numD = numel(DD);
saveStr = ['est-dim_small-world-manhattan_DLim-', num2str(DLims(1)), '-', num2str(DLims(2)), '_', num2str(numRep), 'rep', '_k-', num2str(k), '_', num2str(numD), 'D-', num2str(min(DD)), '-', num2str(max(DD)), '_p-', num2str(p), '_N0-', num2str(N0)];
if ~exist([pwd, '\', saveDataFolder, '\', saveStr, '.mat'], 'file')
    try
    numFailures = 0;
    dimMatDD = NaN(numMeth, numD, numRep);
    cutoffMatDD = NaN(numMeth, numD, numRep);
    D2MatCellDD = cell(numD, numRep);
    for iiRep = 1:numRep
        for iiD = 1:numD
            D = DD(iiD);
            A = [];
            while (numFailures <= maxFailures) && isempty(A)
                try A = small_world_manhattan(N0, k, D, p, sigma, kappa, omega, rewireFlag, lowerAndUpperQuantile); N = size(A, 1);
                catch
                    numFailures = numFailures + 1;
                end
            end
            assert(numFailures <= maxFailures, 'Too many unsuccessful attempts to generate network.');
            [ss, nn] = count_distances(A);
            [DVec, sMaxVec, D2Mat, codeCell, descriptionCell, ~, ~] = est_corr_dim_3(ss, nn, DLims);
            dimMatDD(:, iiD, iiRep) = DVec;
            cutoffMatDD(:, iiD, iiRep) = sMaxVec;
            D2MatCellDD{iiD, iiRep} = D2Mat;
        end
    end

    saveCell = {'numRep', 'N0', 'DD', 'k', 'p', 'dimMatDD', 'cutoffMatDD', 'D2MatCellDD', 'codeCell', 'descriptionCell', 'saveStr', 'saveCell'};
    save([saveDataFolder, '\', saveStr, '.mat'], saveCell{:});
    disp(['Saved ', saveStr, '.']);
    catch ME
        disp(['Error generating ', saveStr, '.']);
        disp(ME.message);
    end
else
    disp(['Skipping ', saveStr, '.']);
end


% versus degree:

N0 = N00; D = D00; k = NaN; p = p00;
kk = 2:2:20; num_k = numel(kk);
saveStr = ['est-dim_small-world-manhattan_DLim-', num2str(DLims(1)), '-', num2str(DLims(2)), '_', num2str(numRep), 'rep', '_', num2str(num_k), 'k-', num2str(min(kk)), '-', num2str(max(kk)), '_D-', num2str(D), '_p-', num2str(p), '_N0-', num2str(N0)];
if ~exist([pwd, '\', saveDataFolder, '\', saveStr, '.mat'], 'file')
    try
    numFailures = 0;
    dimMat_kk = NaN(numMeth, num_k, numRep);
    cutoffMat_kk = NaN(numMeth, num_k, numRep);
    D2MatCell_kk = cell(num_k, numRep);
    for iiRep = 1:numRep
        for ii_k = 1:num_k
            k = kk(ii_k);
            A = [];
            while (numFailures <= maxFailures) && isempty(A)
                try A = small_world_manhattan(N0, k, D, p, sigma, kappa, omega, rewireFlag, lowerAndUpperQuantile); N = size(A, 1);
                catch
                    numFailures = numFailures + 1;
                end
            end
            assert(numFailures <= maxFailures, 'Too many unsuccessful attempts to generate network.');
            [ss, nn] = count_distances(A);
            [DVec, sMaxVec, D2Mat, codeCell, descriptionCell, ~, ~] = est_corr_dim_3(ss, nn, DLims);
            dimMat_kk(:, ii_k, iiRep) = DVec;
            cutoffMat_kk(:, ii_k, iiRep) = sMaxVec;
            D2MatCell_kk{ii_k, iiRep} = D2Mat;
        end
    end

    saveCell = {'numRep', 'N0', 'D', 'kk', 'p', 'dimMat_kk', 'cutoffMat_kk', 'D2MatCell_kk', 'codeCell', 'descriptionCell', 'saveStr', 'saveCell'};
    save([saveDataFolder, '\', saveStr, '.mat'], saveCell{:});
    disp(['Saved ', saveStr, '.']);
    catch ME
        disp(['Error generating ', saveStr, '.']);
        disp(ME.message);
    end
else
    disp(['Skipping ', saveStr, '.']);
end

%versus network size:

N0 = NaN; D = D00; k = k00; p = p00;
NN0 = [1, 2, 5]'*(10.^(2:3)); NN0 = NN0(:)'; NN0 = [NN0, 10^4]; NN0 = unique(NN0); numN = numel(NN0);
saveStr = ['est-dim_small-world-manhattan_DLim-', num2str(DLims(1)), '-', num2str(DLims(2)), '_', num2str(numRep), 'rep', '_k-', num2str(k), '_D-', num2str(D), '_p-', num2str(p), '_', num2str(numN), 'N0-', num2str(min(NN0)), '-', num2str(max(NN0))];
if ~exist([saveDataFolder, '\', saveStr, '.mat'], 'file')
    try
    numFailures = 0;
    NN = NaN(1, numN);
    dimMatNN = NaN(numMeth, numN, numRep);
    cutoffMatNN = NaN(numMeth, numN, numRep);
    D2MatCellNN = cell(numN, numRep);
    for iiRep = 1:numRep
        for iiN = 1:numN
            N0 = NN0(iiN);
            A = [];
            while (numFailures <= maxFailures) && isempty(A)
                try A = small_world_manhattan(N0, k, D, p, sigma, kappa, omega, rewireFlag, lowerAndUpperQuantile); N = size(A, 1);
                catch
                    numFailures = numFailures + 1;
                end
            end
            assert(numFailures <= maxFailures, 'Too many unsuccessful attempts to generate network.');
            N = size(A, 1); NN(iiN) = N;
            [ss, nn] = count_distances(A);
            [DVec, sMaxVec, D2Mat, codeCell, descriptionCell, ~, ~] = est_corr_dim_3(ss, nn, DLims);
            dimMatNN(:, iiN, iiRep) = DVec;
            cutoffMatNN(:, iiN, iiRep) = sMaxVec;
            D2MatCellNN{iiN, iiRep} = D2Mat;
        end
    end

    saveCell = {'numRep', 'NN0', 'D', 'k', 'p', 'NN', 'dimMatNN', 'cutoffMatNN', 'D2MatCellNN', 'codeCell', 'descriptionCell', 'saveStr', 'saveCell'};
    save([saveDataFolder, '\', saveStr, '.mat'], saveCell{:});
    disp(['Saved ', saveStr, '.']);
    catch ME
        disp(['Error generating ', saveStr, '.']);
        disp(ME.message);
    end
else
    disp(['Skipping ', saveStr, '.']);
end

end