% Fit power-law to interval [1, s_max].  
% 
% Consider number of node pairs at a distance.
% 
% compare_corr_dim_est_2 uses the model
% c(s) = s^(D-1)
% where D is correlation dimension and c(s) is correlation at distance s.
%
%
% Associated with 
% "Correlation dimension in empirical networks" 
% by 
% Jack Murdoch Moore, Haiying Wang, Michael Small, Gang Yan, Huijie Yang, 
% and Changgui Gu. 
%
function [d1, d3, d1P, d3P, d1RP, d3RP] = show_fit_func(ss, nn, saveStr, saveFigFlag, varargin)

if (numel(varargin) == 0)
    saveFigFolderStr = '';
else
    saveFigFolderStr = [varargin{1}, '\'];
end

if (numel(varargin) <= 1)
    titleStr = '';
else
    titleStr = varargin{2};
end

M = sum(nn);
num_a = numel(ss);
sR = max(ss) + 1;%Maximum distance

% % pdf:

d1Mat = NaN(1, num_a);
d3Mat = NaN(1, num_a);
scaledNegLogLike1Mat = NaN(1, num_a);
scaledNegLogLike3Mat = NaN(1, num_a);
negLogLike1Mat = NaN(1, num_a);
negLogLike3Mat = NaN(1, num_a);
ku1Mat = NaN(1, num_a);
ks1Mat = NaN(1, num_a);
rmse1Mat = NaN(1, num_a);

scaledEntropyEmp1Mat = NaN(1, num_a);
scaledEntropyFit1Mat = NaN(1, num_a);

% dMin = 0; DMin = dMin + 1;
dMin = -Inf; DMin = dMin + 1;

optionsUnc = optimoptions('fminunc', 'Display', 'off');
for ii_sML = 1:num_a
    ss1 = ss(1:ii_sML);
    nn1 = nn(1:ii_sML);
    fun = @(g) -log_like_3(g, ss1, nn1);
    [d1, ~] = fminunc(fun, 1, optionsUnc);
    if (d1 < dMin); d1 = dMin; end
    negLogLike1 = fun(d1);
    d1Mat(ii_sML) = d1;
    scaledNegLogLike1 = negLogLike1 - log(1 + max(ss1) - min(ss1));
    scaledNegLogLike1 = round(scaledNegLogLike1, 10);%Assume that differences less than 10^-10 are meaningless
    if ~isequal(size(scaledNegLogLike1), [1, 1]); scaledNegLogLike1 = NaN; negLogLike1 = NaN; end
    negLogLike1Mat(ii_sML) = negLogLike1;
    scaledNegLogLike1Mat(ii_sML) = scaledNegLogLike1;
    
    m1 = sum(nn1);
    pp1 = nn1./m1;
    scaledEntropyEmp1 = -pp1*log(pp1') - log(1 + max(ss1) - min(ss1));
    scaledEntropyEmp1Mat(ii_sML) = scaledEntropyEmp1;
    qq1 = ss1.^d1; qq1 = qq1/sum(qq1);
    scaledEntropyFit1 = -qq1*log(qq1') - log(1 + max(ss1) - min(ss1));
    scaledEntropyFit1Mat(ii_sML) = scaledEntropyFit1;
    
    m1 = sum(nn1);
    n1Fit = @(a) a.^d1*(sum(nn1)/sum(ss1.^d1));
    nn1Fit = n1Fit(ss1);
    nn1FitCum = cumsum(nn1Fit);
    nn1Cum = cumsum(nn1);
    ku1 = (max(nn1FitCum - nn1Cum, [], 'includenan') + max(nn1Cum - nn1FitCum, [], 'includenan'))/m1;
    ks1 = max(abs(nn1Cum - nn1FitCum), [], 'includenan')/m1;
    
    if (numel(ss1) < 2); continue; end
    
    ku1Mat(ii_sML) = ku1;
    ks1Mat(ii_sML) = ks1;
    
    pFit = polyfit(log(ss1), log(nn1/m1), 1);
    rmse1 = sqrt(mean((log(nn1/m1) - polyval(pFit, log(ss1))).^2));
    rmse1Mat(ii_sML) = rmse1;
end

[~, ii_sML] = min(fliplr(scaledNegLogLike1Mat)); ii_sML = numel(scaledNegLogLike1Mat) - (ii_sML - 1);
scaledNegLogLike1 = scaledNegLogLike1Mat(ii_sML);
d1 = d1Mat(ii_sML);
ss1 = ss(1:ii_sML);
nn1 = nn(1:ii_sML);
n1Fit = @(a) a.^d1*(sum(nn1)/sum(ss1.^d1));
nn1Fit = n1Fit(ss1);
nn1FitCum = cumsum(nn1Fit);

if isempty(scaledNegLogLike1)
    d1 = NaN;
end

[~, ii_sPR] = max(nn);%Peak of nn to consider for tail dimension

ii_sPL = ii_sML;%To use same sMax for Code rate and Regression methods, have this uncommented

ss1PL = ss(1:ii_sPL);

d1P = d1Mat(ii_sPL);

d1RP = polyfit(log(ss(1:ii_sPL)), log(nn(1:ii_sPL)), 1);

[~, locMinInd] = find_local_minima(ks1Mat);
ii_s_ks = locMinInd(end);
d1_ks = d1Mat(ii_s_ks);
ss1_ks = ss(1:ii_s_ks);
nn1_ks = nn(1:ii_s_ks);

n1Fit_ks = @(a) a.^d1_ks*(sum(nn1_ks)/sum(ss1_ks.^d1_ks));
nn1Fit_ks = n1Fit_ks(ss1_ks);
nn1FitCum_ks = cumsum(nn1Fit_ks);

for ii_sMR = 1:num_a
    ss3 = ss(ii_sMR:num_a);
    nn3 = nn(ii_sMR:num_a);
    fun = @(x) -log_like_3(x(1), sR - ss3, nn3);
    [d3, ~] = fminunc(fun, 1, optionsUnc);
    if (d3 < dMin); d3 = dMin; end
    negLogLike3 = fun(d3);
    d3Mat(ii_sMR) = d3;
    scaledNegLogLike3 = negLogLike3 - log(1 + max(ss3) - min(ss3));
    scaledNegLogLike3 = round(scaledNegLogLike3, 10);%Assume that differences less than 10^-10 are meaningless
    if ~isequal(size(scaledNegLogLike3), [1, 1]); scaledNegLogLike3 = NaN; negLogLike3 = NaN; end
    negLogLike3Mat(ii_sMR) = negLogLike3;
    scaledNegLogLike3Mat(ii_sMR) = scaledNegLogLike3;
end

[~, ii_sMR] = min(scaledNegLogLike3Mat);
scaledNegLogLike3 = scaledNegLogLike3Mat(ii_sMR);
d3 = d3Mat(ii_sMR);


if isempty(scaledNegLogLike3)
    d3 = NaN;
end

d3P = d3Mat(ii_sPR);


d3RP = polyfit(log(sR - ss(ii_sPR:num_a)), log(nn(ii_sPR:num_a)), 1); d3RP = d3RP(1);



% % CDF:

NN = cumsum(nn);

C1_NN = cumsum(nn);
C3_NN = fliplr(cumsum(fliplr(nn)));

D1MAT = NaN(1, num_a);
D3Mat = NaN(1, num_a);
SCALEDNEGLOGLIKE1MAT = NaN(1, num_a);
SCALEDNEGLOGLIKE3MAT = NaN(1, num_a);
NEGLOGLIKE1MAT = NaN(1, num_a);
NEGLOGLIKE3MAT = NaN(1, num_a);
KU1Mat = NaN(1, num_a);
KS1Mat = NaN(1, num_a);
RMSE1Mat = NaN(1, num_a);

for II_SM = 1:num_a
    SS1 = ss(1:II_SM);
    NN1 = C1_NN(1:II_SM);
    FUN = @(g) -log_like_3(g, SS1, NN1);
    [D1, ~] = fminunc(FUN, 1, optionsUnc);
    if (D1 < DMin); D1 = DMin; end
    NEGLOGLIKE1 = FUN(D1);
    D1MAT(II_SM) = D1;
    SCALEDNEGLOGLIKE1 = NEGLOGLIKE1 - log(1 + max(SS1) - min(SS1));
    SCALEDNEGLOGLIKE1 = round(SCALEDNEGLOGLIKE1, 10);%Assume that differences less than 10^-10 are meaningless
    if ~isequal(size(SCALEDNEGLOGLIKE1), [1, 1]); SCALEDNEGLOGLIKE1 = NaN; NEGLOGLIKE1 = NaN; end
    NEGLOGLIKE1MAT(II_SM) = NEGLOGLIKE1;
    SCALEDNEGLOGLIKE1MAT(II_SM) = SCALEDNEGLOGLIKE1;    
    
    AA3 = ss(II_SM:num_a);
    NN3 = C3_NN(II_SM:num_a);
    FUN = @(g) -log_like_3(g, sR - AA3, NN3);
    [D3, ~] = fminunc(FUN, 1, optionsUnc);
    if (D3 < DMin); D3 = DMin; end
    NEGLOGLIKE3 = FUN(D3);
    D3Mat(II_SM) = D3;
    SCALEDNEGLOGLIKE3 = NEGLOGLIKE3 - log(1 + max(AA3) - min(AA3));
    SCALEDNEGLOGLIKE3 = round(SCALEDNEGLOGLIKE3, 10);%Assume that differences less than 10^-10 are meaningless
    if ~isequal(size(SCALEDNEGLOGLIKE3), [1, 1]); SCALEDNEGLOGLIKE3 = NaN; NEGLOGLIKE3 = NaN; end
    NEGLOGLIKE3MAT(II_SM) = NEGLOGLIKE3;
    SCALEDNEGLOGLIKE3MAT(II_SM) = SCALEDNEGLOGLIKE3;
    
    M1 = sum(NN1);
    N1Fit = @(A) A.^D1*(sum(NN1)/sum(SS1.^D1));
    NN1Fit = N1Fit(SS1);
    NN1FitCum = cumsum(NN1Fit);
    NN1Cum = cumsum(NN1);
    KU1 = (max(NN1FitCum - NN1Cum, [], 'includenan') + max(NN1Cum - NN1FitCum, [], 'includenan'))/M1;
    KS1 = max(abs(NN1Cum - NN1FitCum), [], 'includenan')/M1;
    
    if (numel(SS1) < 2); continue; end
    
    KU1Mat(II_SM) = KU1;
    KS1Mat(II_SM) = KS1;
    
    PFit = polyfit(log(SS1), log(NN1/M1), 1);
    RMSE1 = sqrt(mean((log(NN1/M1) - polyval(PFit, log(SS1))).^2));
    RMSE1Mat(II_SM) = RMSE1;
end

[~, II_SML] = min(fliplr(SCALEDNEGLOGLIKE1MAT)); II_SML = numel(SCALEDNEGLOGLIKE1MAT) - (II_SML - 1);

D1RP = polyfit(log(ss(1:ii_sPL)), log(C1_NN(1:ii_sPL)), 1); C1RP = D1RP(2); D1RP = D1RP(1);
N1PRFIT = @(a) exp(C1RP)*a.^D1RP;%Infer multiplicative constant from least squares fit
NN1PRFIT = N1PRFIT(ss1PL);

fontSize = 22; lineWidth = 2; markerSize = 8;

load('colour_scheme.mat', 'colOrder');

trueCol = colOrder(8, :);%Grey
standardCol = colOrder(1, :);%Blue
proposedCol = colOrder(5, :);%Purple
ksCol = colOrder(4, :);%Red
normalisedCol = proposedCol;%
% unnormalisedCol = colOrder(10, :);%Light blue
unnormalisedCol = colOrder(3, :);%Green

legendBoxAlpha = 0.5;

disp(['DNLL', char(9), 'sMax', char(9), 'C(sMax)', ':'])
disp([num2str(d1 + 1), char(9), num2str(ss(ii_sML)), char(9), num2str(NN(ii_sML)/NN(end))]);

disp(['DKSD', char(9), 'sMax', char(9), 'C(sMax)', ':'])
disp([num2str(d1_ks + 1), char(9), num2str(ss(ii_s_ks)), char(9), num2str(NN(ii_s_ks)/NN(end))]);

% % Figures:
legCell = {'\begin{tabular}{l}True\end{tabular}', ['\begin{tabular}{l}Lin. regression\\($\hat{D} = ', sprintf('%.2f', D1RP), '$)\end{tabular}'], ['\begin{tabular}{l}Neg. log-like.\\($\hat{D} = ', sprintf('%.2f', d1 + 1), '$)\end{tabular}']};
legCell2 = {'True', 'Lin. regression', 'Neg. log-like.'};


% % CDF, min. neg. log. likelihood  (like Fig. 1, S1):

f = figure;
p = plot(ss, C1_NN/M, '-', ss1PL, NN1PRFIT/M, '--h', ss1, nn1FitCum/M, ':s', 'LineWidth', lineWidth, 'MarkerSize', markerSize);
leg = legend(legCell, 'FontSize', fontSize, 'Interpreter', 'LaTeX', 'Location', 'NorthWest', 'Color', 'None'); leg.BoxFace.ColorType='truecoloralpha'; leg.BoxFace.ColorData=uint8(255*[1; 1; 1; legendBoxAlpha]);
set(gca, 'XMinorTick', 'off', 'TickLabelInterpreter', 'LaTeX', 'FontSize', fontSize); xlabel('Distance, $s$', 'FontSize', fontSize, 'Interpreter', 'LaTeX'); ylabel('Fraction within dist., $C(s)$', 'FontSize', fontSize, 'Interpreter', 'LaTeX');
titleStr2 = titleStr;
if isequal(titleStr2(19:29), 'Small world')
    titleStr2 = 'Small world';
else
    titleStr2 = titleStr2(19:(strfind(titleStr, '\\') - 1));
end
title(titleStr2, 'Interpreter', 'LaTeX', 'FontSize', fontSize);
p(1).Color = trueCol; p(1).MarkerFaceColor = trueCol;
p(2).Color = standardCol; p(2).MarkerFaceColor = standardCol;
p(3).Color = proposedCol; p(3).MarkerFaceColor = 'w';
xlim([1, max(ss1) + 1]);

if saveFigFlag
    exportgraphics(f, [saveFigFolderStr, 'corr-int_lin-scale_', saveStr, '.png'], 'BackgroundColor', 'None', 'ContentType', 'vector');
    exportgraphics(f, [saveFigFolderStr, 'corr-int_lin-scale_', saveStr, '.pdf'], 'BackgroundColor', 'None', 'ContentType', 'vector');
    savefig(f, [saveFigFolderStr, 'corr-int_lin-scale_', saveStr, '.fig']);
end

f = figure;
p = loglog(ss, C1_NN/M, '-', ss1PL, NN1PRFIT/M, '--h', ss1, nn1FitCum/M, ':s', 'LineWidth', lineWidth, 'MarkerSize', markerSize);
set(gca, 'XMinorTick', 'off', 'TickLabelInterpreter', 'LaTeX', 'FontSize', fontSize); xlabel('Distance, $s$', 'FontSize', fontSize, 'Interpreter', 'LaTeX'); ylabel('Fraction within dist., $C(s)$', 'FontSize', fontSize, 'Interpreter', 'LaTeX');
title(titleStr, 'Interpreter', 'LaTeX', 'FontSize', fontSize);
leg = legend(legCell2, 'FontSize', fontSize, 'Interpreter', 'LaTeX', 'Location', 'SouthEast', 'Color', 'None'); leg.BoxFace.ColorType='truecoloralpha'; leg.BoxFace.ColorData=uint8(255*[1; 1; 1; legendBoxAlpha]);
p(1).Color = trueCol; p(1).MarkerFaceColor = trueCol;
p(2).Color = standardCol; p(2).MarkerFaceColor = standardCol;
p(3).Color = proposedCol; p(3).MarkerFaceColor = 'w';
xticks([1, 2, 5, 10, 20, 50, 100, 200, 500, 1000]);
yLim = ylim; ylim([yLim(1), 1]);
yticks(10.^(-10:0));
xlim([1, max(ss)]);

if saveFigFlag
    exportgraphics(f, [saveFigFolderStr, 'corr-int_log-scale_', saveStr, '.png'], 'BackgroundColor', 'None', 'ContentType', 'vector');
    exportgraphics(f, [saveFigFolderStr, 'corr-int_log-scale_', saveStr, '.pdf'], 'BackgroundColor', 'None', 'ContentType', 'vector');
    savefig(f, [saveFigFolderStr, 'corr-int_log-scale_', saveStr, '.fig']);
end

% % CDF, min. KS distance and min. neg. log-like.  (like Fig. 7):

f = figure; hold on;
% if (max(ss1) < 100)
    p = plot(ss, C1_NN/M, '-', ss1_ks, nn1FitCum_ks/M, ':d', ss1, nn1FitCum/M, ':s', 'LineWidth', lineWidth, 'MarkerSize', markerSize);
    p(1).Color = trueCol; p(1).MarkerFaceColor = trueCol;
    p(2).Color = ksCol; p(2).MarkerFaceColor = 'w';
    p(3).Color = proposedCol; p(3).MarkerFaceColor = 'w';
% else
%     p = plot(ss, C1_NN/M, '-', ss1_ks, nn1FitCum_ks/M, ':', ss1, nn1FitCum/M, ':', 'LineWidth', lineWidth, 'MarkerSize', markerSize);
%     p(1).Color = trueCol; p(1).MarkerFaceColor = trueCol;
%     p(2).Color = ksCol; p(2).MarkerFaceColor = 'w';
%     p(3).Color = proposedCol; p(3).MarkerFaceColor = 'w';
%     delta_s = 50;
%     p = plot(ss1_ks(1:delta_s:max(ss1_ks)), nn1FitCum_ks(1:delta_s:max(ss1_ks))/M, 'd', ss1(1:delta_s:max(ss1)), nn1FitCum(1:delta_s:max(ss1))/M, 's', 'LineWidth', lineWidth, 'MarkerSize', markerSize);
%     p(1).Color = ksCol; p(1).MarkerFaceColor = 'w';
%     p(2).Color = proposedCol; p(2).MarkerFaceColor = 'w';
% end
leg = legend({'True', ['KS distance ($\hat{D} = ', sprintf('%.2f', d1_ks + 1), '$)'], ['Neg. log-like. ($\hat{D} = ', sprintf('%.2f', d1 + 1), '$)']}, 'FontSize', fontSize, 'Interpreter', 'LaTeX', 'Location', 'NorthWest', 'Color', 'None'); leg.BoxFace.ColorType='truecoloralpha'; leg.BoxFace.ColorData=uint8(255*[1; 1; 1; legendBoxAlpha]);
set(gca, 'XMinorTick', 'off', 'TickLabelInterpreter', 'LaTeX', 'FontSize', fontSize); xlabel('Distance, $s$', 'FontSize', fontSize, 'Interpreter', 'LaTeX'); ylabel('Fraction within dist., $C(s)$', 'FontSize', fontSize, 'Interpreter', 'LaTeX');
title(titleStr, 'Interpreter', 'LaTeX', 'FontSize', fontSize);
xlim([1, max([ss1_ks, ss1]) + 1]);
if isequal(titleStr, '\begin{tabular}{c}Small world from 1-dim. lattice\\$\langle k \rangle=2, p=0.01, N=9612$\end{tabular}')
    yLim = ylim; ylim([0, 1.6*yLim(2)]);
elseif isequal(titleStr, '\begin{tabular}{c}Visual cortex: Mouse\\$\langle k \rangle=3.11, N=987$\end{tabular}')
    yLim = ylim; ylim([0, 1.8*yLim(2)]);
else
    yLim = ylim; ylim([0, 1.4*yLim(2)]);
end

if saveFigFlag
    exportgraphics(f, [saveFigFolderStr, 'corr-int_lin-scale_nll_ksd_', saveStr, '.png'], 'BackgroundColor', 'None', 'ContentType', 'vector');
    exportgraphics(f, [saveFigFolderStr, 'corr-int_lin-scale_nll_ksd_', saveStr, '.pdf'], 'BackgroundColor', 'None', 'ContentType', 'vector');
    savefig(f, [saveFigFolderStr, 'corr-int_lin-scale_nll_ksd_', saveStr, '.fig']);
end



% % Fitting process, min. neg. log. likelihood, fit to correlation c(s)  (like Fig. 2b):

f = figure; hold on;
set(gca, 'FontSize', fontSize);
yyaxis left;
ax = gca;
ax.YAxis(1).Color = unnormalisedCol;
ax.YAxis(2).Color = normalisedCol;
yyaxis left;
plot(ss, scaledNegLogLike1Mat + log(ss), '--', 'LineWidth', lineWidth, 'MarkerSize', markerSize, 'Color', unnormalisedCol);
ylabel('Without norm.', 'FontSize', fontSize, 'Interpreter', 'LaTeX');
yyaxis right;
plot(ss, scaledNegLogLike1Mat, ':', 'LineWidth', lineWidth, 'MarkerSize', markerSize, 'Color', normalisedCol);
ylabel('With normalisation', 'FontSize', fontSize, 'Interpreter', 'LaTeX');
yLim = ylim;
plot(ss(ii_sML)*[1, 1], yLim, ':', 'LineWidth', lineWidth, 'MarkerSize', markerSize, 'Color', proposedCol);
set(gca, 'XMinorTick', 'off', 'TickLabelInterpreter', 'LaTeX');
xlabel('Maximum of interval, $s_{\max}$', 'FontSize', fontSize, 'Interpreter', 'LaTeX');
ylim(yLim); xlim([1, max(ss)]);
xlim([1, max(ss)]);
leg = legend({'Without norm.', 'With norm.'}, 'FontSize', fontSize, 'Interpreter', 'LaTeX', 'Location', 'North', 'Color', 'None');
leg.BoxFace.ColorType='truecoloralpha'; leg.BoxFace.ColorData=uint8(255*[1; 1; 1; legendBoxAlpha]);
title('Negative log-likelihood per observation', 'FontSize', fontSize, 'Interpreter', 'LaTeX');

ax.Position = [0.1524, 0.2216, 0.7364, 0.5867];

box on;

ax = gca;
ax.Position(4) = ax.Position(4) - 0.1*ax.Position(2);

if saveFigFlag
    exportgraphics(f, [saveFigFolderStr, 'likelihood_', saveStr, '.pdf'], 'BackgroundColor', 'None', 'ContentType', 'vector');
    savefig(f, [saveFigFolderStr, 'likelihood_', saveStr, '.fig']);
end

end