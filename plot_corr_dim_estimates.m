% Plot previously saved results of benchmarking different approaches to
% estimating correlation dimension.
%
% Associated with 
% "Correlation dimension in empirical networks" 
% by 
% Jack Murdoch Moore, Haiying Wang, Michael Small, Gang Yan, Huijie Yang, 
% and Changgui Gu. 
%
close all;
clear;

saveFigFlag = 1;

loadDataFolderCell = {'results-est-dim-4', 'results-est-dim-3'}; numLoadDataFolder = numel(loadDataFolderCell);
saveFigFolder = 'figures';

load('colour_scheme.mat', 'colOrder');

colMat = colOrder([8, 5, 4, 1, 3], :);
styleCellCell = {{'-', '--s', '--d', '--h', '--<'}, {'-', ':s', ':d', ':h', ':<'}};
nameCell = {'True', 'Code rate', 'KS distance', 'Lin. regression'};

toPlot = [1, 4, 3, 2]; numToPlot = numel(toPlot);

fontSize = 22; lineWidth = 2; markerSize = 8;
legendBoxAlpha = 1;

numRep = 100; D00 = 4; N00 = 10^4; k00 = 2*D00 + 2; p00 = 0;

DLims = [-Inf, Inf];

N0 = N00;

% Dimension (like Fig. 5):

try

f = figure; hold on;

for iiLoadDataFolder = 1:numLoadDataFolder

loadDataFolder = loadDataFolderCell{iiLoadDataFolder};
styleCell = styleCellCell{iiLoadDataFolder};

switch iiLoadDataFolder
    case 2
        loadStr = [loadDataFolder, '\', 'est-dim_small-world-manhattan_DLim-', num2str(DLims(1)), '-', num2str(DLims(2)), '_', num2str(numRep), 'rep', '_k-', num2str(k00), '_5D-1-5', '_p-', num2str(p00), '_N0-', num2str(N00)];
    case 1
        loadStr = [loadDataFolder, '\', 'est-dim-2_small-world-manhattan_DLim-', num2str(DLims(1)), '-', num2str(DLims(2)), '_', num2str(numRep), 'rep', '_k-', num2str(k00), '_5D-1-5', '_p-', num2str(p00), '_N0-', num2str(N00)];
end

load([loadStr, '.mat']);

numMeth = numel(codeCell);
xx = DD;
Y = mean(dimMatDD, 3);%./repmat(DD, [numMeth, 1]);
YStd = std(dimMatDD, 0, 3);%./repmat(DD, [numMeth, 1]);
YErr = YStd/sqrt(numRep - 1);
for iiPlot = 1:numToPlot
    iiMeth = toPlot(iiPlot);
    if (iiMeth == 1)
        if (iiLoadDataFolder == 2); continue; end
        yy = DD;%./DD;
        yyErr = zeros(size(DD));
        yyStd = zeros(size(DD));
    else
        yy = Y(iiMeth - 1, :);
        yyErr = YErr(iiMeth - 1, :);
        yyStd = YStd(iiMeth - 1, :);
    end
    col = colMat(iiMeth, :);
    sty = styleCell{iiMeth};
    if (iiMeth ~= 1) || (iiLoadDataFolder ~= 2)
        p = plot(xx, yy, sty, 'Color', col, 'LineWidth', lineWidth, 'MarkerSize', markerSize);
        if (iiLoadDataFolder == 1)
            p.MarkerFaceColor = col;
        else
            p.MarkerFaceColor = 'w';
        end
    end
end

end

set(gca, 'XTick', DD, 'XMinorTick', 'off', 'TickLabelInterpreter', 'LaTeX', 'FontSize', fontSize);

axis('tight');
daspect([1, 1, 1]);
ylabel('Estimated dimension, $\hat{D}$', 'Interpreter', 'LaTeX', 'FontSize', fontSize);
xlabel('Dimension, $D$', 'Interpreter', 'LaTeX', 'FontSize', fontSize);
if p > 0
    title(['$\langle k \rangle = ', num2str(k00), '$'], 'Interpreter', 'LaTeX', 'FontSize', fontSize);
else
    title(['$k = ', num2str(k00), '$'], 'Interpreter', 'LaTeX', 'FontSize', fontSize);
end
box on;

xLim = xlim; xRange = diff(xLim); xLim = [1, xRange + 1]; xlim(xLim);

nameStr = [saveFigFolder, '\', 'est-dim-1-2_small-world-manhattan_DLim-', num2str(DLims(1)), '-', num2str(DLims(2)), '_N0-', num2str(N00), '_p-', num2str(p00), '_k-', num2str(k00)];

if saveFigFlag
% savefig([nameStr, '.fig']);
% exportgraphics(f, [nameStr, '.png'], 'BackgroundColor', 'None', 'ContentType', 'vector');
exportgraphics(f, [nameStr, '.pdf'], 'BackgroundColor', 'None', 'ContentType', 'vector');
end

catch
    disp('Cannot plot estimated dimension versus true dimension.');
end

% return;

% Mean degree (like Fig. 4):

try

f = figure; hold on;

for iiLoadDataFolder = 1:numLoadDataFolder

loadDataFolder = loadDataFolderCell{iiLoadDataFolder};
styleCell = styleCellCell{iiLoadDataFolder};

switch iiLoadDataFolder
    case 2
        loadStr = [loadDataFolder, '\', 'est-dim_small-world-manhattan_DLim-', num2str(DLims(1)), '-', num2str(DLims(2)), '_', num2str(numRep), 'rep', '_10k-2-20', '_D-', num2str(D00), '_p-', num2str(p00), '_N0-', num2str(N0)];
    case 1
        loadStr = [loadDataFolder, '\', 'est-dim-2_small-world-manhattan_DLim-', num2str(DLims(1)), '-', num2str(DLims(2)), '_', num2str(numRep), 'rep', '_10k-2-20', '_D-', num2str(D00), '_p-', num2str(p00), '_N0-', num2str(N0)];
end
    

load([loadStr, '.mat']);

xx = kk;
set(gca, 'XMinorTick', 'off', 'TickLabelInterpreter', 'LaTeX', 'FontSize', fontSize);
Y = mean(dimMat_kk, 3);
YStd = std(dimMat_kk, 0, 3);
YErr = YStd/sqrt(numRep - 1);
for iiPlot = 1:numToPlot
    iiMeth = toPlot(iiPlot);
    if (iiMeth == 1)
        if (iiLoadDataFolder == 2); continue; end
        yy = D*ones(size(xx));
        yyErr = zeros(size(xx));
        yyStd = zeros(size(xx));
    else
        yy = Y(iiMeth - 1, :);
        yyErr = YErr(iiMeth - 1, :);
        yyStd = YStd(iiMeth - 1, :);
    end
    col = colMat(iiMeth, :);
    sty = styleCell{iiMeth};
    p = plot(xx, yy, sty, 'Color', col, 'LineWidth', lineWidth, 'MarkerSize', markerSize);
    if (iiLoadDataFolder == 1)
        p.MarkerFaceColor = col;
    else
        p.MarkerFaceColor = 'w';
    end
end

end

xlim([2*D, max(xx)]);
yLim = ylim; yMid = D; yHalfRange = max(abs(yLim - yMid));
yHalfRange = max(yHalfRange, 0.2);
yLim = [yMid - yHalfRange, yMid + yHalfRange]; ylim(yLim);
if (diff(yLim) < 0.4)
    yTicks = yticks; yTicks = min(yTicks):0.1:max(yTicks); yticks(yTicks);
end
ylabel('Estimated dimension, $\hat{D}$', 'Interpreter', 'LaTeX', 'FontSize', fontSize);
if p > 0
    xlabel('Mean degree, $\langle k \rangle$', 'Interpreter', 'LaTeX', 'FontSize', fontSize);
else
    xlabel('Degree, $k$', 'Interpreter', 'LaTeX', 'FontSize', fontSize);
end
N = round(N0^(1/D00))^D00;
title(['$D = ', num2str(D00), ', N = ', num2str(N), '$'], 'Interpreter', 'LaTeX', 'FontSize', fontSize);
box on;

nameStr = [saveFigFolder, '\', 'est-dim-1-2-vs-degree', '_N0-', num2str(N00), '_p-', num2str(p00), '_D-', num2str(D00)];

if saveFigFlag
% savefig([nameStr, '.fig']);
% exportgraphics(f, [nameStr, '.png'], 'BackgroundColor', 'None', 'ContentType', 'vector');
exportgraphics(f, [nameStr, '.pdf'], 'BackgroundColor', 'None', 'ContentType', 'vector');
end

catch
    disp('Cannot plot estimated dimension versus degree.');
end



% Network size  (like Fig. 3):

try

f = figure; hold on;

for iiLoadDataFolder = 1:numLoadDataFolder

loadDataFolder = loadDataFolderCell{iiLoadDataFolder};
styleCell = styleCellCell{iiLoadDataFolder};

if (k00 == 2*D00)
    switch iiLoadDataFolder
        case 2
            loadStr = [loadDataFolder, '\', 'est-dim_small-world-manhattan_DLim-', num2str(DLims(1)), '-', num2str(DLims(2)), '_', num2str(numRep), 'rep', '_k-', num2str(k00), '_D-', num2str(D00), '_p-', num2str(p00), '_10N0-100-100000'];
        case 1
            loadStr = [loadDataFolder, '\', 'est-dim-2_small-world-manhattan_DLim-', num2str(DLims(1)), '-', num2str(DLims(2)), '_', num2str(numRep), 'rep', '_k-', num2str(k00), '_D-', num2str(D00), '_p-', num2str(p00), '_10N0-100-100000'];
    end
else
    switch iiLoadDataFolder
        case 2
            loadStr = [loadDataFolder, '\', 'est-dim_small-world-manhattan_DLim-', num2str(DLims(1)), '-', num2str(DLims(2)), '_', num2str(numRep), 'rep', '_k-', num2str(k00), '_D-', num2str(D00), '_p-', num2str(p00), '_7N0-100-10000'];
        case 1
            loadStr = [loadDataFolder, '\', 'est-dim-2_small-world-manhattan_DLim-', num2str(DLims(1)), '-', num2str(DLims(2)), '_', num2str(numRep), 'rep', '_k-', num2str(k00), '_D-', num2str(D00), '_p-', num2str(p00), '_7N0-100-10000'];
    end
end
load([loadStr, '.mat']);

xx = NN;
set(gca, 'XMinorTick', 'off', 'TickLabelInterpreter', 'LaTeX', 'FontSize', fontSize);
Y = mean(dimMatNN, 3);
YStd = std(dimMatNN, 0, 3);
YErr = YStd/sqrt(numRep - 1);
for iiPlot = 1:numToPlot
    iiMeth = toPlot(iiPlot);
    if (iiMeth == 1)
        if (iiLoadDataFolder == 2); continue; end
        yy = D*ones(size(xx));
        yyErr = zeros(size(xx));
        yyStd = zeros(size(xx));
    else
        yy = Y(iiMeth - 1, :);
        yyErr = YErr(iiMeth - 1, :);
        yyStd = YStd(iiMeth - 1, :);
    end
    col = colMat(iiMeth, :);
    sty = styleCell{iiMeth};
    p = plot(xx, yy, sty, 'Color', col, 'LineWidth', lineWidth, 'MarkerSize', markerSize);
    if (iiLoadDataFolder == 1)
        p.MarkerFaceColor = col;
    else
        p.MarkerFaceColor = 'w';
    end
end
end

xlim([min(xx), max(max(xx), N00)]);
yLim = ylim; yMid = D; yHalfRange = max(abs(yLim - yMid));
yHalfRange = max(yHalfRange, 0.2);
yLim = [yMid - yHalfRange, yMid + yHalfRange]; ylim(yLim);
if (diff(yLim) < 0.4)
    yTicks = yticks; yTicks = (round(10*min(yTicks))/10):0.1:max(yTicks); yticks(yTicks);
end
ylabel('Estimated dimension, $\hat{D}$', 'Interpreter', 'LaTeX', 'FontSize', fontSize);
xlabel('Network size, $N$', 'Interpreter', 'LaTeX', 'FontSize', fontSize);
set(gca, 'XScale', 'Log');
set(gca, 'XTick', [10, 100, 1000, 10000, 100000]);
if p > 0
    title(['$D = ', num2str(D00), ', \langle k \rangle = ', num2str(k00), '$'], 'Interpreter', 'LaTeX', 'FontSize', fontSize);
else
    title(['$D = ', num2str(D00), ', k = ', num2str(k00), '$'], 'Interpreter', 'LaTeX', 'FontSize', fontSize);
end
box on;

ax = gca; ax.Position(3) = 0.77;

nameStr = [saveFigFolder, '\', 'est-dim-1-2-vs-network-size', '_p-', num2str(p00), '_D-', num2str(D00), '_k-', num2str(k00)];

if saveFigFlag
% savefig([nameStr, '.fig']);
% exportgraphics(f, [nameStr, '.png'], 'BackgroundColor', 'None', 'ContentType', 'vector');
exportgraphics(f, [nameStr, '.pdf'], 'BackgroundColor', 'None', 'ContentType', 'vector');
end

catch
    disp('Cannot plot estimated dimension versus network size.');
end



% Rate of deletion  (like Fig. 6):

try

f = figure; hold on;

for iiLoadDataFolder = 1:numLoadDataFolder

loadDataFolder = loadDataFolderCell{iiLoadDataFolder};
styleCell = styleCellCell{iiLoadDataFolder};

switch iiLoadDataFolder
    case 2
        loadStr = [loadDataFolder, '\', 'est-dim_small-world-manhattan_DLim-', num2str(DLims(1)), '-', num2str(DLims(2)), '_', num2str(numRep), 'rep', '_k-', num2str(k00), '_D-', num2str(D00), '_p-0', '_N0-', num2str(N0), '_21f-0-1', '_lcc'];
    case 1
        loadStr = [loadDataFolder, '\', 'est-dim-2_small-world-manhattan_DLim-', num2str(DLims(1)), '-', num2str(DLims(2)), '_', num2str(numRep), 'rep', '_k-', num2str(k00), '_D-', num2str(D00), '_p-0', '_N0-', num2str(N0), '_21f-0-1', '_lcc'];
end

load(loadStr);

xx = ff;
Y = mean(dimMat_ff, 3, 'omitnan');
YStd = std(dimMat_ff, 0, 3, 'omitnan');
YErr = YStd./sqrt(sum(~isnan(dimMat_ff), 3) - 1);
for iiPlot = 1:numToPlot
    iiMeth = toPlot(iiPlot);
    if (iiMeth == 1)
        if (iiLoadDataFolder == 2); continue; end
        yy = D*ones(size(xx));
        yyErr = zeros(size(xx));
        yyStd = zeros(size(xx));
    else
        yy = Y(iiMeth - 1, :);
        yyErr = YErr(iiMeth - 1, :);
        yyStd = YStd(iiMeth - 1, :);
    end
    col = colMat(iiMeth, :);
    sty = styleCell{iiMeth};    p = plot(xx, yy, sty, 'Color', col, 'LineWidth', lineWidth, 'MarkerSize', markerSize);
    if (iiLoadDataFolder == 1)
        p.MarkerFaceColor = col;
    else
        p.MarkerFaceColor = 'w';
    end
end

end

set(gca, 'XMinorTick', 'off', 'TickLabelInterpreter', 'LaTeX', 'FontSize', fontSize);

xlim([min(xx), max(xx)]);
yLim = ylim; yMid = D; yHalfRange = max(abs(yLim - yMid)); yLim = [yMid - yHalfRange, yMid + yHalfRange]; ylim(yLim);
ylabel('Estimated dimension, $\hat{D}$', 'Interpreter', 'LaTeX', 'FontSize', fontSize);
xlabel('Link deletion rate, $f$', 'Interpreter', 'LaTeX', 'FontSize', fontSize);
N = round(N0^(1/D00))^D00;
if p > 0
    title(['$D = ', num2str(D00), ', \langle k \rangle = ', num2str(k00), ', N = ', num2str(N), '$'], 'Interpreter', 'LaTeX', 'FontSize', fontSize);
else
    title(['$D = ', num2str(D00), ', k = ', num2str(k00), ', N = ', num2str(N), '$'], 'Interpreter', 'LaTeX', 'FontSize', fontSize);
end
% set(gca, 'XScale', 'Log');
set(gca, 'XTick', 0:0.2:1);
box on;

nameStr = [saveFigFolder, '\', 'est-dim-1-2-vs-deletion-rate', '_N0-', num2str(N00), '_D-', num2str(D00), '_k-', num2str(k00)];

if saveFigFlag
% savefig([nameStr, '.fig']);
% exportgraphics(f, [nameStr, '.png'], 'BackgroundColor', 'None', 'ContentType', 'vector');
exportgraphics(f, [nameStr, '.pdf'], 'BackgroundColor', 'None', 'ContentType', 'vector');
end

catch
    disp('Cannot plot estimated dimension versus rate of deletion.');
end




% Make legends:

% Legend:

f = figure; hold on;
% True:
p = plot(NaN, NaN, '-', 'Color', colMat(1, :), 'LineWidth', lineWidth, 'MarkerSize', markerSize);
% Fitting to corr. int.:
p = plot(NaN, NaN, '--o', 'Color', 'k', 'LineWidth', lineWidth, 'MarkerSize', markerSize, 'MarkerFaceColor', 'k');
% Fitting to corr.:
p = plot(NaN, NaN, ':o', 'Color', 'k', 'LineWidth', lineWidth, 'MarkerSize', markerSize, 'MarkerFaceColor', 'w');
fr = 0.5;
% Lin. regression:
p = plot(NaN, NaN, '-.h', 'Color', colMat(4, :), 'LineWidth', lineWidth, 'MarkerSize', markerSize, 'MarkerFaceColor', (1 - fr) + fr*colMat(4, :));
% KS distance:
p = plot(NaN, NaN, '-.d', 'Color', colMat(3, :), 'LineWidth', lineWidth, 'MarkerSize', markerSize, 'MarkerFaceColor', (1 - fr) + fr*colMat(3, :));
% Code rate:
p = plot(NaN, NaN, '-.s', 'Color', colMat(2, :), 'LineWidth', lineWidth, 'MarkerSize', markerSize, 'MarkerFaceColor', (1 - fr) + fr*colMat(2, :));

leg = legend({'True', 'Fitting $C(s) \propto s^D$', 'Fitting $c(s) \propto s^{D-1}$', 'Lin. regression', 'KS distance', 'Neg. log-like. per obs.'}, 'Interpreter', 'LaTeX', 'FontSize', fontSize, 'Location', 'Best', 'Orientation', 'Horizontal', 'NumColumns', 3);
leg.BoxFace.ColorType='truecoloralpha'; leg.BoxFace.ColorData=uint8(255*[1; 1; 1; legendBoxAlpha]);

axis off;

nameStr = [saveFigFolder, '\', 'est-dim-1-2_legend_3col'];

f.Position([1, 3]) = [50, 1000];

if saveFigFlag
% savefig([nameStr, '.fig']);
% exportgraphics(f, [nameStr, '.png'], 'BackgroundColor', 'None', 'ContentType', 'vector');
exportgraphics(f, [nameStr, '.pdf'], 'BackgroundColor', 'None', 'ContentType', 'vector');
end




% Legend:

f = figure; hold on;
% True:
p = plot(NaN, NaN, '-', 'Color', colMat(1, :), 'LineWidth', lineWidth, 'MarkerSize', markerSize);
% Fitting to corr. int.:
p = plot(NaN, NaN, '--o', 'Color', 'k', 'LineWidth', lineWidth, 'MarkerSize', markerSize, 'MarkerFaceColor', 'k');
% Fitting to corr.:
p = plot(NaN, NaN, ':o', 'Color', 'k', 'LineWidth', lineWidth, 'MarkerSize', markerSize, 'MarkerFaceColor', 'w');
fr = 0.5;
% Lin. regression:
p = plot(NaN, NaN, '-.h', 'Color', colMat(4, :), 'LineWidth', lineWidth, 'MarkerSize', markerSize, 'MarkerFaceColor', (1 - fr) + fr*colMat(4, :));
% KS distance:
p = plot(NaN, NaN, '-.s', 'Color', colMat(3, :), 'LineWidth', lineWidth, 'MarkerSize', markerSize, 'MarkerFaceColor', (1 - fr) + fr*colMat(3, :));
% Code rate:
p = plot(NaN, NaN, '-.d', 'Color', colMat(2, :), 'LineWidth', lineWidth, 'MarkerSize', markerSize, 'MarkerFaceColor', (1 - fr) + fr*colMat(2, :));

% leg = legend({'True', 'Fitting $C(s) \propto s^D$', 'Fitting $c(s) \propto s^{D-1}$', 'Lin. regression', 'KS distance', 'Code rate'}, 'Interpreter', 'LaTeX', 'FontSize', fontSize, 'Location', 'Best', 'Orientation', 'Vertical', 'NumColumns', 1);
leg = legend({'True', 'Fitting $C(s) \propto s^D$', 'Fitting $c(s) \propto s^{D-1}$', 'Lin. regression', 'KS distance', 'Neg. log-like. per obs.'}, 'Interpreter', 'LaTeX', 'FontSize', fontSize, 'Location', 'Best', 'Orientation', 'Vertical', 'NumColumns', 1);
leg.BoxFace.ColorType='truecoloralpha'; leg.BoxFace.ColorData=uint8(255*[1; 1; 1; legendBoxAlpha]);
% leg.Box = 'Off';

f.Position([1, 3]) = [50, 1000];

leg.Position([2, 4]) = [0, 1];

axis off;

nameStr = [saveFigFolder, '\', 'est-dim-1-2_legend_1col'];

if saveFigFlag
% savefig([nameStr, '.fig']);
% exportgraphics(f, [nameStr, '.png'], 'BackgroundColor', 'None', 'ContentType', 'vector');
exportgraphics(f, [nameStr, '.pdf'], 'BackgroundColor', 'None', 'ContentType', 'vector');
end
