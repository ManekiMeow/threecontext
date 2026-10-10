%% PerkelData.m
% MATLAB port of PerkelData.nb: phase screening, Hardy-type paradox,
% violation of the NCHV inequality with bootstrapped errors, amplitude
% correction for the level-2 vectors, and the figures.
%
% Requires PerkelTheory.m to run (it is run automatically if its variables
% are not in the workspace). Reads data1a.mat, data2a.mat and data3a.mat
% from dataDir (default: the Perkel/ folder above this one).
%
% Set makePlots = false before running to skip the figures.

if ~exist('makePlots', 'var'), makePlots = true; end
dataScriptDir = fileparts(mfilename('fullpath'));
if isempty(dataScriptDir), dataScriptDir = pwd; end
if ~exist('stateposlv2', 'var') || ~exist('errpos', 'var')
    run(fullfile(dataScriptDir, 'PerkelTheory.m'));
end
if ~exist('dataDir', 'var'), dataDir = fullfile(dataScriptDir, '..'); end

data1 = loadMathematicaMat(fullfile(dataDir, 'data1a.mat'));
data2 = loadMathematicaMat(fullfile(dataDir, 'data2a.mat'));
data3 = loadMathematicaMat(fullfile(dataDir, 'data3a.mat'));
data3(1, :, 1) = 0;

fprintf('Mean/@data2[[;;, ;;, 1]] = %s\n', mat2str(mean(data2(:, :, 1), 2).', 6));

%% Phase extraction and screening
data1f = filterPhase(data1);
fprintf('(Length/@data1f)//Min = %d\n', min(cellfun(@(c) size(c, 1), data1f)));
data2f = filterPhase(data2);
fprintf('(Length/@data2f)//Min = %d\n', min(cellfun(@(c) size(c, 1), data2f)));
data3f = filterPhase(data3);
fprintf('(Length/@data3f)//Min = %d\n', min(cellfun(@(c) size(c, 1), data3f)));

% Join[data1, data2, data3] has a different number of repetitions per
% level, so keep the three levels in a cell array.
allData = {data1, data2, data3};
if makePlots
    th = pi/24;
    stateMeans = cell2mat(cellfun(@(d) squeeze(mean(d(:, :, 2:3), 2)), allData(:), 'UniformOutput', false));
    avg = mean(stateMeans, 1);
    pts = cell2mat(cellfun(@(d) squeeze(d(:, 7, 2:3)), allData(:), 'UniformOutput', false)) - avg;
    figure('Name', 'Phase screening'); hold on; axis equal;
    rectangle('Position', [-.2 -.2 .4 .4], 'FaceColor', [0.85 0.85 0.85], 'EdgeColor', 'none');
    rectangle('Position', [-th -th 2*th 2*th], 'FaceColor', 'w', 'EdgeColor', 'none');
    plot([-.2 .2; -.2 .2].', [-th -th; th th].', 'Color', [0.5 0.5 0.5], 'LineWidth', 2);
    plot([-th -th; th th].', [-.2 .2; -.2 .2].', 'Color', [0.5 0.5 0.5], 'LineWidth', 2);
    inPlot = all(abs(pts) < .196, 2);
    inside = all(abs(pts) < th, 2);
    plot(pts(inPlot & inside, 1), pts(inPlot & inside, 2), '.', 'Color', [0 105 169]/255, 'MarkerSize', 8);
    plot(pts(inPlot & ~inside, 1), pts(inPlot & ~inside, 2), '.', 'Color', [243 90 93]/255, 'MarkerSize', 8);
    rectangle('Position', [-.2 -.2 .4 .4], 'EdgeColor', 'k', 'LineWidth', 2);
    xlim([-.2 .2]); ylim([-.2 .2]);
    set(gca, 'XTick', (-2:2)*pi/36, 'XTickLabel', 5*(-2:2), 'YTick', (-2:2)*pi/36, 'YTickLabel', 5*(-2:2));
    xlabel('LO phase / Degree'); ylabel('Convolution phase / Degree');
end

allPhase = @(c) cell2mat(cellfun(@(d) reshape(d(:, :, c), [], 1), allData(:), 'UniformOutput', false));
fprintf('Mean LO phase / Degree = %.15g\n', mean(allPhase(2)) * 180/pi);
fprintf('Mean convolution phase / Degree = %.15g\n', mean(allPhase(3)) * 180/pi);
fprintf('StdDev LO phase / Degree = %.15g\n', std(allPhase(2)) * 180/pi);
fprintf('StdDev convolution phase / Degree = %.15g\n', std(allPhase(3)) * 180/pi);

%% Hardy-type paradox
meanIntensity = @(cells) cellfun(@(c) mean(c(:, 1)), cells);
chunks = {1:37, 38:74, 75:111};
hardyFun = @(m) cellfun(@(ix) sum(m(ix(1:19)).^2) / sum(m(ix).^2), chunks);
m2all = meanIntensity(data2f);
m2 = m2all(stateposlv2);
hardyprob = hardyFun(m2);
fprintf('hardyprob = %s\n', mat2str(hardyprob, 15));
fprintf('Total@hardyprob = %.15g\n', sum(hardyprob));
hardyterms = cell2mat(cellfun(@(ix) (m2(ix(1:19)).^2 / sum(m2(ix).^2)).', chunks.', 'UniformOutput', false));

%% Violation of NCHV inequality
estnorm = sum(m2) / 57;
% Means are taken per measured state first, then picked out by errpos.
penaltyFun = @(m3) sum((sum(m3(errpos), 2) / estnorm).^2) / 19 / 2;
penalty = penaltyFun(meanIntensity(data3f));
fprintf('penalty = %.15g\n', penalty);
fprintf('Total@hardyprob - penalty = %.15g\n', sum(hardyprob) - penalty);
nComplementEdges = nnz(gmat == 0) / 2;      % EdgeCount[GraphComplement@Perkel]
% The notebook evaluates this with the penalty and its error typed in.
fprintf('per-edge penalty = %s\n', mat2str([0.650940342556335, 0.040063858396451185] * 19 * 2 / nComplementEdges, 15));

%% Bootstrapping for error estimation
% SeedRandom[8023] seeds Mathematica's generator; MATLAB's rng(8023) gives
% different draws, so these errors agree statistically, not digit by digit.
rng(8023);
rnd = zeros(200, 3);
for j = 1:200
    seq = cellfun(@(c) c(randperm(size(c, 1), 20), :), data2f, 'UniformOutput', false);
    mseq = meanIntensity(seq);
    rnd(j, :) = hardyFun(mseq(stateposlv2));
end
hardyerr = std(rnd);
fprintf('hardyerr = %s\n', mat2str(hardyerr, 6));

rng(8023);
rnd = zeros(200, 1);
for j = 1:200
    seq = cellfun(@(c) c(randperm(size(c, 1), min(size(c, 1), 9)), :), data3f, 'UniformOutput', false);
    rnd(j) = penaltyFun(meanIntensity(seq));
end
penaltyerr = std(rnd);
fprintf('penaltyerr = %.6g\n', penaltyerr);

fprintf('Total@hardyerr + penaltyerr = %.6g\n', sum(hardyerr) + penaltyerr);
fprintf('Violation / sigma = %.6g\n', (sum(hardyprob) - penalty - 2) / (sum(hardyerr) + penaltyerr));
fprintf('per-edge penalty (this run) = %s\n', mat2str([penalty, penaltyerr] * 19 * 2 / nComplementEdges, 6));

%% Amplitude correction for LV2 vectors
m1 = meanIntensity(data1f);
dataq = m1(stateposlv1);                 % 111 x 6
dataq(:, 4:end) = 0;
prlv1new = zeros(111, 6);
flatChunks = {1:222, 223:444, 445:666};
dataqFlat = reshape(dataq.', [], 1);     % Flatten (row by row)
for c = 1:3
    prlv1new(chunks{c}, :) = dataq(chunks{c}, :) / sqrt(sum(dataqFlat(flatChunks{c}).^2));
end
cc = corrcoef(sum(prlv1new, 2), sum(prlv1, 2));
fprintf('Correlation = %.15g\n', cc(1, 2));

%% Plots
% Mathematica's indexed colour schemes ColorData[75], [55], [51], [24] and
% "SunsetColors" have no MATLAB equivalent; similar palettes are used.
pal75 = [0.263 0.341 0.588; 0.235 0.475 0.690; 0.282 0.604 0.745; 0.380 0.718 0.765; ...
         0.541 0.808 0.773; 0.722 0.867 0.796; 0.471 0.553 0.737; 0.341 0.420 0.655];
pal55 = [0.851 0.373 0.008; 0.906 0.541 0.180; 0.945 0.678 0.373; 0.969 0.792 0.565; ...
         0.808 0.443 0.098; 0.737 0.333 0.035];
pal51 = [0.400 0.651 0.118; 0.553 0.749 0.267; 0.686 0.831 0.427; 0.808 0.902 0.608; ...
         0.302 0.553 0.090; 0.471 0.706 0.200; 0.616 0.792 0.349; 0.745 0.867 0.522];
sunset = [0 0 0 0; .25 0.33 0.10 0.45; .5 0.80 0.18 0.25; .75 0.98 0.56 0.13; 1 1 0.96 0.75];

if makePlots
    h = .23;
    acc = [zeros(3, 1), cumsum(hardyterms, 2)];
    stack = [0; hardyprob(1); hardyprob(1) + hardyprob(2)] + acc;
    figure('Name', 'Hardy-type paradox'); hold on;
    rectangle('Position', [0 -.1 2 .35 + h], 'FaceColor', [0.85 0.85 0.85], 'EdgeColor', 'none');
    for k = 1:19
        rectangle('Position', [stack(1, k) 0 stack(1, k + 1) - stack(1, k) h], 'FaceColor', pal75(mod(k - 1, 8) + 1, :), 'EdgeColor', 'none');
        rectangle('Position', [stack(2, k) 0 stack(2, k + 1) - stack(2, k) h], 'FaceColor', pal55(mod(k + 2, 6) + 1, :), 'EdgeColor', 'none');
        rectangle('Position', [stack(3, k) 0 stack(3, k + 1) - stack(3, k) h], 'FaceColor', pal51(mod(k - 1, 8) + 1, :), 'EdgeColor', 'none');
    end
    plot([2 2], [-.1 .25 + h], '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 2);
    rectangle('Position', [0 -.1 3 .35 + h], 'EdgeColor', 'k', 'LineWidth', 2);
    lo = sum(hardyprob) - sum(hardyerr) - penaltyerr; hi = sum(hardyprob) + sum(hardyerr) + penaltyerr;
    plot([lo hi], [h/2 h/2], 'r', [lo lo], h/2 + [-.08 .08], 'r', [hi hi], h/2 + [-.08 .08], 'r', 'LineWidth', 2);
    rectangle('Position', [1.98 h + .08 .04 .1], 'FaceColor', [0.1 0.1 0.1], 'EdgeColor', 'none');
    plot([2 1.7], [h + .13 h + .13], 'Color', [0.1 0.1 0.1], 'LineWidth', 2); text(1.48, h + .13, 'NCHV', 'HorizontalAlignment', 'center');
    purple = [73 48 146]/255;
    plot([2 2] + penalty, [-.1 .25 + h], '--', 'Color', (purple + 1)/2, 'LineWidth', 2);
    plot([2 2] + penalty, [.38 .2] + h, 'Color', (purple + 1)/2, 'LineWidth', 2);
    rectangle('Position', [1.98 + penalty h + .32 .04 .1], 'FaceColor', purple, 'EdgeColor', 'none');
    plot([2 1.7] + penalty, [h + .37 h + .37], 'Color', purple, 'LineWidth', 2);
    text(0.9 + penalty, h + .37, 'Exclusivity-corrected NCHV', 'HorizontalAlignment', 'center', 'Color', purple);
    text(3, h + .37, 'Q_{max}', 'HorizontalAlignment', 'center');
    for x = [0 2 3], text(x, -.2, num2str(x), 'HorizontalAlignment', 'center'); end
    axis off;

    m3 = meanIntensity(data3f);
    err = (sum(m3(errpos), 2) / estnorm).^2;
    [oc, orr] = find(gmat.' == 0);
    orthoreduced = [orr, oc]; orthoreduced = (orthoreduced > 19) .* (orthoreduced - 19);
    nCell = size(orthoreduced, 1);
    verts = zeros(4*nCell, 2); faces = reshape(1:4*nCell, 4, []).';
    for k = 1:nCell
        x = orthoreduced(k, 1); y = orthoreduced(k, 2);
        verts(4*k - 3:4*k, :) = [x - 1, y - 1; x, y - 1; x, y; x - 1, y];
    end
    figure('Name', 'Exclusivity errors'); hold on; axis equal off;
    patch('Faces', faces, 'Vertices', verts, 'FaceVertexCData', blendColor(sunset, 1 - err/0.05), 'FaceColor', 'flat', 'EdgeColor', 'none');
    rectangle('Position', [0 0 38 38], 'EdgeColor', 'k', 'LineWidth', 2);
    for k = 1:38
        plot([k - .5 k - .5], [0 .5], 'k'); plot([0 .5], [k - .5 k - .5], 'k');
    end
    for k = 11:10:31
        plot([k - .5 k - .5], [0 1], 'k'); plot([0 1], [k - .5 k - .5], 'k');
    end
    text(19, -5, 'Preparation', 'HorizontalAlignment', 'center');
    text(-5, 19, 'Measurement', 'Rotation', 90, 'HorizontalAlignment', 'center');
    for v = [20 30 40 50 57]
        text(v - 19, -1.5, num2str(v), 'HorizontalAlignment', 'center');
        text(-1, v - 19, num2str(v), 'HorizontalAlignment', 'right');
    end
    text(52, 35, '\times10^{-2}');
    cb = blendColor(sunset, 1 - ((1:120) - 1)/99);
    for k = 1:120
        rectangle('Position', [48 8.8 + .2*k 2 .2], 'FaceColor', cb(k, :), 'EdgeColor', 'none');
    end
    rectangle('Position', [48 9 2 24], 'EdgeColor', 'k', 'LineWidth', 2);
    for k = 1:9, plot([50 49.5], [9 + 2*k, 9 + 2*k], 'k'); end
    for v = [0 2.5 5], text(51.5, 9 + 4*v, sprintf('%.1f', v)); end
    text(50, 0, '$\square = \not\perp$', 'Interpreter', 'latex', 'HorizontalAlignment', 'center');
end

%% Detector response
pwr = [0.0721, 0.320, 0.602, 0.916, 1.266, 1.655, 1.873, 2.100, 2.329, 2.563, 2.815, 3.050, 3.306, 3.560, 3.824, 4.090, 4.363, 4.643, 4.930, 5.17, 5.45, 5.74, 6.04, 6.33, 6.63, 6.92, 7.22, 7.54];
volt = [18.6, 96.5, 186., 281.4, 397.7, 520.9, 587.2, 651.2, 730.2, 808.1, 874.4, 954.7, 1028., 1106., 1187., 1264., 1340., 1426., 1506., 1578., 1657., 1717., 1750., 1795., 1820., 1833., 1840., 1850.];

pfit = polyfit(pwr(1:20), volt(1:20), 1);       % Fit[..., {1, x}, x]
fprintf('Fit: %.15g + %.15g x\n', pfit(2), pfit(1));

if makePlots
    figure('Name', 'Detector response'); hold on; box on;
    rectangle('Position', [-1 -1 4.05 4], 'FaceColor', [0.85 0.85 0.85], 'EdgeColor', 'none');
    rectangle('Position', [3.0 1.85 .1 .3], 'FaceColor', [0.4 0.4 0.4], 'EdgeColor', 'none');
    plot([3.05 2.05], [2 2], 'Color', [0.4 0.4 0.4], 'LineWidth', 2);
    allI = cell2mat(cellfun(@(d) reshape(d(:, :, 1), [], 1), allData(:), 'UniformOutput', false));
    edges = 0:.5:3.5;
    counts = arrayfun(@(k) sum(abs(allI)/10 >= edges(k) & abs(allI)/10 < edges(k + 1)), 1:7);
    pal24 = [0.18 0.31 0.56; 0.25 0.45 0.66; 0.36 0.58 0.73; 0.52 0.70 0.78; 0.70 0.80 0.80; 0.85 0.85 0.75; 0.95 0.85 0.60];
    for k = 1:7
        if counts(k) > 0
            rectangle('Position', [.5*k - .5, -1, .5, .6*log10(counts(k)) + 1], 'FaceColor', pal24(k, :), 'EdgeColor', 'none');
        end
    end
    text(3.25, 2.2, 'Exp. Max.', 'VerticalAlignment', 'bottom');
    blue = [0 105 169]/255;
    plot(pwr, volt/1000, '.', 'Color', blue, 'MarkerSize', 18);
    xs = linspace(-.2, 8.2, 2);
    plot(xs, pfit(2)/1000 + pfit(1)/1000*xs, 'Color', 0.75*blue + 0.25, 'LineWidth', 2);
    xlim([-.2 8.2]); ylim([-.2 2.5]);
    xlabel('Power difference / \muW'); ylabel('BHD response / V');
    yyaxis right; ylim([-.2 2.5]);
    set(gca, 'YTick', [0 0.6 1.2 1.8 2.4], 'YTickLabel', {'1', '10', '10^2', '10^3', '10^4'});
    ylabel('Events registered');
end

fprintf('180/7.5 = %g\n', 180/7.5);
