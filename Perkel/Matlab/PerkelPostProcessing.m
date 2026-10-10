%% PerkelPostProcessing.m
% MATLAB port of PerkelPostProcessing.nb: extract (intensity, LO phase,
% convolution phase) from every oscilloscope trace and write data1a.mat,
% data2a.mat and data3a.mat.
%
% The raw scope traces (D:\pcdata\data2 ... data5) are not in the
% repository. Set rawDir to their location before running. Output files
% are written to outDir (default: this Matlab/ folder), so the original
% data*.mat files in Perkel/ are not overwritten. They use the same layout
% as the Mathematica export (variables Expression1, Expression2, ...), so
% PerkelData.m reads either.
%
% Set makePlots = false before running to skip the figure.

if ~exist('makePlots', 'var'), makePlots = true; end
if ~exist('rawDir', 'var'), rawDir = 'D:\pcdata'; end
if ~exist('outDir', 'var')
    outDir = fileparts(mfilename('fullpath'));
    if isempty(outDir), outDir = pwd; end
end
if ~exist(rawDir, 'dir')
    error('PerkelPostProcessing:noRawData', ...
          'Raw scope traces not found in %s; set rawDir to their location.', rawDir);
end
traceFile = @(sub, run, i, j) fullfile(rawDir, sub, sprintf('%s___%d_%d.csv', run, i, j));
readTrace = @(f) readmatrix(f);

%% How to dig the data out of the waveform
sampledataraw = readTrace(traceFile('data3', 'data6557', 5, 4));
sampledataraw = sampledataraw(50:9950, 1) + 33;
if makePlots
    figure('Name', 'Sample trace');
    subplot(1, 2, 1); plot(sampledataraw); axis tight;
    w = zeros(300, 1);
    for k = 1:300
        seg = sampledataraw(1 + 33*(k - 1):33*k);
        w(k) = max(seg) + min(seg);
    end
    subplot(1, 2, 2); plot(w); hold on; plot(diff(w)); axis tight;
end

%% Parse function (see parseWave.m)
disp(parseWave(readTrace(traceFile('data3', 'data6557', 278, 1))));

%% Generate files
% lv1: states 1..277, ten repetitions in each of two runs
datalv1 = cat(2, parseRun(traceFile, readTrace, 'data2', 'data9157', 0, 277), ...
                 parseRun(traceFile, readTrace, 'data3', 'data6557', 0, 277));
saveLikeMathematica(fullfile(outDir, 'data1a.mat'), datalv1);

% lv2: states 278..317, three runs
datalv2 = cat(2, parseRun(traceFile, readTrace, 'data3', 'data6557', 277, 40), ...
                 parseRun(traceFile, readTrace, 'data4', 'data8491', 277, 40), ...
                 parseRun(traceFile, readTrace, 'data5', 'data9340', 277, 40));
saveLikeMathematica(fullfile(outDir, 'data2a.mat'), datalv2);

% lv3: states 318..1516, one run
datalv3 = parseRun(traceFile, readTrace, 'data4', 'data8491', 317, 1199);
saveLikeMathematica(fullfile(outDir, 'data3a.mat'), datalv3);
