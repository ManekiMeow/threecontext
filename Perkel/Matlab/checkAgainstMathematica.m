%% checkAgainstMathematica.m
% Runs PerkelTheory.m and PerkelData.m without figures and compares their
% results with the outputs saved in the Mathematica notebooks and with the
% pulse tables impulse.txt / pmpulse.txt in Perkel/.
%
% Bootstrapped errors use a different random generator from Mathematica's,
% so they are only checked loosely.

clear; makePlots = false;
here = fileparts(mfilename('fullpath'));
if isempty(here), here = pwd; end
run(fullfile(here, 'PerkelData.m'));
orig = fullfile(here, '..');

checks = {
    'MatrixRank[gmat]',          rank(gmat),                    37,                  0
    'Length@statelistlv1',       size(statelistlv1, 1),         277,                 0
    'Length@statelistlv2',       size(statelistlv2, 1),         40,                  0
    'Length@errlist',            size(errlist, 1),              1199,                0
    'Max@pulseerr',              max(pulseerr(:)),              0.1666666666666667,  1e-12
    'eff(2)',                    eff(2),                        0.7830582405624588,  1e-9
    'eff(12)',                   eff(12),                       0.055095317128106516, 1e-9
    'min Length data1f',         min(cellfun(@(c) size(c, 1), data1f)), 12,          0
    'min Length data2f',         min(cellfun(@(c) size(c, 1), data2f)), 24,          0
    'min Length data3f',         min(cellfun(@(c) size(c, 1), data3f)), 6,           0
    'mean LO phase / deg',       mean(allPhase(2)) * 180/pi,    2.993568060040112,   1e-9
    'mean conv. phase / deg',    mean(allPhase(3)) * 180/pi,    4.462341563435254,   1e-9
    'std LO phase / deg',        std(allPhase(2)) * 180/pi,     3.9358044877797393,  1e-9
    'std conv. phase / deg',     std(allPhase(3)) * 180/pi,     2.7375711153450344,  1e-9
    'hardyprob(1)',              hardyprob(1),                  0.993895039720542,   1e-12
    'hardyprob(2)',              hardyprob(2),                  0.9980003072997136,  1e-12
    'hardyprob(3)',              hardyprob(3),                  0.9983409513934605,  1e-12
    'penalty',                   penalty,                       0.650940342556335,   1e-12
    'Total@hardyprob - penalty', sum(hardyprob) - penalty,      2.339295955857381,   1e-12
    'Correlation',               cc(1, 2),                      0.9957714337509325,  1e-12
    'Fit intercept',             pfit(2),                       7.767027518670593,   1e-9
    'Fit slope',                 pfit(1),                       306.5939863832945,   1e-9
    'penaltyerr (bootstrap)',    penaltyerr,                    0.040063858396451185, 0.005
    'violation / sigma (bootstrap)', (sum(hardyprob) - penalty - 2) / (sum(hardyerr) + penaltyerr), 8.057124326891413, 0.5
    };

label = {'FAIL', 'ok'};
ok = true;
for k = 1:size(checks, 1)
    [name, got, want, tol] = checks{k, :};
    pass = abs(got - want) <= tol;
    ok = ok && pass;
    fprintf('%-32s %-6s got %.15g, Mathematica %.15g\n', name, label{pass + 1}, got, want);
end

A = load(fullfile(orig, 'impulse.txt')); B = load(fullfile(here, 'impulse.txt'));
pass = isequal(size(A), size(B)) && max(abs(A(:) - B(:))) < 1e-4;
ok = ok && pass;
fprintf('%-32s %s\n', 'impulse.txt identical', label{pass + 1});
A = load(fullfile(orig, 'pmpulse.txt')); B = load(fullfile(here, 'pmpulse.txt'));
pass = isequal(A, B);
ok = ok && pass;
fprintf('%-32s %s\n', 'pmpulse.txt identical', label{pass + 1});

if ok, disp('All checks passed.'); else, error('Some checks failed.'); end
