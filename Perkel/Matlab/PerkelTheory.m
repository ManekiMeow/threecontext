%% PerkelTheory.m
% MATLAB port of PerkelTheory.nb: orthogonal representation of the
% complement of the Perkel graph, state lists for the three experiment
% levels, dispersion/homodyne efficiency, and the AFG pulse tables
% (impulse.txt, pmpulse.txt).
%
% PerkelData.m needs the variables this script leaves in the workspace
% (gmat, prlv1, stateposlv1, stateposlv2, errpos, ...), just as
% PerkelData.nb required PerkelTheory.nb to have been evaluated.
%
% Set makePlots = false before running to skip the figures.

if ~exist('makePlots', 'var'), makePlots = true; end
theoryDir = fileparts(mfilename('fullpath'));
if isempty(theoryDir), theoryDir = pwd; end

%% Parameters of Perkel Graph
% Gram matrix of Complement(Perkel), see the next section. Entry 1/3 marks
% an edge of the Perkel graph, entry 0 an edge of its complement
% (orthogonal rays).
gmat = [ ...
    3 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 1 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 1 0 0 0 1 0 0 0;
    0 3 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 1 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 1 0 0 0 1 0 0;
    0 0 3 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 1 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 1 0 0 0 1 0;
    0 0 0 3 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 1 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 1 0 0 0 1;
    0 0 0 0 3 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 1 0 1 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 1 0 0 0 0 0 1 0 0 0;
    0 0 0 0 0 3 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 1 0 1 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 1 0 0 0 0 0 1 0 0;
    0 0 0 0 0 0 3 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 1 0 1 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 1 0 0 0 0 0 1 0;
    0 0 0 0 0 0 0 3 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 1 0 1 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 1 0 0 0 0 0 1;
    0 0 0 0 0 0 0 0 3 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 1 0 1 0 0 0 1 0 0 0 1 0 0 0 0 0 0 0 0 1 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 3 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 1 0 1 0 0 0 1 0 0 0 1 0 0 0 0 0 0 0 0 1 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 3 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 1 0 1 0 0 0 1 0 0 0 1 0 0 0 0 0 0 0 0 1 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 3 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 1 0 1 0 0 0 1 0 0 0 1 0 0 0 0 0 0 0 0 1 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 3 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 1 0 0 0 0 0 1 0 0 0 1 0 0 0 0 0 0 0 0 1 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 3 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 1 0 0 0 0 0 1 0 0 0 1 0 0 0 0 0 0 0 0 1;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 3 0 0 0 0 1 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 1 0 0 0 0 0 1 0 0 0 1 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 3 0 0 0 0 1 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 1 0 0 0 0 0 1 0 0 0 1 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 3 0 0 0 0 1 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 1 0 0 0 0 0 1 0 0 0 1 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 3 0 1 0 0 1 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 1 0 0 0 1 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 3 0 1 0 0 1 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 1 0 0 0 1 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 1 0 1 0 0 1 0 3 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 1 1 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 1 0 0 1 0 3 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 1 1 0 0;
    1 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 1 0 0 0 0 3 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 1 1 0;
    0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 1 0 0 0 0 3 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 1 1;
    0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 1 0 0 0 0 3 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 1;
    1 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 3 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0;
    0 1 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 3 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0;
    1 0 1 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 3 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0;
    0 1 0 1 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 3 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 0 0 0 0 0 0 0 0 0 0 1 0 0 0;
    0 0 1 0 1 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 3 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 0 0 0 0 0 0 0 0 0 0 1 0 0;
    0 0 0 1 0 1 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 3 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 0 0 0 0 0 0 0 0 0 0 1 0;
    0 0 0 0 1 0 1 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 3 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 0 0 0 0 0 0 0 0 0 0 1;
    0 0 0 0 0 1 0 1 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 3 0 0 0 0 0 0 1 0 0 0 0 0 0 1 1 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 1 0 1 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 3 0 0 0 0 0 0 1 0 0 0 0 0 0 1 1 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 1 0 1 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 3 0 0 0 0 0 0 1 0 0 0 0 0 0 1 1 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 1 0 1 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 3 0 0 0 0 0 0 1 0 0 0 0 0 0 1 1 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 1 0 1 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 3 0 0 0 0 0 0 1 0 0 0 0 0 0 1 1 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 1 0 1 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 3 0 0 0 0 0 0 1 0 0 0 0 0 0 1 1 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 1 0 1 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 3 0 0 0 0 0 0 1 0 0 0 0 0 0 1 1 0 0 0 0;
    0 0 0 0 1 0 0 0 1 0 0 0 0 0 1 0 0 0 0 0 0 0 0 1 1 0 0 0 0 0 0 1 0 0 0 0 0 0 3 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 1 0 0 0 1 0 0 0 0 0 1 0 0 0 0 0 0 0 0 1 1 0 0 0 0 0 0 1 0 0 0 0 0 0 3 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 1 0 0 0 1 0 0 0 0 0 1 0 0 0 0 0 0 0 0 1 1 0 0 0 0 0 0 1 0 0 0 0 0 0 3 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 1 0 0 0 1 0 0 0 0 0 1 0 0 0 0 0 0 0 0 1 1 0 0 0 0 0 0 1 0 0 0 0 0 0 3 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 1 0 0 0 1 0 0 0 0 0 1 0 0 0 0 0 0 0 0 1 1 0 0 0 0 0 0 1 0 0 0 0 0 0 3 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    1 0 0 0 0 0 0 0 0 1 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 0 0 0 0 0 0 1 0 0 0 0 0 0 3 0 0 0 0 0 0 0 0 0 0 0 0 0;
    0 1 0 0 0 0 0 0 0 0 1 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 0 0 0 0 0 0 1 0 0 0 0 0 0 3 0 0 0 0 0 0 0 0 0 0 0 0;
    0 0 1 0 0 0 0 0 0 0 0 1 0 0 0 1 0 0 0 1 0 0 0 0 0 0 0 0 0 0 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 3 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 1 0 0 0 0 0 0 0 0 1 0 0 0 1 0 0 0 1 0 0 0 0 0 0 0 0 0 0 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 3 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 1 0 0 0 0 0 0 0 0 1 0 0 0 1 0 0 0 1 0 0 0 0 0 0 0 0 0 0 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 3 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 1 0 0 0 0 0 0 0 0 1 0 0 0 1 0 0 0 1 0 0 0 0 0 0 0 0 0 0 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 3 0 0 0 0 0 0 0 0;
    1 0 0 0 0 0 1 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 3 0 0 0 0 0 0 0;
    0 1 0 0 0 0 0 1 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 3 0 0 0 0 0 0;
    0 0 1 0 0 0 0 0 1 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 3 0 0 0 0 0;
    0 0 0 1 0 0 0 0 0 1 0 0 0 0 0 0 0 0 1 1 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 3 0 0 0 0;
    1 0 0 0 1 0 0 0 0 0 1 0 0 0 0 0 0 0 0 1 1 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 3 0 0 0;
    0 1 0 0 0 1 0 0 0 0 0 1 0 0 0 0 0 0 0 0 1 1 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 3 0 0;
    0 0 1 0 0 0 1 0 0 0 0 0 1 0 0 0 0 0 0 0 0 1 1 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 3 0;
    0 0 0 1 0 0 0 1 0 0 0 0 0 1 0 0 0 0 0 0 0 0 1 1 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 3 ...
    ] / 3;

g = (gmat == 1/3);                 % adjacency matrix of the Perkel graph
gtest = (gmat == 0);               % adjacency matrix of GraphComplement@g
nV = size(g, 1);

% Independence number of \bar G = clique number of G. G is triangle-free
% (trace(A^3) = 0) and has edges, so the largest clique has 2 vertices.
assert(trace(double(g)^3) == 0 && any(g(:)));
alphaBar = 2;
% Chromatic number of G: the three blocks of 19 vertices are independent
% sets, so 3 colours suffice; G has odd cycles (girth 5), so 2 do not.
colourOf = repelem(1:3, 19);
assert(~any(any(g & (colourOf(:) == colourOf))));
assert(~isBipartite(g));
chiG = 3;
fprintf('alpha(\\bar G)=%d, chi(G)=%d\n', alphaBar, chiG);
% Also, we know from SDP the Lovasz theta of Complement(Perkel) is
% theta(\bar G)=3.

%% Drawing the Perkel graph
if makePlots
    [ei, ej] = find(triu(g));
    kk = repelem((1:3).', 19); ii = repmat((1:19).', 3, 1);   % k outer, i inner
    pos = [sin(2*pi*ii/19) + 2*sin(2*pi*kk/3), cos(2*pi*ii/19) + 2*cos(2*pi*kk/3)];
    figure('Name', 'Perkel graph'); hold on; axis equal off;
    plot([pos(ei, 1) pos(ej, 1)].', [pos(ei, 2) pos(ej, 2)].', 'Color', [0.5 0.5 0.5], 'LineWidth', 1);
    plot(pos(1:19, 1), pos(1:19, 2), '.', 'Color', [0 165 180]/255, 'MarkerSize', 20);
    plot(pos(20:38, 1), pos(20:38, 2), '.', 'Color', [0 85 170]/255, 'MarkerSize', 20);
    plot(pos(39:end, 1), pos(39:end, 2), '.', 'Color', 'k', 'MarkerSize', 20);

    conn = [1 2;1 3;1 4;1 5;1 6;1 10;1 12;1 14;1 15;2 3;2 4;2 5;2 6;2 9;2 11;2 13;2 16;3 4;3 7;3 8;3 10;3 12;3 13;3 16;4 7;4 8;4 9;4 11;4 14;4 15;5 6;5 7;5 8;5 9;5 12;5 14;5 16;6 7;6 8;6 10;6 11;6 13;6 15;7 8;7 10;7 11;7 14;7 16;8 9;8 12;8 13;8 15;9 10;9 11;9 12;9 13;9 14;10 11;10 12;10 13;10 14;11 12;11 15;11 16;12 15;12 16;13 14;13 15;13 16;14 15;14 16;15 16];
    kk = repelem((0.5:1:3.5).', 4); ii = repmat((1:4).', 4, 1);
    pos = [sin(2*pi*(ii + kk/2)/4) + 2.5*sin(2*pi*kk/4), cos(2*pi*(ii + kk/2)/4) + 2.5*cos(2*pi*kk/4)];
    figure('Name', '16-vertex example'); hold on; axis equal off;
    plot([pos(conn(:, 1), 1) pos(conn(:, 2), 1)].', [pos(conn(:, 1), 2) pos(conn(:, 2), 2)].', 'k', 'LineWidth', 1);
    cols = {[0 165 180]/255, [0 85 170]/255, [0 0 0], [1 0.5 0]};
    for c = 1:4
        plot(pos(4*c-3:4*c, 1), pos(4*c-3:4*c, 2), '.', 'Color', cols{c}, 'MarkerSize', 20);
    end

    conn = [1 2;2 3;3 4;4 5;1 5];
    pos = [sin(2*pi*(1:5).'/5), cos(2*pi*(1:5).'/5)];
    figure('Name', 'Pentagon'); hold on; axis equal off;
    plot([pos(conn(:, 1), 1) pos(conn(:, 2), 1)].', [pos(conn(:, 1), 2) pos(conn(:, 2), 2)].', 'k', 'LineWidth', 1);
    plot(pos([1 2], 1), pos([1 2], 2), '.', 'Color', [0 165 180]/255, 'MarkerSize', 20);
    plot(pos([3 4], 1), pos([3 4], 2), '.', 'Color', [0 85 170]/255, 'MarkerSize', 20);
    plot(pos(5, 1), pos(5, 2), '.', 'Color', 'k', 'MarkerSize', 20);
end

%% Gram matrix of Complement(Perkel)
fprintf('MatrixRank[gmat] = %d\n', rank(gmat));
% We need 37-dimensional vectors to realize the orthogonal representation
% of \bar Perkel.

%% Orthogonal representation of Complement(Perkel)
% Use Cholesky decomposition to obtain the first 37 linearly independent
% vectors. chol returns upper-triangular U with gmat(1:37,1:37) = U'*U, as
% CholeskyDecomposition does; the rays are the (normalized) columns of U.
U = chol(gmat(1:37, 1:37));
vec37 = U.';
vec37 = vec37 ./ sqrt(sum(vec37.^2, 2));
if makePlots
    figure('Name', 'SVD / vec37');
    subplot(1, 2, 1); [~, Sg, ~] = svd(gmat); imagesc(Sg); axis image; title('SingularValueDecomposition[gmat][[2]]');
    subplot(1, 2, 2); imagesc(vec37.'); axis image; title('vec37^T');
end

% Linear combination for the remaining 20 vectors:
% solve vec37 . t == gmat(k, 1:37) for each k = 38..57.
vec57 = vec37;
for k = 38:57
    t = vec37 \ gmat(k, 1:37).';
    vec57(k, :) = t.';
end
vec57 = chopArray(vec57);

% We can check the set vec57 gives exactly the desired Gram matrix.
fprintf('max |vec57*vec57'' - gmat| = %.3g\n', max(max(abs(chopArray(vec57*vec57.') - gmat))));
if makePlots
    figure('Name', 'Gram matrix check');
    subplot(1, 2, 1); imagesc(chopArray(vec57*vec57.')); axis image; title('vec57 \cdot vec57^T');
    subplot(1, 2, 2); imagesc(gmat); axis image; title('gmat');
end

%% Gram-Schmidt for orthonormal basis
% Each block of 19 rays is an orthonormal set (an independent set of the
% Perkel graph); complete it to a basis of R^37 with e_20 ... e_37.
vec111 = chopArray([gramsch(vec57(1:19, :)); gramsch(vec57(20:38, :)); gramsch(vec57(39:57, :))]);
fprintf('Dimensions[vec111] = {%d, %d}\n', size(vec111));

% The vectors themselves look ugly but they should work! Every vector is
% represented in one column, rows correspond to vector entries.
if makePlots
    figure('Name', 'vec57^T'); imagesc(vec57.'); axis image; colorbar;
    figure('Name', 'vec111^T'); imagesc(vec111.'); axis image; colorbar;
    figure('Name', 'Gram matrix of vec111'); imagesc(chopArray(vec111*vec111.')); axis image;
end
overlapInit = chopArray(vec111 * [ones(19, 1)/sqrt(19); zeros(18, 1)]);
fprintf('Overlaps with initial state (first 3): %s\n', mat2str(overlapInit(1:3).', 17));

%% Partition into subspaces
% vec42(r, j, :) is the j-th of six subspace blocks of ray r, padded to 7.
blocks = {1:7, 8:14, 15:19, 20:25, 26:31, 32:37};
vec42 = zeros(111, 6, 7);
for j = 1:6
    vec42(:, j, 1:numel(blocks{j})) = vec111(:, blocks{j});
end
% (The notebook's first ortholist, built from EdgeList@gtest, is
% overwritten below before it is used, so it is not reproduced here.)
prlv1 = chopArray([sum(vec42(:, 1:3, :), 3) / sqrt(19), zeros(111, 3)]);
prlv2 = chopArray(sum(prlv1, 2));
fprintf('prlv2 (first 3): %s\n', mat2str(prlv2(1:3).', 17));

% Orthogonal pairs (zero entries of gmat) in Position[] order, i.e. row by
% row, mapped from 1..57 to the ray index 1..111.
[oc, orr] = find(gmat.' == 0);
ortholist = [orr, oc];
ortholist = ortholist + 18*(ortholist > 19) + 18*(ortholist > 38);
erlv0 = chopArray(vec42(ortholist(:, 1), :, :) .* vec42(ortholist(:, 2), :, :));
erlv1 = sum(erlv0, 3);
erlv2 = chopArray(sum(erlv1, 2));
fprintf('Dimensions[erlv0] = {%d, %d, %d}\n', size(erlv0));
fprintf('Tally[erlv2]: %d values, all zero: %d\n', numel(erlv2), all(erlv2 == 0));

%% Dispersion and homodyne efficiency
deltalambda0 = 6e-9;
fprintf('Exp[...] at 6 nm = %g\n', exp(-(6e-9)^2 / (deltalambda0^2/log(2))));
% Pulse in time domain is ps wide
fprintf('pulse duration = %.4g s\n', (1560e-9)^2 / (3e8*deltalambda0));
% Pulse in spatial domain has 0.4 mm coherence length
l0 = (1560e-9)^2 / deltalambda0;
fprintf('l0 = %.4g m\n', l0);
% Length of fiber in small loop is roughly 8 m
deltaL = 3e8 / 75.914e6 * 3 / 1.468;
fprintf('deltaL = %.15g m\n', deltaL);
% Estimated coherence length growth per round trip
deltal = 17e-6 * deltaL * deltalambda0 * 3e8;
fprintf('deltal = %.6g m\n', deltal);

% Wavefunction overlapping incl. 20 mum misalignment, 80% loop efficiency,
% and dispersion in both signal/LO.
deltalambda = @(n) (1560e-9)^2 ./ (l0 + (n + 3)*deltal);
r = .794;
effRaw = zeros(1, 12);
for k = 1:12
    effRaw(k) = real(deltalambdacorr(k, 5, deltalambda, deltalambda0, l0)) * r^k;
end
eff = effRaw / effRaw(1);
fprintf('eff = %s\n', mat2str(eff, 6));

%% Electronic pulse
[s, p] = genMajorSeq([1 4 2 2 -3 7]);
fprintf('genMajorSeq[{1,4,2,2,-3,7}] = {%s, %s}\n', mat2str(s), mat2str(p));

sortedLv1 = zeros(111*6, 7);           % Flatten[..., 1]: ray outer, block inner
for rr = 1:111
    for j = 1:6
        sortedLv1((rr - 1)*6 + j, :) = genMajorSeq(squeeze(vec42(rr, j, :)));
    end
end
[statelistlv1, ic] = tallyRows(sortedLv1);
fprintf('Length@statelistlv1 = %d\n', size(statelistlv1, 1));
stateposlv1 = reshape(ic, 6, 111).';

revEff = fliplr(eff(1:7));
pulselistlv1 = statelistlv1 ./ revEff;
fprintf('Max@pulselistlv1 = %g\n', max(pulselistlv1(:)));

sortedLv2 = zeros(111, 7);
for rr = 1:111
    sortedLv2(rr, :) = genMajorSeq([prlv1(rr, :), 0]);
end
[statelistlv2, stateposlv2] = tallyRows(sortedLv2);
fprintf('Length@statelistlv2 = %d\n', size(statelistlv2, 1));
fprintf('stateposlv2 = %s\n', mat2str(stateposlv2.'));

pulselistlv2 = statelistlv2 ./ revEff;
fprintf('Max@pulselistlv2 = %.17g\n', max(pulselistlv2(:)));
pulselistlv2 = pulselistlv2 / max(pulselistlv2(:));
fprintf('Max@pulselistlv2 = %g\n', max(pulselistlv2(:)));

errpool = erlv0;
nPair = size(errpool, 1);
sortedErr = zeros(nPair*6, 7);
sortedErrRounded = zeros(nPair*6, 7);
for k = 1:nPair
    for j = 1:6
        v = squeeze(errpool(k, j, :));
        sortedErr((k - 1)*6 + j, :) = genMajorSeq(v);
        sortedErrRounded((k - 1)*6 + j, :) = genMajorSeq(round(v * 1e6) / 1e6);
    end
end
errlist = tallyRows(sortedErr);
% The level-3 states were measured in the order of Mathematica's errlist,
% which has 1199 entries: its Tally kept two copies of the state
% {0,0,0,0,0,0,-1/9} (entries 3 and 28) because the copies differed in the
% last bits of their floating-point values. MATLAB's rounding differs, so
% a tolerant tally finds 1198 distinct states. Re-insert the copy so that
% row k of errlist is state k of data3a.mat (and row 317+k of impulse.txt).
% errpos below always points at entry 3, never at the copy, as in the
% notebook.
assert(size(errlist, 1) == 1198 && isequal(round(errlist(3, :) * 9), [0 0 0 0 0 0 -1]));
errlist = [errlist(1:27, :); errlist(3, :); errlist(28:end, :)];
fprintf('Length@errlist = %d\n', size(errlist, 1));
% errpos: position of each (rounded) sorted block in Round[errlist, 10^-6].
% Position[...][[1, 1]] takes the first match; ismember's choice among
% repeated rows differs between MATLAB and Octave, so pick it explicitly.
errKeys = round(errlist * 1e6);
[uKeys, ~, jKeys] = unique(errKeys, 'rows');
firstOf = accumarray(jKeys, (1:size(errKeys, 1)).', [], @min);
[found, where] = ismember(round(sortedErrRounded * 1e6), uKeys, 'rows');
assert(all(found));
errpos = reshape(firstOf(where), 6, nPair).';

pulseerr = errlist ./ revEff;
fprintf('Max@pulseerr = %.16g\n', max(pulseerr(:)));

%% Generate AFG waveform
volIM = @(pulse) [ones(1, 15), -1, ones(1, 9), zeros(1, 13), pulse, zeros(1, 10), ones(1, 15)];
volPM = @() [repmat(-.9, 1, 22), repmat(.9, 1, 10), repmat(-.9, 1, 38)];
allPulses = [pulselistlv1; pulselistlv2; pulseerr];
impulse = zeros(size(allPulses, 1), 70);
for k = 1:size(allPulses, 1)
    impulse(k, :) = volIM(allPulses(k, :));
end
pmpulse = repmat(volPM(), size(impulse, 1), 1);
% SetPrecision[#, 5] -> 5 significant digits in the exported table
writeTable(fullfile(theoryDir, 'impulse.txt'), impulse);
writeTable(fullfile(theoryDir, 'pmpulse.txt'), pmpulse);

%% Electronic pulse representation
% Needs the oscilloscope trace samplewave.txt, which is not in the
% repository; set sampleWaveFile to its location to draw this figure.
if ~exist('sampleWaveFile', 'var')
    sampleWaveFile = 'E:\GooglePlex\Codes\Mathematica\2205 PerkelContext\raw\samplewave.txt';
end
samplepulse = [0, roundSignificant(pulselistlv1(112, :), 3), 0];
fprintf('samplepulse = %s\n', mat2str(samplepulse));
if makePlots && exist(sampleWaveFile, 'file')
    samplewaveraw = load(sampleWaveFile);
    samplewaveraw = reshape(samplewaveraw.', [], 1);
    sw = conv(samplewaveraw, ones(99, 1), 'valid');
    sw = sw(2750:8000);
    samplewave = interp1(1:numel(sw), sw, linspace(1, numel(sw), 1601)).' + 940;
    plotElectronicPulse(samplewave, samplepulse);
elseif makePlots
    fprintf('samplewave.txt not found; skipping the electronic pulse figure.\n');
end

%% Colour plots of vec57 and gmat
if makePlots
    stops = [0.,       0.260487, 0.356,    0.891569;
             0.166667, 0.230198, 0.499962, 0.848188;
             0.333333, 0.392401, 0.658762, 0.797589;
             0.499999, 0.964837, 0.982332, 0.98988;
             0.5,      1,        1,        1;
             0.500001, 0.95735,  0.957281, 0.896269;
             0.666667, 0.913252, 0.790646, 0.462837;
             0.833333, 0.860243, 0.558831, 0.00695811;
             1.,       1.,       0.42,     0.];
    cmap = blendColor(stops, linspace(0, 1, 256));

    figure('Name', 'vec57'); image(1:57, 1:37, reshape(blendColor(stops, 1/2 + vec57.'/2), 37, 57, 3));
    axis image; xlabel('Ray #'); ylabel('Entry #');
    set(gca, 'XTick', [1 10 20 30 40 50 57], 'YTick', [1 10 20 30 37]);

    figure('Name', 'gmat'); image(1:57, 1:57, reshape(blendColor(stops, 1/2 + gmat/2), 57, 57, 3));
    axis image; xlabel('Ray #1'); ylabel('Ray #2');
    set(gca, 'XTick', [1 10 20 30 40 50 57], 'YTick', [1 10 20 30 40 50 57]);

    figure('Name', 'colour bar'); colormap(cmap); caxis([-1 1]); colorbar; axis off;
end
