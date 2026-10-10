function [u, ic] = tallyRows(M, tol)
%TALLYROWS  Distinct rows of M in order of first appearance.
%   Stand-in for Mathematica's (Tally@list)[[;;, 1]] on lists of machine
%   reals. Rows are compared after rounding to multiples of TOL (default
%   1e-12), which mimics Tally's tolerance for floating-point noise.
%   U holds the distinct rows; IC maps every row of M to its row in U, so
%   IC(k) is what Position[u, M[[k]]][[1, 1]] returns in the notebooks.

if nargin < 2, tol = 1e-12; end
keys = round(M / tol);
[~, ~, j] = unique(keys, 'rows');            % sorted order
n = size(M, 1);
firstIdx = accumarray(j(:), (1:n).', [], @min);
[~, ord] = sort(firstIdx);                   % order of first appearance
rank = zeros(numel(ord), 1);
rank(ord) = 1:numel(ord);
ic = rank(j(:));
u = M(firstIdx(ord), :);
end
