function [sorted, pos] = genMajorSeq(x)
%GENMAJORSEQ  Port of genMajorSeq[] from PerkelTheory.nb.
%   Mathematica:
%     genMajorSeq[x_] := Module[{sorted = Sort[x, Abs[#1] < Abs[#2] &]},
%       {sorted, Flatten[Position[x, #] & /@ DeleteDuplicates@sorted, 2]}];
%
%   SORTED is X ordered by increasing absolute value. POS lists, for each
%   distinct value of SORTED (in order of appearance), every index of X that
%   holds it. Example: genMajorSeq([1 4 2 2 -3 7]) gives
%   sorted = [1 2 2 -3 4 7] and pos = [1 3 4 5 2 6].
%
%   Tie handling matters for the pulse tables: Mathematica's Sort with the
%   strict test Abs[#1] < Abs[#2] puts elements of equal magnitude in the
%   reverse of their original order, and the rays contain entries such as
%   +a and -a whose magnitudes differ only by rounding error. Magnitudes
%   that agree to TOL (relative) are therefore treated as equal, and ties
%   are ordered last-first. This reproduces the notebook's impulse.txt.

tol = 1e-9;
x = x(:).';
n = numel(x);
[a, o] = sort(abs(x));
grp = zeros(1, n);
grp(o) = cumsum([1, diff(a) > tol * max(1, a(2:end))]);
[~, idx] = sortrows([grp.', -(1:n).']);
sorted = x(idx);

if nargout > 1
    % DeleteDuplicates keeps first occurrences, in order
    [~, first] = unique(sorted, 'first');
    vals = sorted(sort(first));
    pos = [];
    for v = vals
        pos = [pos, find(x == v)]; %#ok<AGROW>
    end
end
end
