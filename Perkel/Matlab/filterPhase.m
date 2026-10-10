function out = filterPhase(data, th)
%FILTERPHASE  Port of filterPhase[] from PerkelData.nb.
%   filterPhase[data_, th_:Pi/24] := Module[
%     {avg = Mean@Mean@data[[;;, ;;, {2, 3}]], temp},
%     temp = Map[# - Join[{0}, avg] &, data, {2}];
%     Select[#, (Abs[#[[2]]] < th && Abs[#[[3]]] < th &)] & /@ temp];
%
%   DATA is S-by-R-by-3. Both phases are centred on their mean over all
%   states and repetitions, then only repetitions with both phases within
%   TH of zero are kept. OUT is an S-by-1 cell array of (kept)-by-3
%   matrices.
if nargin < 2, th = pi/24; end
avg = squeeze(mean(mean(data(:, :, 2:3), 1), 2)).';
out = cell(size(data, 1), 1);
for k = 1:size(data, 1)
    temp = squeeze(data(k, :, :));
    if size(data, 2) == 1, temp = temp.'; end
    temp = temp - [0, avg];
    out{k} = temp(abs(temp(:, 2)) < th & abs(temp(:, 3)) < th, :);
end
end
