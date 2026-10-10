function data = loadMathematicaMat(fileName)
%LOADMATHEMATICAMAT  Import["dataNa.mat"] as Mathematica returns it.
%   The data files were written by Mathematica's Export, which stores the
%   list {m1, m2, ...} as variables Expression1, Expression2, ... Each mk
%   is an R-by-3 matrix (one row per repetition: intensity, LO phase,
%   convolution phase). DATA is the S-by-R-by-3 array, DATA(k, :, :) = mk.
s = load(fileName);
names = fieldnames(s);
idx = cellfun(@(nm) sscanf(nm, 'Expression%d'), names);
[~, order] = sort(idx);
first = s.(names{order(1)});
data = zeros(numel(names), size(first, 1), size(first, 2));
for k = 1:numel(order)
    data(k, :, :) = s.(names{order(k)});
end
end
