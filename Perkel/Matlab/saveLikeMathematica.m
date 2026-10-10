function saveLikeMathematica(fileName, data)
%SAVELIKEMATHEMATICA  Export[fileName, data] for an S-by-R-by-3 array, in
%   the layout Mathematica uses: one R-by-3 variable ExpressionK per state.
s = struct();
for k = 1:size(data, 1)
    s.(sprintf('Expression%d', k)) = reshape(data(k, :, :), size(data, 2), size(data, 3));
end
save(fileName, '-struct', 's');
end
