function rgb = blendColor(stops, t)
%BLENDCOLOR  Port of Mathematica's Blend[{{t1, col1}, {t2, col2}, ...}, t].
%   STOPS is an n-by-4 matrix [t r g b] with increasing t. T may be an array;
%   RGB is numel(T)-by-3. Values outside the stop range are clamped.
t = min(max(t(:), stops(1, 1)), stops(end, 1));
rgb = [interp1(stops(:, 1), stops(:, 2), t), ...
       interp1(stops(:, 1), stops(:, 3), t), ...
       interp1(stops(:, 1), stops(:, 4), t)];
end
