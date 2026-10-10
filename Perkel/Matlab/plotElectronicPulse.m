function plotElectronicPulse(samplewave, samplepulse)
%PLOTELECTRONICPULSE  Figure from the "Electronic pulse representation"
%   section of PerkelTheory.nb: oscilloscope trace, IM and PM drive pulses.
figure('Name', 'Electronic pulse representation'); hold on; axis equal off;
rectangle('Position', [-10 -48 10 96], 'FaceColor', [0.85 0.85 0.85], 'EdgeColor', 'none');
rectangle('Position', [120 -48 30 96], 'FaceColor', [0.85 0.85 0.85], 'EdgeColor', 'none');
rectangle('Position', [119 -26 2 8], 'FaceColor', [0.4 0.4 0.4], 'EdgeColor', 'none');
plot([120 136], [-22 -22], 'Color', [0.4 0.4 0.4], 'LineWidth', 2);
blue = [0 105 169]/255;
n = numel(samplewave);
xs = (0:n - 1)/10 - 10;
ys = min(samplewave(:).'/15, 48);
plot(xs, ys, 'Color', blue, 'LineWidth', 2);
fill([-10, (1:n)/10 - 10, 150, -10], [0, ys, 0, 0], blue, 'FaceAlpha', .3, 'EdgeColor', 'none');
green = [0.560181, 0.691569, 0.194885];      % ColorData[97, 3]
plot([150 150], [-48 48], 'Color', green, 'LineWidth', 2);
plot([150 148; 150 148].', [-30 -30; 30 30].', 'Color', green, 'LineWidth', 2);
text(168, 0, 'PM phase', 'Rotation', -90, 'HorizontalAlignment', 'center', 'Color', green);
text(158, 30, '\pi/2', 'Color', green); text(158, -30, '0', 'Color', green);
plot([-10 21 21 51 51 150], [-30 -30 30 30 -30 -30], 'Color', green, 'LineWidth', 2);
orange = [0.9728, 0.6216, 0.0734];
plot([-10 -10], [-48 48], 'Color', orange, 'LineWidth', 2);
for y = -40:40:40, plot([-10 -8], [y y], 'Color', orange, 'LineWidth', 2); end
text(-30, 4, 'IM voltage / a.u.', 'Rotation', 90, 'HorizontalAlignment', 'center', 'Color', orange);
text(-14, 40, '1', 'HorizontalAlignment', 'right'); text(-14, 0, '0', 'HorizontalAlignment', 'right');
text(-14, -40, '-1', 'HorizontalAlignment', 'right');
plot([-10 0 0 3 3 30 30 67.5], [40 40 -40 -40 40 40 0 0], 'Color', orange, 'LineWidth', 2);
for k = 1:numel(samplepulse) - 1
    plot([64.5 + 3*k, 66 + 3*k, 66 + 3*k, 67.5 + 3*k], ...
         40*[samplepulse(k), samplepulse(k), samplepulse(k + 1), samplepulse(k + 1)], 'Color', orange, 'LineWidth', 2);
end
plot([91.5 120 120 150], [0 0 40 40], 'Color', orange, 'LineWidth', 2);
text(136, -15, '98.5 \mus', 'HorizontalAlignment', 'center'); text(135, -7, 'lock', 'HorizontalAlignment', 'center');
plot([-10 150; -10 150].', [-48 -48; 48 48].', 'k', 'LineWidth', 2);
for x = -8:2:148, plot([x x], [-48 -47], 'k'); end
for x = 0:10:140, plot([x x], [-48 -46], 'k'); end
text(70, -67, 'Time / \mus', 'HorizontalAlignment', 'center');
for k = 0:3, text(40*k, -52, num2str(.5*k), 'HorizontalAlignment', 'center'); end
end
