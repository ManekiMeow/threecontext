function out = parseWave(raw)
%PARSEWAVE  Port of parseWave[] from PerkelPostProcessing.nb.
%   RAW is one oscilloscope trace (first column of a scope CSV). Every 33
%   samples are reduced to (max + min); the falling edge near point 80 is
%   used as the time reference. Returns [data, lophase, circphase]:
%   the homodyne intensity, the LO phase, and the convolution phase
%   (corrected for 0.4654 times the LO phase).
raw = raw(:, 1);
wave = zeros(300, 1);
for k = 1:300
    seg = raw(1 + 33*(k - 1):33*k);
    wave(k) = max(seg) + min(seg);
end
% ListConvolve[{1, -1}, w] is diff(w); First@Position[#, Min@#] is the
% first index of the minimum, as MATLAB's min returns.
[~, iMin] = min(diff(wave(75:85)));
snap = 77 + iMin;
offset = -mean(wave(snap + 55:snap + 64));
data = wave(snap + 88) + offset;
cal = mean(offset + wave(end - 49:end));
% ArcTan[x, y] is atan2(y, x)
lophase = atan2(64 + mean(wave(snap + 24:snap + 28)), cal);
circphase = atan2(256 + 4*mean(wave(snap + 36:snap + 44)), cal) - .4654*lophase;
out = [data, lophase, circphase];
end
