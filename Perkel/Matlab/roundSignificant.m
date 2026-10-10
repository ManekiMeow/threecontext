function y = roundSignificant(x, digits)
%ROUNDSIGNIFICANT  Round to DIGITS significant digits, as SetPrecision[x, digits] displays.
y = x;
nz = x ~= 0;
scale = 10.^(digits - 1 - floor(log10(abs(x(nz)))));
y(nz) = round(x(nz) .* scale) ./ scale;
end
