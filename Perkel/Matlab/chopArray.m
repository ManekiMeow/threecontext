function x = chopArray(x, delta)
%CHOPARRAY  Port of Mathematica's Chop: set entries with |x| < delta to 0.
if nargin < 2, delta = 1e-10; end
x(abs(x) < delta) = 0;
end
