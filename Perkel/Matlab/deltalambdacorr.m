function c = deltalambdacorr(m, n, deltalambda, deltalambda0, l0)
%DELTALAMBDACORR  Port of deltalambdacorr[m, n] from PerkelTheory.nb.
%   Normalised overlap of two Gaussian spectra of widths deltalambda(m) and
%   deltalambda(n), the second shifted by n*deltalambda0*2e-5/l0 (20 um
%   misalignment). The integration variable is x/deltalambda0; that scale
%   factor cancels in the ratio.
wm = deltalambda(m) / deltalambda0;
wn = deltalambda(n) / deltalambda0;
shift = n * (2e-5) / l0;
num = integral(@(u) exp(-u.^2 / (2*wm^2/log(2))) .* exp(-(u - shift).^2 / (2*wn^2/log(2))), -Inf, Inf);
dm = integral(@(u) exp(-u.^2 / (wm^2/log(2))), -Inf, Inf);
dn = integral(@(u) exp(-u.^2 / (wn^2/log(2))), -Inf, Inf);
c = num / (sqrt(dm) * sqrt(dn));
end
