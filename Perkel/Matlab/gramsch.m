function raypool = gramsch(vecpool)
%GRAMSCH  Port of gramsch[] from PerkelTheory.nb.
%   Appends e_20 ... e_37 to the 19 orthonormal rays in VECPOOL, each
%   orthogonalised against every ray already in the pool (classical
%   Gram-Schmidt) and normalised, giving an orthonormal basis of R^37.
raypool = vecpool;
for k = 1:18
    e = zeros(1, 37); e(19 + k) = 1;
    vectemp = e - (raypool * e.').' * raypool;
    raypool(end + 1, :) = vectemp / norm(vectemp); %#ok<AGROW>
end
end
