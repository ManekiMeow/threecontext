function writeTable(fileName, M)
%WRITETABLE  Export[file, M, "Table"]: tab-separated rows, 5 significant digits
%   (the notebook applies SetPrecision[#, 5] before exporting).
fid = fopen(fileName, 'w');
fmt = [repmat('%.5g\t', 1, size(M, 2) - 1), '%.5g\n'];
fprintf(fid, fmt, M.');
fclose(fid);
end
