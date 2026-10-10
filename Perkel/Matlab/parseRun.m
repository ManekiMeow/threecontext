function out = parseRun(traceFile, readTrace, sub, run, first, nStates)
%PARSERUN  Parse the ten repetitions of states first+1 ... first+nStates of
%   one measurement run with parseWave. OUT is nStates-by-10-by-3, i.e.
%   Map[parseWave, Table[Import[...][[;;, 1]], {i, nStates}, {j, 10}], {2}].
out = zeros(nStates, 10, 3);
for i = 1:nStates
    for j = 1:10
        raw = readTrace(traceFile(sub, run, i + first, j));
        out(i, j, :) = parseWave(raw(:, 1));
    end
end
end
