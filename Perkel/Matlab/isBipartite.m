function tf = isBipartite(A)
%ISBIPARTITE  True if the graph with adjacency matrix A has no odd cycle.
%   Used by PerkelTheory.m to confirm chi(Perkel) > 2.
n = size(A, 1); side = zeros(n, 1);
for s = 1:n
    if side(s), continue; end
    side(s) = 1; queue = s;
    while ~isempty(queue)
        v = queue(1); queue(1) = [];
        for w = find(A(v, :))
            if side(w) == 0
                side(w) = -side(v); queue(end + 1) = w; %#ok<AGROW>
            elseif side(w) == side(v)
                tf = false; return;
            end
        end
    end
end
tf = true;
end
