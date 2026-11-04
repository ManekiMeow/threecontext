"""Small independent checker for the positive-hint (RUP) subset of text LRAT.

This deliberately rejects RAT steps. Each added clause is verified by unit
propagation under its negation, using only the named earlier clauses. The
checker must derive an empty clause. CNF inputs still require an independent
audit that they encode the claimed mathematical problem and certified cuts.
No SAT solver, proof generator, external checker or third-party module is used.
"""
from pathlib import Path
import argparse
import json

if not __debug__:
    raise RuntimeError('LRAT verification requires assertions; run Python without -O')


def read_cnf(path):
    clauses, pending, variables, count = {}, [], None, None
    for line in Path(path).read_text().splitlines():
        fields = line.split()
        if not fields or fields[0] == 'c':
            continue
        if fields[0] == 'p':
            assert variables is None and len(fields) == 4 and fields[1] == 'cnf'
            variables, count = map(int, fields[2:])
            assert variables >= 0 and count >= 0
            continue
        assert variables is not None
        for token in fields:
            literal = int(token)
            if literal:
                assert abs(literal) <= variables
                pending.append(literal)
            else:
                clauses[len(clauses)+1] = frozenset(pending)
                pending = []
    assert variables is not None and not pending and len(clauses) == count
    return variables, clauses


def rup(clause, hints, clauses):
    assigned = {-literal for literal in clause}
    if any(-literal in assigned for literal in assigned):
        return True  # The proposed clause itself is tautological.
    for identifier in hints:
        assert identifier > 0, 'RAT hints are outside this checker\'s supported fragment'
        antecedent = clauses[identifier]
        if antecedent & assigned:
            raise ValueError('LRAT hint names a satisfied, non-unit antecedent')
        remaining = antecedent - {-literal for literal in assigned}
        if not remaining:
            return True
        if len(remaining) != 1:
            raise ValueError('LRAT hint is neither unit nor conflicting')
        assigned.update(remaining)
    return False


def verify(cnf_path, proof_path):
    variables, clauses = read_cnf(cnf_path)
    original = len(clauses)
    empty_derived = any(not clause for clause in clauses.values())
    additions = deletions = 0
    seen = set(clauses)
    for line in Path(proof_path).read_text().splitlines():
        fields = line.split()
        if not fields or fields[0] == 'c':
            continue
        identifier = int(fields[0])
        assert identifier > 0
        if fields[1] == 'd':
            removed = list(map(int, fields[2:]))
            assert removed and removed[-1] == 0 and all(i > 0 for i in removed[:-1])
            for i in removed[:-1]:
                clauses.pop(i, None)
                deletions += 1
            continue
        values = list(map(int, fields[1:]))
        zero = values.index(0)
        clause, hints = frozenset(values[:zero]), values[zero+1:]
        assert identifier not in seen and all(0 < abs(literal) <= variables for literal in clause)
        assert hints and hints[-1] == 0 and all(i > 0 for i in hints[:-1])
        if not rup(clause, hints[:-1], clauses):
            raise ValueError('LRAT hints did not derive a conflict')
        clauses[identifier] = clause
        seen.add(identifier)
        additions += 1
        empty_derived |= not clause
    assert empty_derived, 'No empty clause was derived'
    return {'variables': variables, 'input_clauses': original,
            'verified_rup_additions': additions, 'deletions': deletions,
            'empty_clause_verified': True}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('cnf', type=Path)
    parser.add_argument('proof', type=Path)
    args = parser.parse_args()
    print(json.dumps(verify(args.cnf, args.proof), indent=2))
