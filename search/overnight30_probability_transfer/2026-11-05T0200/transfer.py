"""Exact transfer of saved probability certificates to a new allowed support.

The saved inequalities remain valid if their target-graph geometry checks.
Their weighted sum and context equations may now have residual coefficients r.
For each context let lambda=max(r_i). Add (lambda-r_i) times c_i<1/9 and
subtract lambda times that context's normalization equation. Coefficients
cancel exactly; the right side changes by (lambda-sum(r_i))/9 because the
context has ten events. A negative resulting right side, or zero with positive
strict multiplier, certifies an exclusion. Renormalize all nonnegative row
multipliers to sum to one and use the existing independent exact verifier.

No LP, SDP, symmetry assumption or new certificate format is used.
"""
if not __debug__:
    raise RuntimeError('Exact transfer requires assertions; do not use -O.')
from copy import deepcopy
from fractions import Fraction as F
from pathlib import Path
import sys

ROOT=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(ROOT/'search'))
from event30_vertex_probability_lp import geometry, constraint, verify_certificate


def prepare_template(proof):
    verify_certificate(proof)
    assert 'degree_two_closed_subset' not in proof
    data=geometry(proof['graph'],proof['mean_bound_certificates'])
    return {'proof':proof,'bounds':data[3],
            'multipliers':[F(term['multiplier']) for term in proof['inequalities']],
            'equalities':list(map(F,proof['equalities']))}


def transfer(prepared,graph,base=None):
    source=prepared['proof']
    base=base or geometry(graph)
    data=base[:3]+(prepared['bounds'],set())
    rows={}
    # Reject incompatible subset geometry before computing the longer sum.
    try:
        for i,term in enumerate(source['inequalities']):
            if term['constraint']['kind']=='disjoint_subset':
                rows[i]=constraint(graph,data,term['constraint'])
    except (AssertionError,KeyError):
        return None,'incompatible_subset_geometry'
    residual=[F(0)]*30
    for i,(term,multiplier) in enumerate(zip(source['inequalities'],prepared['multipliers'])):
        row,bound,strict=rows.get(i) or constraint(graph,data,term['constraint'])
        for v,value in enumerate(row):
            if value:
                residual[v]+=multiplier*value
    equalities=prepared['equalities'][:]
    for part,multiplier in zip(graph['parts'],equalities):
        for v in part:
            residual[v]+=multiplier
    rhs=-F(source['contradiction_margin'])
    correction=[F(0)]*30
    for k,part in enumerate(graph['parts']):
        lam=max(residual[v] for v in part)
        equalities[k]-=lam
        rhs+=(lam-sum(residual[v] for v in part))/9
        for v in part:
            correction[v]=lam-residual[v]
            assert correction[v]>=0
    extra=sum(correction)
    strict_mass=F(source['strict_multiplier_mass'])+extra
    if rhs>0 or (rhs==0 and strict_mass==0):
        return None,'insufficient_margin'
    mass=1+extra
    proof=deepcopy(source)
    proof['graph']=graph
    for term in proof['inequalities']:
        term['multiplier']=str(F(term['multiplier'])/mass)
    proof['inequalities'] += [{'constraint':{'kind':'upper','orbit':v},'multiplier':str(value/mass)}
                             for v,value in enumerate(correction) if value]
    proof['equalities']=list(map(str,(value/mass for value in equalities)))
    proof['contradiction_margin']=str(-rhs/mass)
    proof['strict_multiplier_mass']=str(strict_mass/mass)
    verify_certificate(proof)
    return proof,{'added_upper_rows':sum(value>0 for value in correction),
                  'renormalization_mass':str(mass),'uncorrected_residual':list(map(str,residual)),
                  'corrected_margin_before_normalization':str(-rhs)}
