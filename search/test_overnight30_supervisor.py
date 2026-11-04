"""Network-free lifecycle and publication regression checks for the supervisor.

Run: python search/test_overnight30_supervisor.py

Uses temporary mock workers and a fake Git adapter; no repository or remote is
modified. Includes deadline/forced termination, locking, restart, compression,
startup-failure protection, and push-before-next-commit checks. Takes about 20s.
"""
import concurrent.futures
from datetime import datetime, timedelta, timezone
import gzip
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import time

SCRIPT = Path(__file__).with_name('overnight30_supervisor.py')
TEMPORARY = tempfile.TemporaryDirectory(prefix='ghz30-supervisor-test-')
BASE = Path(TEMPORARY.name)
WORKER = '''import argparse,json,os,signal,subprocess,sys,time
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('--seconds',type=float);p.add_argument('--seed',type=int);p.add_argument('--output-dir',type=Path);p.add_argument('--family');p.add_argument('--resume-from',type=Path);a=p.parse_args()
kind='symmetry' if a.family else 'asymmetric'
d=a.output_dir/a.family/f'seed_{a.seed}' if a.family else a.output_dir
d.mkdir(parents=True,exist_ok=True)
if a.resume_from: assert a.resume_from.exists(),a.resume_from
state={'seed':a.seed,'resumed':bool(a.resume_from)}
name='state.json' if a.family else f'asymmetric_seed{a.seed}.json'
(d/name).write_text(json.dumps(state));(d/'models.jsonl').write_text(json.dumps(state)+'\\n')
mode=Path(__file__).parent.parent.name
if mode=='failed': raise SystemExit(2)
if mode=='orphan':
 child=subprocess.Popen([sys.executable,'-c','import signal,time; signal.signal(signal.SIGTERM,signal.SIG_IGN); time.sleep(60)'])
 (d/'child.pid').write_text(str(child.pid))
if mode=='hanging':
 signal.signal(signal.SIGTERM,signal.SIG_IGN)
 time.sleep(60)
else: time.sleep(min(a.seconds,0.5))
'''
for mode in ('normal','failed','hanging','orphan'):
    d=BASE/mode/'search';d.mkdir(parents=True)
    for worker in ('symmetry','asymmetric'):
        (d/f'overnight30_{worker}.py').write_text(WORKER)

def command(mode, repo, seconds=8):
    until=(datetime.now(timezone.utc)+timedelta(seconds=seconds)).isoformat()
    return [sys.executable,str(SCRIPT),'--until',until,'--data-repo',str(repo),'--code-root',str(BASE/mode),'--workers','3','--chunk-seconds','1','--first-chunk-seconds','0.5','--families','c2_pairs,c3_fixed','--no-publish']

def check(mode):
    repo=BASE/f'data_{mode}';started=time.monotonic()
    process=subprocess.Popen(command(mode,repo),stdout=subprocess.PIPE,stderr=subprocess.PIPE,text=True)
    root=repo/'search/overnight30_runs/2026-11-04'
    if mode=='normal':
        while not (root/'status.json').exists(): time.sleep(.05)
        duplicate=subprocess.run(command(mode,repo),capture_output=True,text=True,timeout=5)
        assert duplicate.returncode and 'holds this run' in duplicate.stderr,duplicate
    out,err=process.communicate(timeout=20)
    assert process.returncode==0,(mode,out,err)
    status=json.loads((root/'status.json').read_text());manifest=json.loads((root/'manifest.json').read_text())
    assert status['stopped'] and not status['active'],status
    expected='repeated_worker_startup_failure' if mode=='failed' else 'deadline'
    assert status['stop_reason']==expected,status
    assert time.monotonic()-started<13,(mode,time.monotonic()-started)
    for path in root.rglob('*.jsonl.gz'):
        for line in gzip.open(path,'rt'): json.loads(line)
    assert not list(root.rglob('*.jsonl'))
    for path in root.rglob('child.pid'):
        status_path=Path('/proc')/path.read_text()/'status'
        if status_path.exists():
            state=next(line for line in status_path.read_text().splitlines() if line.startswith('State:'))
            assert 'Z (zombie)' in state,state
    if mode=='normal':
        old=len(manifest['jobs'])
        resumed=subprocess.run(command(mode,repo),capture_output=True,text=True,timeout=20)
        assert resumed.returncode==0,(resumed.stdout,resumed.stderr)
        jobs=json.loads((root/'manifest.json').read_text())['jobs']
        assert len(jobs)>old and len({j['id'] for j in jobs})==len(jobs)
        assert any(j['resume_job'] is not None for j in jobs[old:]),jobs
    return mode, status['jobs_finished'], status['stop_reason']
with concurrent.futures.ThreadPoolExecutor(max_workers=4) as executor:
    results=list(executor.map(check,('normal','failed','hanging','orphan')))
print(json.dumps({'lifecycle_tests':results}))

import argparse,importlib.util,json,tempfile
from pathlib import Path
import subprocess,sys
spec=importlib.util.spec_from_file_location('overnight',str(SCRIPT))
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
class Fake(m.Supervisor):
 def git(self,*arguments,timeout=30):
  self.calls=getattr(self,'calls',[]);self.calls.append(arguments)
  text='';code=0
  if arguments==('rev-parse','--show-toplevel'):text=str(self.repo)
  elif arguments==('symbolic-ref','--short','HEAD'):text=m.BRANCH
  elif arguments[:2]==('remote','get-url'):text='https://github.com/ManekiMeow/threecontext.git'
  elif arguments==('diff','--cached','--quiet'):code=1
  elif arguments==('rev-parse','HEAD'):text='testsha'
  elif arguments[0]=='push':code=1 if getattr(self,'fail_push',False) else 0
  elif arguments[0]=='commit':self.fail_push=True
  return subprocess.CompletedProcess(arguments,code,text)
a=argparse.Namespace(data_repo=BASE/'git_mock',code_root=SCRIPT.parent.parent,run_date='test',until='2026-11-05T05:00:00Z',no_publish=False,branch=m.BRANCH,stale_host_seconds=1200)
s=Fake(a);s.queue_directory(str(s.relative/'chunks/00000'));assert not s.publish();assert s.control['push_pending']
commits=sum(c[0]=='commit' for c in s.calls);assert commits==1
assert not s.publish();assert sum(c[0]=='commit' for c in s.calls)==commits
ci=next(i for i,c in enumerate(s.calls) if c[0]=='commit');assert s.calls[ci+1][0]=='push'
s.fail_push=False;assert s.publish();assert not s.control['push_pending']
print(json.dumps({'push_failure_preserves_commit':True,'no_second_commit_before_successful_push':True,'immediate_push_after_commit':True}))

TEMPORARY.cleanup()
