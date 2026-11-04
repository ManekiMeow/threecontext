import sys,pathlib,json,gzip,itertools,subprocess,time,collections,hashlib
sys.path.insert(0,str(pathlib.Path(__file__).resolve().parents[2]))
from event30_probability_lp import certify,verify_certificate
from overnight30_screen import fast_screen,normalize_graph
import argparse
parser=argparse.ArgumentParser(description='Reproduce the bounded published C5 quotient review.')
parser.add_argument('--data-repo',type=pathlib.Path,default=pathlib.Path('/workspace/ghz30-overnight'))
parser.add_argument('--seconds',type=float,default=45)
parser.add_argument('--limit',type=int,default=150)
args=parser.parse_args()
root=args.data_repo
commit='f37683f803842de5b3a8d7fcffcbdd781bcff492'
files=subprocess.check_output(['git','ls-tree','-r','--name-only',commit],cwd=root,text=True).splitlines()
files=[p for p in files if '_c5_double/' in p and p.endswith(('models.jsonl','models.jsonl.gz'))]
perms=[]
for order in itertools.permutations(range(3)):
 for flips in itertools.product(range(2),repeat=3):
  perms.append([2*k+(a^flips[k]) for k in order for a in range(2)])
cache={}; pool={}; scored=0
for path in files:
 with (gzip.open(root/path,'rt') if path.endswith('.gz') else (root/path).open()) as stream:
  for line in stream:
   r=json.loads(line)
   if r.get('status')!='scored_numerically':continue
   scored+=1
   data=r['graph6'].encode(); q=[[0]*6 for _ in range(6)]
   for v in range(1,30):
    for u in range(v):
     bit=v*(v-1)//2+u
     if (data[1+bit//6]-63)&(1<<(5-bit%6)):
      q[u//5][v//5]+=1;q[v//5][u//5]+=1
   assert all(x%5==0 for row in q for x in row)
   q=[[x//5 for x in row]for row in q]; raw=tuple(sum(q,[]))
   if raw not in cache:cache[raw]=min(tuple(q[p[i]][p[j]] for i in range(6) for j in range(6))for p in perms)
   key=cache[raw]
   if key not in pool or r['theta_numeric']>pool[key]['theta_numeric']:pool[key]={**r,'source':path,'quotient':q}
print(json.dumps({'files':len(files),'scored':scored,'quotient_classes':len(pool)}),flush=True)
out=pathlib.Path(__file__).resolve().parent;out.mkdir(parents=True,exist_ok=True)
start=time.monotonic();results=[]
for index,r in enumerate(sorted(pool.values(),key=lambda r:-r['theta_numeric'])[:args.limit]):
 if time.monotonic()-start>args.seconds:break
 cert=certify(r['quotient'])
 screen=fast_screen(r)
 item={'theta':r['theta_numeric'],'graph6':r['graph6'],'source':r['source'],'quotient':r['quotient'],'old_screen':screen.get('certificate',{}).get('method'),'exact_lp_excluded':cert is not None}
 if cert:
  digest=hashlib.sha256(json.dumps(cert['quotient']).encode()).hexdigest()[:20]
  name=digest+'.json';(out/name).write_text(json.dumps(cert,indent=2)+'\n');item['certificate']=name;item['margin']=str(verify_certificate(cert))
 results.append(item)
 if cert and screen['status']!='excluded':print(json.dumps({'index':index,'theta':r['theta_numeric'],'margin':item['margin'],'new_exclusion':True}),flush=True)
summary={'published_commit':commit,'source_files':files,'scored_records':scored,'quotient_classes_in_scored_pool':len(pool),'reviewed':len(results),'certified':sum(r['exact_lp_excluded']for r in results),'new_beyond_fast_screen':sum(r['exact_lp_excluded'] and r['old_screen']is None for r in results),'scope':'Selected quotient classes only. A certificate excludes every support admitting the stated actual group-orbit quotient; no global nonexistence claim.','results':results}
(out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
print(json.dumps({k:v for k,v in summary.items() if k not in ('source_files','results')}),flush=True)
