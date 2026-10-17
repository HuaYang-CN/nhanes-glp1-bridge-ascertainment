"""Standalone controlled reconstruction from public XPTs. No downloads or installs."""
from pathlib import Path
import argparse,csv,json,hashlib,subprocess,sys,datetime,platform,shutil,os
sys.dont_write_bytecode=True
RELEASE=Path(__file__).resolve().parent
def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for b in iter(lambda:f.read(1048576),b''):h.update(b)
 return h.hexdigest()
def inside(p,root):return p==root or p.is_relative_to(root)
def validate(config):
 def path(key,required=True):
  value=config.get(key,'')
  if required and (not isinstance(value,str) or not value.strip()):raise ValueError('Missing/empty '+key)
  if not value:return None
  p=Path(value).expanduser()
  lexical=p if p.is_absolute() else Path(config["PROJECT_ROOT"])/p
  for ancestor in [lexical]+list(lexical.parents):
   if ancestor.exists() and (ancestor.is_symlink() or (hasattr(ancestor,"is_junction") and ancestor.is_junction())):raise ValueError("Symlink/junction path rejected before resolution: "+str(ancestor))
  if key!='PROJECT_ROOT' and not p.is_absolute():p=Path(config['PROJECT_ROOT'])/p
  return p.resolve()
 root=path('PROJECT_ROOT');config['PROJECT_ROOT']=str(root);raw=path('RAW_DATA_DIR');out=path('OUTPUT_DIR');lib=path('R_LIBRARY_DIR',False)
 if not root.is_dir() or not raw.is_dir():raise ValueError('PROJECT_ROOT and RAW_DATA_DIR must already exist')
 if inside(out,raw) or inside(raw,out) or inside(out,RELEASE) or inside(RELEASE,out) or out==root or out==Path.cwd().resolve():raise ValueError('OUTPUT_DIR overlaps inputs/code/project root/current working directory')
 if out.exists() and (not out.is_dir() or any(out.iterdir())):raise ValueError('OUTPUT_DIR must be new or empty; existing outputs are never overwritten')
 for p in [root,raw,out,RELEASE]+([lib] if lib else []):
  for parent in [p]+list(p.parents):
   if parent.exists() and (parent.is_symlink() or (hasattr(parent,'is_junction') and parent.is_junction())):raise ValueError('Symlink/junction path rejected: '+str(parent))
 if lib and not lib.is_dir():raise ValueError('Configured library does not exist')
 rscript=config.get('R_SCRIPT') or shutil.which('Rscript')
 if not rscript or not Path(rscript).is_file():raise ValueError('Supply an existing R_SCRIPT executable; no installation is performed')
 return root,raw,out,lib,Path(rscript).resolve()
def run(config_path,figures=False):
 config=json.loads(Path(config_path).read_text(encoding='utf-8-sig'));root,raw,out,lib,rscript=validate(config)
 manifest=json.loads((RELEASE/'release_manifest.json').read_text(encoding='utf-8'))
 for row in manifest['files']:
  p=RELEASE/row['path']
  if not p.is_file() or sha(p)!=row['sha256']:raise ValueError('Candidate manifest mismatch: '+row['path'])
 with (RELEASE/'schemas/public_components.csv').open(encoding='utf-8-sig',newline='') as f:components=list(csv.DictReader(f))
 inputrows=[]
 for r in components:
  p=raw/r['filename']
  if not p.is_file():raise ValueError('Missing public component: '+r['filename'])
  h=sha(p)
  if h!=r['sha256'] or str(p.stat().st_size)!=r['bytes']:raise ValueError('Public input checksum/size differs from verified release: '+r['filename'])
  inputrows.append(dict(filename=r['filename'],path=str(p),bytes=p.stat().st_size,sha256=h))
 out.mkdir(parents=True,exist_ok=True)
 for name in ['logs','qa','provenance']:(out/name).mkdir()
 (out/'provenance/input_hashes_before.json').write_text(json.dumps(inputrows,indent=2),encoding='utf-8')
 (out/'provenance/executed_release_manifest.json').write_text(json.dumps(manifest,indent=2),encoding='utf-8')
 version=dict(candidate='v0.1.0-rc1',python=sys.version,platform=platform.platform(),release_root=str(RELEASE),output_dir=str(out),config_source=str(Path(config_path).resolve()))
 (out/'provenance/execution_environment.json').write_text(json.dumps(version,indent=2),encoding='utf-8')
 resolved={'PROJECT_ROOT':str(root),'RAW_DATA_DIR':str(raw),'OUTPUT_DIR':str(out),'R_LIBRARY_DIR':str(lib) if lib else '', 'INPUT_SHA256':''}
 runtime_config=out/'provenance/runtime_config.csv'
 def saveconfig():
  with runtime_config.open('w',encoding='utf-8',newline='') as f:
   w=csv.writer(f);w.writerow(['key','value']);w.writerows(resolved.items())
 saveconfig();events=[];started=datetime.datetime.now(datetime.timezone.utc).isoformat();status=0
 rmodules=['input_preflight.R','b1/01_build_analysis_dataset.R','b1/02_reproduce_primary_results_v4.R','b1/03_sensitivity_fasting_glucose.R','b1/04_ascertainment_gap_conditional.R','membership_bridge.R','b2/01_design.R','b2/02_diq.R','b2/03_denominators.R','b2/04_age.R','b2/05_glycemic.R','b2/06_paired.R','b2/07_subgroups.R','a3/01_preflight_and_design.R','a3/02_generate_table1.R','a3/03_generate_completeness.R','a3/04_validate_outputs.R']
 def invoke(name,command,script):
  tag=name.replace('/','_').replace('.','_');start=datetime.datetime.now(datetime.timezone.utc).isoformat();stdout=out/'logs'/(tag+'_stdout.log');stderr=out/'logs'/(tag+'_stderr.log')
  with stdout.open('wb') as so,stderr.open('wb') as se:
   child_env=os.environ.copy()
   if Path(command[0]).resolve()==rscript:
    child_env.update(LANG="English_United States.utf8",LC_ALL="English_United States.utf8")
   result=subprocess.run(command,cwd=root,stdout=so,stderr=se,shell=False,env=child_env)
  row=dict(module=name,command=command,script_sha256=sha(script),start=start,end=datetime.datetime.now(datetime.timezone.utc).isoformat(),exit_code=result.returncode,stdout=str(stdout),stderr=str(stderr),output_root=str(out),status='PASS' if result.returncode==0 else 'FAIL')
  events.append(row);(out/'logs/execution_events.json').write_text(json.dumps(events,indent=2),encoding='utf-8');print(name,'exit',result.returncode,flush=True)
  if result.returncode:raise subprocess.CalledProcessError(result.returncode,command)
 try:
  for module in rmodules:
   invoke(module,[str(rscript),'--vanilla',str(RELEASE/'src/run_module.R'),module,str(runtime_config),str(RELEASE)],RELEASE/'src'/module)
   if module=='b1/01_build_analysis_dataset.R':resolved['INPUT_SHA256']=sha(out/'b1/diagnostics/rebuilt/bridge_analysis_ready.rds');saveconfig()
  invoke('registry_adapter',[sys.executable,'-B',str(RELEASE/'src/registry_adapter.py'),str(RELEASE),str(out)],RELEASE/'src/registry_adapter.py')
  if figures:invoke('presentation/03_draw_figures.R',[str(rscript),'--vanilla',str(RELEASE/'src/run_module.R'),'presentation/03_draw_figures.R',str(runtime_config),str(RELEASE)],RELEASE/'src/presentation/03_draw_figures.R')
  invoke('concordance',[sys.executable,'-B',str(RELEASE/'tests/concordance.py'),str(RELEASE),str(out)],RELEASE/'tests/concordance.py')
 except subprocess.CalledProcessError as e:status=e.returncode
 finally:
  after=[dict(r,sha256=sha(raw/r['filename']),bytes=(raw/r['filename']).stat().st_size) for r in inputrows]
  inputs_unchanged=all(a['sha256']==b['sha256'] and a['bytes']==b['bytes'] for a,b in zip(inputrows,after))
  (out/'provenance/input_hashes_after.json').write_text(json.dumps(after,indent=2),encoding='utf-8')
  if not inputs_unchanged:status=3
  summary=dict(start=started,end=datetime.datetime.now(datetime.timezone.utc).isoformat(),exit_code=status,modules_executed=len(events),input_hashes_unchanged=inputs_unchanged,scientific_test='PASS' if status==0 else 'FAIL',optional_figures_requested=figures,native_R1_Word_layout='NOT_TESTED')
  (out/'qa/run_summary.json').write_text(json.dumps(summary,indent=2),encoding='utf-8');print(json.dumps(summary),flush=True)
 return status
if __name__=='__main__':
 parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--config',required=True);parser.add_argument('--figures',action='store_true');args=parser.parse_args()
 try:raise SystemExit(run(args.config,args.figures))
 except (ValueError,FileNotFoundError) as e:print('PRE_EXECUTION_FAILURE:',str(e),file=sys.stderr);raise SystemExit(4)
