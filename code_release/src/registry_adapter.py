"""Map freshly calculated outputs into aggregate registries. No test fixtures imported."""
from pathlib import Path
import csv,json,hashlib,sys
from decimal import Decimal,ROUND_HALF_UP

def read(p):return list(csv.DictReader(Path(p).open(encoding='utf-8-sig',newline='')))
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def write(p,rows):
 p=Path(p)
 if p.exists():raise FileExistsError(p)
 p.parent.mkdir(parents=True,exist_ok=True)
 fields=list(dict.fromkeys(k for r in rows for k in r))
 with p.open('w',encoding='utf-8',newline='') as f:
  w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(rows)
def fmt(v,d=1):
 if v in ['',None,'NA','NaN']:return 'Not estimable'
 return format(Decimal(v).quantize(Decimal(1).scaleb(-d),rounding=ROUND_HALF_UP),f'.{d}f')
def make(release,out):
 release=Path(release).resolve();out=Path(out).resolve();cache={}
 def data(path):
  if path not in cache:cache[path]=read(out/path)
  return cache[path]
 mapping=json.loads((release/'schemas/registry_mapping.json').read_text(encoding='utf-8'));registry=[]
 for m in mapping:
  source=out/m['source_file'];rows=data(m['source_file']);r=rows[m['source_row']-1]
  if not all(r.get(k)==str(v) for k,v in m['selector'].items()):raise ValueError(f"Source row/selector mismatch: {m['result_id']}")
  z={k:m.get(k,'') for k in ['result_id','domain','estimand','definition','unit','weight','inference_scope']}
  keys=['estimate','se','lower','upper','n','universe_n','positive_n','comparison_n','conditional_domain_n','missing_n','ci_df','parent_df','represented_domain_df','p_value','p_value_Holm']
  for k in keys:
   sk=m['field_map'].get(k,k)
   if k=='p_value_Holm' and sk not in r:sk='p_value_Holm_two_contrasts'
   z[k]=r.get(sk,'')
  if r.get('estimand')=='conditional_measurement_only':z['conditional_domain_n']=r['n']
  if z['domain']=='PAIR':z['universe_n']=r['n']
  z['missing_n']=z['missing_n'] or r.get('missing_outcome_n','')
  z['represented_domain_df']=z['represented_domain_df'] or r.get('domain_df',r.get('design_df',''))
  z['ci_df']=z['ci_df'] or r.get('design_df','')
  z['ci_method']=r.get('ci_method') or m['ci_method']
  z['weight']=r.get('weight_variable') or m['weight']
  z['df_policy']=r.get('df_policy') or m['df_policy']
  z['analytical_status']=r.get('status','FAITHFULLY_REPRODUCED')
  z['reliability_flags']=r.get('reliability_flags','')
  z['effective_nonmissing_n']=z['n'];z['source_path']=m['source_file'];z['source_row']=str(m['source_row']);z['source_sha256']=sha(source)
  z['source_selector']=json.dumps(m['selector'],ensure_ascii=False);z['field_map']=json.dumps(m['field_map']);z['source_record_json']=json.dumps(r,ensure_ascii=False)
  if z['domain']=='TRACKA':
   z['universe_n']=z['n'];z['ci_df']=z['ci_df'] or str(int(data('b2/results/design_comparison.csv')[0]['parent_df']))
   if m['source_file'].endswith('KEY_measurement_only_among_measurement_positive.csv'):
    q=next(x for x in data('b2/results/diq_scenario_estimates.csv') if x['scenario']=='D0' and x['scope']=='scenario_specific' and x['estimand']=='conditional_measurement_only')
    z.update(n=q['n'],universe_n=q['comparison_n'],conditional_domain_n=q['n'],positive_n=q['positive_n'])
   elif m['source_file'].endswith('fasting_glucose_conditional_gain.csv'):
    q=next(x for x in data('b2/results/glycemic_sensitivity.csv') if x['scenario']=='G0_fasting' and x['scope']=='scenario_specific' and x['estimand']=='HbA1c_or_FPG')
    gain=next(x for x in data('b2/results/glycemic_sensitivity.csv') if x['scenario']=='G0_fasting' and x['scope']=='scenario_specific' and x['estimand']=='gained_by_FPG')
    z.update(n=q['positive_n'],universe_n=q['n'],positive_n=gain['positive_n'])
  z['display_estimate']=fmt(z['estimate']);z['display_ci']=f"{fmt(z['lower'])} to {fmt(z['upper'])}" if z['lower'] not in ['',None,'NA'] else 'Not estimable'
  registry.append(z)
 if len(registry)!=len({r['result_id'] for r in registry}):raise ValueError('Duplicate result ID')
 write(out/'registry/estimands.csv',registry)
 t1=read(out/'a3/diagnostics/statistics_full_precision.csv');comp=read(out/'a3/diagnostics/completeness_full_precision.csv')
 for r in t1:
  d=2 if r['variable']=='LBXSCR' else 1;ci=f"{fmt(r['lower_ci'],d)} to {fmt(r['upper_ci'],d)}" if r['status']=='ESTIMATED' else 'Not estimable'
  r['formatted_value']=f"{fmt(r['estimate'],d)} ({ci})" if r['statistic']=='weighted_mean' else f"{int(r['category_n'])}; {fmt(r['estimate'])}% ({ci})"
  r['source_artifact']='a3/diagnostics/statistics_full_precision.csv';r['verification_status']='COMPUTED_WITH_INDEPENDENT_ESTIMATOR_CHECKS'
 for r in comp:
  r['formatted_missing_percent']=fmt(r['unweighted_missing_percent']);r['source_artifact']='a3/diagnostics/completeness_full_precision.csv';r['verification_status']='COMPUTED_WITH_INDEPENDENT_COUNTS'
 write(out/'registry/table1.csv',t1);write(out/'registry/completeness.csv',comp)
 idx={r['result_id']:r for r in registry};tables={
  'Table1_Source.csv':t1,'Table_S1_Source.csv':comp,
  'Table_S2_Source.csv':[idx[f'C02-DENOM-{k:03d}'] for k in range(18,25)],
  'Table_S3_Source.csv':[r for r in registry if r['domain']=='DIQ'],
  'Table_S4_Source.csv':[idx[f'C02-AGE-{k:03d}'] for k in list(range(1,10))+list(range(19,28))]+[idx[f'C02-GLYC-{k:03d}'] for k in [1,2,3,4,9,10,11,12,13,14]],
  'Table_S5_Source.csv':[idx[f'C02-SUB-{k:03d}'] for k in range(3,46,3)],
  'Table_S6_Source.csv':[idx['C02-PAIR-001'],idx['C02-PAIR-002']],
  'Table_S7_Source.csv':[idx['C02-DENOM-025'],idx['C02-DENOM-026']]}
 for name,rows in tables.items():write(out/'source_tables'/name,rows)
 panels=[];flow=data('b2/results/denominator_audit.csv')
 for m in json.loads((release/'schemas/panel_mapping.json').read_text(encoding='utf-8')):
  if m['result_id'].startswith('C03B-SUPPORT-FLOW-'):
   i=int(m['result_id'].rsplit('-',1)[1])-1;r=flow[i]
   if r['record_type']!='selection_flow':raise ValueError('Flow row mapping changed')
   z=dict(estimate=r['nonmissing_n'],n=r['nonmissing_n'],lower='',upper='',universe_n=r['universe_n'],missing_n=r['population_removed_n'],positive_n='',comparison_n='',conditional_domain_n='',unit='count',weight='None (unweighted selection count)',source_path='b2/results/denominator_audit.csv',source_row=str(i+1))
  else:z=idx[m['result_id']].copy()
  z.update(m);z['displayed_estimate']=fmt(z['estimate']);z['displayed_ci']=z.get('display_ci','Not applicable');panels.append(z)
 write(out/'presentation/results/panel_plot_data.csv',panels)
 write(out/'source_tables/Figure_Panel_Source.csv',panels)
 for d in ['presentation/figures/main/panels','presentation/figures/supplementary/panels','presentation/logs']:(out/d).mkdir(parents=True,exist_ok=True)
 print(f'Fresh registry: {len(registry)} estimands, {len(t1)} Table1 rows, {len(comp)} completeness rows, {len(panels)} figure-source rows, S1-S7 source CSVs')

if __name__=='__main__':make(sys.argv[1],sys.argv[2])
