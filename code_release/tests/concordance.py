"""Independent reference comparisons. Fixtures are never calculation inputs."""
from pathlib import Path
import csv,json,sys,math
from decimal import Decimal
def read(p):return list(csv.DictReader(Path(p).open(encoding='utf-8-sig',newline='')))
def verify(release,out):
 release=Path(release);out=Path(out);tol=json.loads((release/'tests/fixtures/tolerances.json').read_text())
 actual={r['result_id']:r for r in read(out/'registry/estimands.csv')+read(out/'registry/table1.csv')+read(out/'registry/completeness.csv')}
 checks=[]
 for cell in json.loads((release/"schemas/table_layout.json").read_text(encoding="utf-8")):
  if cell["result_id"]=="SECTION":continue
  ok=cell["result_id"] in actual
  checks.append(dict(reference_file="table cell source coverage",result_id=cell["result_id"],field=cell["table"]+":"+cell["row"]+":"+cell["column"],expected_full_precision="Fresh result ID present",calculated_full_precision="PRESENT" if ok else "MISSING",absolute_difference="",prespecified_tolerance=0,status="PASS" if ok else "FAIL"))
 for filename in ['locked_numerical_reference.csv','table1_completeness_reference.csv']:
  for r in read(release/'tests/fixtures'/filename):
   rid=r['result_id'];field=r['field'];expected=r['expected'];got=actual.get(rid,{}).get(field,'MISSING_FIELD')
   count=field.endswith('_n') or field in ['n','df','ci_df','represented_domain_df','parent_df','universe_n']
   threshold=0.0 if count else tol.get(field,tol['other_numeric']);difference=''
   if expected in ['NA','NaN','']:
    ok=got in ['NA','NaN','']
   else:
    try:
     difference=abs(Decimal(got)-Decimal(expected));ok=difference<=Decimal(str(threshold))
    except:ok=got==expected
   checks.append(dict(reference_file=filename,result_id=rid,field=field,expected_full_precision=expected,calculated_full_precision=got,absolute_difference=str(difference),prespecified_tolerance=threshold,status='PASS' if ok else 'FAIL',reference_source_sha256=r.get('reference_source_sha256','')))
 # Explicitly preserve non-estimability and classification/precision policy.
 for r in read(out/'b2/results/subgroup_reliability.csv'):
  if r['implementation']=='B0_full_parent_domain_subgroupdf' and r['domain_df']=='0':
   ok=r['status']=='NOT_ESTIMABLE' and r['lower']=='NA' and r['upper']=='NA'
   checks.append(dict(reference_file='semantic identities',result_id=r['characteristic'],field='df0 CI status',expected_full_precision='NOT_ESTIMABLE with NA limits',calculated_full_precision=f"{r['status']}; {r['lower']}; {r['upper']}",absolute_difference='',prespecified_tolerance=0,status='PASS' if ok else 'FAIL'))
 for file,field in [('paired_covariance_and_identity_checks.csv','identity_pass'),('glycemic_flag_verification.csv','matches'),('DIQ_raw_input_value_integrity.csv','raw_aligned_values_and_missingness_exact')]:
  for r in read(out/'b2/diagnostics'/file):
   checks.append(dict(reference_file=file,result_id=r.get('contrast',r.get('flag','raw DIQ160')),field=field,expected_full_precision='TRUE',calculated_full_precision=r[field],absolute_difference='',prespecified_tolerance=0,status='PASS' if r[field]=='TRUE' else 'FAIL'))
 p=out/'qa/numerical_concordance.csv';p.parent.mkdir(exist_ok=True)
 fields=list(dict.fromkeys(k for r in checks for k in r))
 with p.open('w',encoding='utf-8',newline='') as f:
  w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(checks)
 failures=[r for r in checks if r['status']=='FAIL'];summary=dict(comparisons=len(checks),passed=len(checks)-len(failures),failed=len(failures),tolerances_frozen=True)
 (out/'qa/concordance_summary.json').write_text(json.dumps(summary,indent=2),encoding='utf-8')
 print(json.dumps(summary))
 for r in failures[:30]:print(json.dumps(r,ensure_ascii=False))
 if failures:raise SystemExit(2)
if __name__=='__main__':verify(sys.argv[1],sys.argv[2])
