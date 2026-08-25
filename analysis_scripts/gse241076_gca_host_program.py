from pathlib import Path
import os
import gzip, itertools, json
import numpy as np
import pandas as pd
from scipy.stats import ttest_ind

ROOT=Path(os.environ.get("MICROBIOME_GATE_ROOT", Path(__file__).resolve().parents[1]))
RAW=ROOT/"01_raw_data"/"crc_axis_external"/"GSE241076"/"GSE241076_counts.txt.gz"
OUT=ROOT/"06_results"/"crc_axis_external_yachida2019"/"host_program_GSE241076"; OUT.mkdir(parents=True,exist_ok=True)
CONTROL=["con_1","con_2","con_3"]; GCA=["fa_1","fa_2","fa_3"]
LOCKED_PROGRAM=["SOX14","ZDHHC9","CD274"]
DESCRIPTIVE=["NR1H4","SOX14","ZDHHC9","CD274"]

def bh(p):
 p=np.asarray(p,float); n=len(p); o=np.argsort(p); z=p[o]*n/np.arange(1,n+1); z=np.minimum.accumulate(z[::-1])[::-1]; q=np.empty(n); q[o]=np.minimum(z,1); return q

with gzip.open(RAW,"rt") as f: counts=pd.read_csv(f,sep="\t").set_index("gene_id")
counts=counts[CONTROL+GCA].apply(pd.to_numeric,errors="coerce").fillna(0)
lib=counts.sum(axis=0); cpm=counts.div(lib,axis=1)*1e6; logcpm=np.log2(cpm+.5)
keep=(counts>=10).sum(axis=1)>=3
rows=[]
for gene,x in logcpm.loc[keep].iterrows():
 a=x[CONTROL].to_numpy(float); b=x[GCA].to_numpy(float); stat,p=ttest_ind(b,a,equal_var=False)
 rows.append(dict(gene=gene,mean_logCPM_control=a.mean(),mean_logCPM_GCA=b.mean(),log2FC_GCA_vs_control=b.mean()-a.mean(),welch_t=stat,p_welch=p))
de=pd.DataFrame(rows); de["q_bh"]=bh(de.p_welch.fillna(1)); de=de.sort_values("p_welch"); de.to_csv(OUT/"gene_level_logCPM_welch.csv",index=False)

Z=logcpm.loc[[g for g in LOCKED_PROGRAM if g in logcpm.index]].T
Z=(Z-Z.mean(axis=0))/Z.std(axis=0,ddof=1)
score=Z.mean(axis=1)
obs=float(score[GCA].mean()-score[CONTROL].mean())
samples=CONTROL+GCA; diffs=[]
for gca_idx in itertools.combinations(range(6),3):
 g=[samples[i] for i in gca_idx]; c=[s for s in samples if s not in g]
 diffs.append(float(score[g].mean()-score[c].mean()))
perm_p=float(sum(abs(d)>=abs(obs)-1e-12 for d in diffs)/len(diffs))
pd.DataFrame({"sample":samples,"condition":["control"]*3+["GCA"]*3,"locked_program_score":[score[s] for s in samples]}).to_csv(OUT/"locked_program_sample_scores.csv",index=False)

gene_summary=[]
for g in DESCRIPTIVE:
 if g in counts.index:
  rec=de.loc[de.gene.eq(g)].iloc[0].to_dict() if g in set(de.gene) else {"gene":g,"log2FC_GCA_vs_control":float(logcpm.loc[g,GCA].mean()-logcpm.loc[g,CONTROL].mean()),"p_welch":None,"q_bh":None}
  rec.update({"raw_counts_control":counts.loc[g,CONTROL].astype(int).tolist(),"raw_counts_GCA":counts.loc[g,GCA].astype(int).tolist()}); gene_summary.append(rec)
summary={"accession":"GSE241076","design":"MC38 colorectal cancer cells, GCA versus control, n=3 per group","normalization":"log2(CPM+0.5)","gene_filter":"count >=10 in at least 3 samples","tested_genes":int(keep.sum()),"bh_q_lt_0_10":int((de.q_bh<.10).sum()),"locked_program_genes_detected":list(Z.columns),"locked_program_difference_GCA_minus_control":obs,"locked_program_exact_two_sided_permutation_p":perm_p,"exact_permutations":len(diffs),"mechanism_gene_summary":gene_summary,"interpretation_boundary":"This external perturbation dataset tests a GCA-responsive host program. It does not establish that Peptostreptococcus stomatis produces GCA or that the three layers occur in the same participants."}
(OUT/"host_program_summary.json").write_text(json.dumps(summary,indent=2),encoding="utf-8")
print(json.dumps(summary,indent=2))
