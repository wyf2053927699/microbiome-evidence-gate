from pathlib import Path
import os
import gzip,json
import numpy as np,pandas as pd

ROOT=Path(os.environ.get("MICROBIOME_GATE_ROOT", Path(__file__).resolve().parents[1]))
RAW=ROOT/"01_raw_data"/"crc_axis_external"/"GSE241076"; OUT=ROOT/"06_results"/"crc_axis_external_yachida2019"/"host_program_GSE241076"
samples=["con_1","con_2","con_3","fa_1","fa_2","fa_3"]
def read(name):
 with gzip.open(RAW/name,"rt") as f: return pd.read_csv(f,sep="\t").set_index("gene_id")[samples].apply(pd.to_numeric,errors="coerce").fillna(0)
C=read("GSE241076_counts.txt.gz"); F=read("GSE241076_FPKM.txt.gz")
lib=C.sum(); cpm=C.div(lib,axis=1)*1e6; L=np.log2(cpm+.5); keep=(C>=10).sum(axis=1)>=3
corr=L.loc[keep].corr(); corr.to_csv(OUT/"sample_logCPM_correlation.csv")
X=L.loc[keep].T; X=X-X.mean(axis=0); u,s,v=np.linalg.svd(X.to_numpy(),full_matrices=False); pcs=u[:,:2]*s[:2]
pd.DataFrame({"sample":samples,"condition":["control"]*3+["GCA"]*3,"PC1":pcs[:,0],"PC2":pcs[:,1]}).to_csv(OUT/"sample_pca_scores.csv",index=False)
within=[]; between=[]
for i,a in enumerate(samples):
 for j,b in enumerate(samples[:i]):
  (within if (a.startswith("con")==b.startswith("con")) else between).append(float(corr.loc[a,b]))
mechanism=[g for g in ["NR1H4","SOX14","ZDHHC9","CD274"] if g in C.index]
cross=[]
for g in mechanism:
 lc=np.log2(cpm.loc[g]+.5); lf=np.log2(F.loc[g]+.5)
 cross.append({"gene":g,"counts_logCPM_log2FC":float(lc[["fa_1","fa_2","fa_3"]].mean()-lc[["con_1","con_2","con_3"]].mean()),"FPKM_log2FC":float(lf[["fa_1","fa_2","fa_3"]].mean()-lf[["con_1","con_2","con_3"]].mean())})
summary={"library_sizes":{k:int(v) for k,v in lib.items()},"filtered_genes":int(keep.sum()),"median_within_group_correlation":float(np.median(within)),"median_between_group_correlation":float(np.median(between)),"pc1_variance_fraction":float(s[0]**2/(s**2).sum()),"pc2_variance_fraction":float(s[1]**2/(s**2).sum()),"counts_FPKM_mechanism_concordance":cross,"quality_boundary":"Strong condition separation is compatible with a large perturbation response but cannot distinguish treatment from an unrecorded batch because condition and batch are not independently represented."}
(OUT/"quality_audit_summary.json").write_text(json.dumps(summary,indent=2),encoding="utf-8"); print(json.dumps(summary,indent=2))
