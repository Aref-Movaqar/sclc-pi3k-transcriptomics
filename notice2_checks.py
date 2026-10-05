import pandas as pd, numpy as np
from scipy import stats
rng=np.random.default_rng(20261005)
c1=pd.read_csv('results/Table_S1_DE_PIK75_vs_DMSO.csv').set_index('gene'); c2=pd.read_csv('results/Table_S1_DE_TGX221_vs_DMSO.csv').set_index('gene')
hl=pd.read_csv('ext/human_halflife_PC1.csv').drop_duplicates('gene').set_index('gene')['halflife_PC1']
print('=== POINT 3: stratified half-life test')
d=c1.join(hl,how='inner'); d['pct']=d.halflife_PC1.rank(pct=True)*100
for lab,m in [('shortest 2%',d.pct<=2),('shortest 5%',d.pct<=5),('shortest 10%',d.pct<=10),('10-50%',(d.pct>10)&(d.pct<=50)),('longest 50%',d.pct>50),('longest 10%',d.pct>90)]:
    s=d[m]; print(f'{lab:13s} n={len(s):5d} mean log2FC={s.logFC.mean():+.3f} median={s.logFC.median():+.3f}  %down>=2fold={100*((s.logFC<=-1)&(s["adj.P.Val"]<.05)).mean():.1f}  %up>=2fold={100*((s.logFC>=1)&(s["adj.P.Val"]<.05)).mean():.2f}')
u=stats.mannwhitneyu(d[d.pct<=5].logFC,d[d.pct>5].logFC); print('shortest 5pct vs rest: Mann-Whitney p=%.2g'%u.pvalue)
dn=d[(d.logFC<=-1)&(d['adj.P.Val']<.05)]; print('Two-fold decreased genes with half-life data: n=%d, median half-life percentile=%.0f (50 expected if unrelated); share in shortest 10%%: %.1f%%'%(len(dn),dn.pct.median(),100*(dn.pct<=10).mean()))
print('Half-life percentile of named genes:'); print(d.loc[[g for g in ['SESN3','GADD45A','CITED2','DDIT4','TXNIP','CCNE2','HEXIM1','MYC','FOS','JUN','IER2','MCL1','BRCA1','HMGCS1','SQLE','HMGCR','LDLR','SCD','TRIB3'] if g in d.index],['logFC','pct']].round(2).T.to_string())
d2=c2.join(hl,how='inner'); d2['pct']=d2.halflife_PC1.rank(pct=True)*100
print('TGX-221: shortest 5%% mean %.3f, rest %.3f'%(d2[d2.pct<=5].logFC.mean(), d2[d2.pct>5].logFC.mean()))

print('\n=== POINT 4: size-matched random-gene-set null for mean log2FC')
eff=pd.read_csv('results/Table_S16_set_effect_sizes.csv'); lf=c1.logFC.values; lf2=c2.logFC.values
rows=[]
for _,r in eff.iterrows():
    n=int(r.n); nul=np.array([lf[rng.choice(len(lf),n,replace=False)].mean() for _ in range(10000)]); nul2=np.array([lf2[rng.choice(len(lf2),n,replace=False)].mean() for _ in range(10000)])
    p=(1+(np.abs(nul-nul.mean())>=abs(r.PIK_mean_log2FC-nul.mean())).sum())/10001; p2=(1+(np.abs(nul2-nul2.mean())>=abs(r.TGX_mean_log2FC-nul2.mean())).sum())/10001
    rows.append((r.set[:38],n,r.PIK_mean_log2FC,round((r.PIK_mean_log2FC-nul.mean())/nul.std(),1),p,r.PIK_roast,r.TGX_mean_log2FC,round((r.TGX_mean_log2FC-nul2.mean())/nul2.std(),1),p2))
t=pd.DataFrame(rows,columns=['set','n','PIK_mean','PIK_z_vs_random','PIK_emp_p','PIK_roast','TGX_mean','TGX_z','TGX_emp_p']); pd.set_option('display.width',250); print(t.to_string(index=False))

print('\n=== POINT 7: are the Table 2 SREBP genes cherry-picked? anchored sets')
import subprocess,os
gm={}
for lib in ['Reactome_2022','MSigDB_Hallmark_2020']:
    fn=f'ext/{lib}.gmt'
    if not os.path.exists(fn): subprocess.run(['curl','-sS','-m','120',f'https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName={lib}','-o',fn])
    for l in open(fn):
        f=l.rstrip('\n').split('\t'); gm[f[0]]=[x for x in f[2:] if x]
t2=['SQLE','HMGCS1','HMGCR','LDLR','SCD','LSS','FDFT1','ACACA','INSIG1','FASN']
for k in gm:
    if ('SREBP' in k or 'SREBF' in k or 'Cholesterol Biosynthesis R-HSA' in k or k=='Cholesterol Homeostasis'):
        g=[x for x in gm[k] if x in c1.index]; s=c1.loc[g]
        print(f'{k[:75]:75s} n={len(g):3d} mean={s.logFC.mean():+.2f} median={s.logFC.median():+.2f} down>=2fold&FDR: {((s.logFC<=-1)&(s["adj.P.Val"]<.05)).sum()} ; any-decrease FDR<.05: {((s.logFC<0)&(s["adj.P.Val"]<.05)).sum()} ; up FDR<.05: {((s.logFC>0)&(s["adj.P.Val"]<.05)).sum()} ; Table2 genes in set: {len(set(t2)&set(g))}/10 ; TGX mean={c2.loc[g].logFC.mean():+.2f}')
k=[x for x in gm if 'Activation Of Gene Expression By SREBF' in x][0]; g=[x for x in gm[k] if x in c1.index]
print('All genes of "%s":'%k); print(pd.DataFrame({'PIK':c1.loc[g].logFC.round(2),'FDR':c1.loc[g]['adj.P.Val'].map(lambda v:'%.1g'%v),'TGX':c2.loc[g].logFC.round(2)}).sort_values('PIK').T.to_string())
rank=c1.logFC.rank(pct=True)*100; print('Percentile rank of Table 2 SREBP genes among all 17,458 log2FC (0 = most decreased):', {g:round(rank[g],1) for g in t2})

print('\n=== POINT 3b: LINCS L1000 consensus signatures (Enrichr library)')
up={};down={}
for l in open('ext/lincs_consensus.gmt'):
    f=l.rstrip('\n').split('\t'); name,dirn=f[0].rsplit(' ',1); (up if dirn=='Up' else down)[name]=[x for x in f[2:] if x]
def conn(c):
    out={}
    for n in up:
        if n not in down: continue
        a=[x for x in up[n] if x in c.index]; b=[x for x in down[n] if x in c.index]
        if len(a)<20 or len(b)<20: continue
        out[n]=c.loc[a,'t'].mean()-c.loc[b,'t'].mean()
    s=pd.Series(out); return (s-s.mean())/s.std()
s1=conn(c1); s2=conn(c2); pr1=s1.rank(pct=True)*100; pr2=s2.rank(pct=True)*100
print('compounds scored:',len(s1))
cls={'PIK-75 itself':['PIK-75'],'PI3Kalpha/pan-PI3K':['Alpelisib','A-66','Taselisib','GDC-0941','Buparlisib','ZSTK-474','LY-294002','Wortmannin'],'PI3K/mTOR dual & mTOR':['GDC-0980','GSK-2126458','PI-103','AZD-8055','Torin-1','Torin-2','Sirolimus','Everolimus','Temsirolimus'],'AKT':['MK-2206'],'PI3Kbeta':['TGX-221','TGX-115','AZD-6482'],'CDK9/transcription':['Alvocidib','Dinaciclib','AT-7519','Actinomycin-D','Dactinomycin','Triptolide'],'DNA-PK':['NU-7441','KU-0060648']}
for k,v in cls.items():
    v=[x for x in v if x in s1.index]
    print(f'{k:24s} median z (PIK-75 contrast)={s1[v].median():+.2f} median percentile={pr1[v].median():.0f} | TGX contrast z={s2[v].median():+.2f} pct={pr2[v].median():.0f} :: '+', '.join(f'{x} {s1[x]:+.1f}' for x in v))
print('Top 25 most similar compounds to H69 PIK-75 response:'); print(', '.join(f'{n} ({z:+.1f})' for n,z in s1.sort_values(ascending=False).head(25).items()))
print('Top 12 for TGX-221 contrast:'); print(', '.join(f'{n} ({z:+.1f})' for n,z in s2.sort_values(ascending=False).head(12).items()))
pd.DataFrame({'z_PIK75_contrast':s1.round(3),'pct_PIK75':pr1.round(1),'z_TGX221_contrast':s2.round(3),'pct_TGX221':pr2.round(1)}).sort_values('z_PIK75_contrast',ascending=False).to_csv('results/Table_S22_LINCS_connectivity.csv')
t.to_csv('results/Table_S23_random_set_null.csv',index=False)
