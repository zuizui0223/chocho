from pathlib import Path
import pandas as pd,numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
r=Path.cwd();p=pd.read_csv('outputs/first_record_lag_event_summary.csv');p=p.sort_values(['species','region']).reset_index(drop=True)
plt.rcParams['font.family']='Meiryo'
f,ax=plt.subplots(figsize=(12,9))
for i,a in p.iterrows():
 if pd.notna(a.earliest_host_lag_lower):
  lo,hi=a.earliest_host_lag_lower,a.earliest_host_lag_upper
  ax.plot([lo,hi],[i,i],color='#416a8c',lw=3);ax.scatter((lo+hi)/2,i,color='#416a8c',s=25)
  ax.text(hi+5,i,f'{lo:g}' if lo==hi else f'{lo:g}～{hi:g}',va='center',fontsize=9)
 else:ax.text(5,i,'年代不足：順序を判定しない',color='#777',va='center',fontsize=9)
ax.set_yticks(range(len(p)),p.species+' / '+p.region);ax.invert_yaxis();ax.axvline(0,color='#777',linestyle='--');ax.set_xlim(-100,270);ax.set_xlabel('チョウ初記録 − 対応できた食草の最も早い初記録（年）\n正：食草が先に記録された');ax.spines[['top','right']].set_visible(False)
f.suptitle('食草とチョウの初記録の時間差：年代付き部分集合',fontsize=14)
f.text(.03,.025,f'候補225種中{p.species.nunique()}種・{len(p)}種×地域（地理不整合を除外）。線は原資料の年代幅であり、信頼区間ではない。\n導入・到来時期ではなく、同一地域で対応できた初記録の差。未収録の食草・在来食草・栽培履歴は含まれない。\n出典：Alien Species First Records v3.1 (Seebens et al.; DOI:10.5281/zenodo.10039630) と既存食草カタログ。',fontsize=9)
f.tight_layout(rect=[0,.13,1,.95]);f.savefig('outputs/first_record_lag_intervals.png',dpi=160)
