from pathlib import Path
import json,numpy as np
from scipy.stats import spearmanr
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import PercentFormatter
P=Path(__file__).parent; d=json.loads((P/'plot_data.json').read_text())
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.labelsize':10,'axes.titlesize':11,'axes.titleweight':'bold','axes.spines.top':False,'axes.spines.right':False,'axes.edgecolor':'#444444','xtick.labelsize':9,'ytick.labelsize':9,'savefig.facecolor':'white','figure.facecolor':'white','axes.facecolor':'white','pdf.fonttype':42,'svg.fonttype':'none'})
N='#286777'; I='#BA6438'; G='#788997'; DARK='#1F3038'
def data(key):return np.array([[float(y) for y in x.split(',')] for x in d[key].split(';')])
s=data('scatter'); o=data('observations'); assert s.shape==(239,2) and s[:,1].sum()==5741 and s[:,0].sum()==26530 and abs(spearmanr(*s.T).statistic-.5788964068803357)<1e-12
assert o.shape==(84,3) and np.all(o[:,1:].sum(1)>0) and tuple(o[:,1:].sum(0))==(916,193)
assert [int(s[j:k,1].sum()) for j,k in [(0,60),(60,120),(120,180),(180,239)]]==[644,1077,1749,2271]
def save(f,name):
 for ext in ['png','pdf']:f.savefig(P/f'{name}.{ext}',dpi=450,bbox_inches='tight',pad_inches=.12)
 plt.close(f)
def title(a,t):a.set_title(t,loc='left',pad=10)
def grid(a,axis='y'):a.grid(axis=axis,color='#e1e5e7',lw=.65);a.set_axisbelow(True)
def err(a,est,y,c=N,marker='o'):
 x,l,h=est
 if l is None:a.plot(x,y,marker,color=c,ms=5);return
 a.errorbar(x,y,xerr=[[x-l],[h-x]],fmt=marker,mfc='white' if marker=='s' else c,mec=c,color=c,ms=5,capsize=3,lw=1.25)
def forest(a,rows,labs,cols=None):
 for y,r in enumerate(rows):err(a,r,len(rows)-1-y,cols[y] if cols else N)
 a.set_yticks(range(len(rows)-1,-1,-1),labs);a.set_ylim(-.65,len(rows)-.35);a.axvline(0,color=G,ls='--',lw=.9);grid(a,'x')
# Fig 1: exact counts, no model intervals on descriptive pooled totals
f,ax=plt.subplots(1,3,figsize=(7.6,3.0),layout='constrained');b=np.array(d['baseline']);T,n,ad=b.T
for a,t,v,lim in zip(ax,['(a) Native coverage','(b) Gaps filled','(c) Relative gain'],[100*n/T,100*ad/(T-n),100*ad/n],[50,50,140]):
 title(a,t);a.bar([0,1],v,color=[N,G],width=.62);a.set_xticks([0,1],['Island','Mainland']);a.set_ylim(0,lim);a.set_ylabel('Percent (%)');grid(a)
 for x,z in enumerate(v):a.text(x,z+.025*lim,f'{z:.1f}%',ha='center',fontsize=9)
save(f,'fig1')
# Fig2: only main results; morphology reported as numerical result in S3.
f,ax=plt.subplots(1,2,figsize=(7.6,2.9),layout='constrained');a,b=ax;title(a,'(a) Ports and plant richness');forest(a,d['port'],['Introduced','Native','Difference'],[I,N,G]);a.set_xlim(-.16,.58);a.set_xlabel('Coefficient per log₂(1 + ports)');title(b,'(b) Floristic similarity');b.barh([1,0],d['jaccard'],color=[I,N],height=.5);b.set_yticks([1,0],['Introduced','Native']);b.set_xlim(0,.060);b.set_xlabel('Median Jaccard similarity');grid(b,'x')
for y,v in zip([1,0],d['jaccard']):b.text(v+.001,y,f'{v:.4f}',va='center',fontsize=9)
save(f,'fig2')
# Fig3 exact point predictions only. No interpolation, no fabricated CIs.
f,ax=plt.subplots(2,2,figsize=(7.6,5.9),layout='constrained');a,b,c,e=ax.ravel();xx=[100,1000,4000]
for aa in [a,b]:aa.set_xscale('log');aa.set_xticks(xx,['100','1,000','4,000']);aa.minorticks_off();aa.set_xlim(65,6500);aa.set_xlabel('Mainland distance (km)');grid(aa)
title(a,'(a) Mapped host coverage')
for vals,col,l in zip(d['count_predictions'],[N,I],['Native only','Native + introduced']):a.scatter(xx,vals,c=col,s=29,label=l,zorder=3)
a.set_ylim(0,82);a.set_ylabel('Fitted number of butterflies');a.legend(frameon=False,loc='upper left',fontsize=8)
title(b,'(b) Introduced-only share');b.scatter(xx,100*np.array(d['shares']),color=I,s=29,zorder=3);b.set_ylim(0,100);b.set_ylabel('Fitted share (%)')
for x,y in zip(xx,100*np.array(d['shares'])):b.annotate(f'{y:.1f}%',(x,y),xytext=(0,6),textcoords='offset points',ha='center',fontsize=9)
title(c,'(c) Native host taxonomic rank')
for y,r in enumerate(d['hierarchy']):err(c,r[:3],2-y+.10,N);err(c,r[3:],2-y-.10,N,'s')
c.set_yticks([2,1,0],['Family','Genus','Species']);c.set_ylim(-.55,2.55);c.set_xlim(-.33,.16);c.axvline(0,color=G,ls='--',lw=.8);grid(c,'x');c.set_xlabel('Standardized distance coefficient');c.plot([],[],'o',color=N,label='Base');c.plot([],[],'s',mfc='white',mec=N,label='MEM10');c.legend(frameon=False,fontsize=8,loc='lower right')
title(e,'(d) Contemporary − native slope');err(e,d['paired_shift'],0,G);e.axvline(0,color=G,ls='--',lw=.8);e.set_xlim(-.11,.48);e.set_ylim(-.7,.7);e.set_yticks([]);grid(e,'x');e.set_xlabel('Paired standardized slope difference');e.text(.186755,.25,'+0.187',ha='center');e.text(.186755,-.27,'95% CI: −0.017 to +0.391\nClustered p = 0.068',ha='center',fontsize=9)
save(f,'fig3')
# Fig4: all 239 source coordinate pairs; only exact cumulative anchors.
f=plt.figure(figsize=(7.6,5.7),layout='constrained');gs=f.add_gridspec(2,2,height_ratios=[1.35,.9]);a=f.add_subplot(gs[0,0]);b=f.add_subplot(gs[0,1]);c=f.add_subplot(gs[1,:]);title(a,'(a) Host-plant contributions');z=np.array(d['top_k']);a.plot(z[:,0],100*z[:,1],'o:',color=I,lw=1.4,ms=4);a.set_xscale('log');a.set_ylim(0,105);a.set_xticks([1,10,38,100,670],['1','10','38','100','670']);a.minorticks_off();a.set_xlabel('Top k contributing plants');a.set_ylabel('Cumulative share (%)');a.axhline(50,color=G,ls='--',lw=.8);a.axvline(38,color=G,ls='--',lw=.8);a.text(.04,.92,'38 plants reach 50%',transform=a.transAxes,fontsize=8,va='top');grid(a)
title(b,'(b) Native extent and island gains');b.scatter(*s[60:180].T,s=15,c=G,alpha=.5,linewidths=0);b.scatter(*s[:60].T,s=17,c=N,alpha=.8,linewidths=0,label='Narrowest 60');b.scatter(*s[180:].T,s=17,c=I,alpha=.8,linewidths=0,label='Broadest 59');b.set_xlabel('Global regions with native hosts');b.set_ylabel('New island host regions');b.set_xlim(-5,320);b.set_ylim(-2,100);b.text(.03,.97,'Spearman ρ = 0.579',transform=b.transAxes,va='top',fontsize=9);b.legend(frameon=False,fontsize=7.5,loc='upper right',bbox_to_anchor=(1,.88));grid(b)
title(c,'(c) Island gaps filled by native-resource extent');q=np.array(d['gap_groups']);v=100*q[:,1]/q[:,2];c.bar(range(4),v,color=[N,G,G,I],width=.6);c.set_xticks(range(4),['Q1 (60)','Q2 (60)','Q3 (60)','Q4 (59)']);c.set_ylim(0,55);c.set_ylabel('Gaps filled (%)');c.set_xlabel('Native-resource extent quartile (number of butterflies)');grid(c)
for x,y in enumerate(v):c.text(x,y+1.2,f'{y:.1f}%',ha='center',fontsize=9)
save(f,'fig4')
# Fig5: source-defined logit fits, exact observed share dots; log1p x-axis.
f,a=plt.subplots(figsize=(7.6,4.25),layout='constrained');dist,n,i=o.T;a.scatter(np.log1p(dist),100*i/(n+i),s=9+1.9*(n+i),c=I,alpha=.57,edgecolors='white',linewidths=.35,label='Recorded regional shares');fit=d['fit'];x=np.linspace(np.log1p(fit['min_km']),np.log1p(fit['max_km']),200)
for k,col,ls,lab in [('potential',N,'--','Mapped potential fit'),('observed',I,'-','Recorded fit')]:
 mu,b=fit[k];y=100/(1+np.exp(-(mu+b*(x-fit['mean_logdistance'])/fit['sd_logdistance'])));a.plot(x,y,color=col,ls=ls,lw=1.8,label=lab)
a.set_xticks(np.log1p([0,1,10,100,1000,6000]),['0','1','10','100','1,000','6,000']);a.set_xlim(-.1,np.log1p(6500));a.set_ylim(-3,104);a.set_xlabel('Mainland distance (km; log1p scale)');a.set_ylabel('Introduced-only share (%)');a.legend(frameon=False,loc='upper left',fontsize=9);grid(a);save(f,'fig5')
# SI: same-unit port origin contrasts, all intervals include zero
f,a=plt.subplots(figsize=(7.6,3.2),layout='constrained');forest(a,d['port_ratio'],['126: area + division','126: + distance, latitude','98: area, division, effort','98: + geography, no effort','98: + geography, effort'],[N,N,I,I,I]);a.set_xlim(-.14,.64);a.set_xlabel('Port coefficient: log(1 + introduced) − log(1 + native)');save(f,'figS33')
# SI marine circuit, same source and approximate original one-island intervals.
f,a=plt.subplots(figsize=(7.6,2.9),layout='constrained');forest(a,d['marine'],['One unseen: introduced','One unseen: native','Both unseen: introduced','Both unseen: native'],[I,N,I,N]);a.set_xlim(-2.25,6.7);a.set_xlabel('Relative reduction in held-out MSE (%)');save(f,'figS32')
# SI transfer: display all adverse results and distinguish bootstrap vs MC quantiles.
f=plt.figure(figsize=(7.6,6.0),layout='constrained');gs=f.add_gridspec(3,1,height_ratios=[1,1,1.25]);a,b,c=[f.add_subplot(gs[i]) for i in range(3)]
for aa,tt,rr,ll in [(a,'(a) Regional holdouts',d['community'],['Unthinned (bootstrap)','Dataset-days (MC)','Record counts (MC)']),(b,'(b) Unseen botanical divisions',d['transfer'],['Unthinned (no interval)','Dataset-days (MC)','Record counts (MC)']),(c,'(c) South America holdout',d['south_america'],['Within: dataset-days','Cross-world: dataset-days','Within: record counts','Cross-world: record counts'])]:
 title(aa,tt);forest(aa,rr,ll,[G,I,I] if len(rr)==3 else [N,I,N,I]);aa.set_xlabel('Relative reduction in held-out MSE (%)')
a.set_xlim(-.18,.85);b.set_xlim(-.18,1.8);c.set_xlim(-8.6,2);save(f,'figS51')
print('PASS: source coordinate counts, sums, Spearman correlation, quartile totals and GBIF category counts')
