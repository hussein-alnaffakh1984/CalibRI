# Publication figures sized to the A4 text width (6.3 in), 600 dpi, no in-figure titles (captions carry them)
import pandas as pd, numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':8,'axes.spines.top':False,'axes.spines.right':False,
  'axes.edgecolor':'#52514e','xtick.color':'#3a3a38','ytick.color':'#3a3a38','axes.titlesize':8.5,'axes.labelsize':8,
  'xtick.labelsize':7.5,'ytick.labelsize':7.5,'legend.fontsize':7.5,'savefig.dpi':600})
OUT='figures/'; W_IN=6.3
C={'refineR':'#2a78d6','reflimR':'#eb6834','RINet':'#1baf7a'}
R=pd.read_csv('work/ana/disjoint_results.csv'); R=R[R.tol=='tol']
R['method']=R.method.replace({'refineR(UM-scale)':'refineR'})
AN=['TSH','Tg','TT4','TT3','HbA1c']; ANL={'TSH':'TSH','Tg':'Tg','TT4':'TT4','TT3':'TT3','HbA1c':'HbA1c'}
# Fig 1
nat=R[R.scheme=='native'].groupby(['method','analyte','limit','n'])['cov'].mean().reset_index()
fig,axes=plt.subplots(1,3,figsize=(W_IN,3.9),sharey=True,gridspec_kw={'wspace':0.08})
cmap=LinearSegmentedColormap.from_list('b',['#f2f6fc','#9fc3ee','#2a78d6','#0b3b78'])
rows=[f'{ANL[a]} {l}' for a in AN for l in ('LRL','URL')]
for ax,me in zip(axes,['reflimR','RINet','refineR']):
    M=np.array([[nat[(nat.method==me)&(nat.analyte==a)&(nat.limit==l)&(nat.n==n)]['cov'].values[0] for n in (300,500,1000)] for a in AN for l in ('LRL','URL')])
    im=ax.imshow(M,cmap=cmap,vmin=0,vmax=1,aspect='auto')
    for i in range(M.shape[0]):
        for j in range(3):
            ax.text(j,i,f'{M[i,j]:.2f}',ha='center',va='center',fontsize=7,color='#0b0b0b' if M[i,j]<0.55 else 'white')
    ax.set_xticks(range(3)); ax.set_xticklabels(['300','500','1000']); ax.set_xlabel('Sample size (n)')
    ax.set_yticks(range(len(rows))); ax.set_yticklabels(rows); ax.tick_params(length=0)
    ax.set_title(me+('\n(uncertainty margin*)' if me=='refineR' else '\n(95% CI)'))
    for s in ax.spines.values(): s.set_visible(False)
cb=fig.colorbar(im,ax=axes,fraction=0.025,pad=0.02); cb.set_label('Coverage'); cb.outline.set_visible(False)
fig.savefig(OUT+'Figure1.png',bbox_inches='tight'); plt.close()
# Fig 2
st=['native','sim-only','sim+otherHalf(all analytes)','sim+otherHalf(other analytes)']
lab=['Native','Conformal:\nsimulation only','Conformal:\nsim + other-half\nreal (all analytes)','Conformal:\nsim + other-half\nreal (other analytes)']
fig,ax=plt.subplots(figsize=(W_IN,3.6)); rng=np.random.default_rng(1)
for i,s in enumerate(st):
    for k,me in enumerate(['reflimR','RINet','refineR']):
        v=R[(R.scheme==s)&(R.method==me)]['cov'].values; c=i+(k-1)*0.25
        ax.scatter(c+rng.uniform(-0.07,0.07,len(v)),v,s=7,color=C[me],alpha=0.75,linewidths=0,
                   label=(me+(' (UM scale)' if me=='refineR' else ' (95% CI)')) if i==0 else None)
        ax.plot([c-0.1,c+0.1],[np.median(v)]*2,color='#0b0b0b',lw=1.4)
ax.axhline(0.95,color='#52514e',lw=0.9,ls='--')
ax.set_xticks(range(4)); ax.set_xticklabels(lab); ax.set_xlim(-0.5,3.5)
ax.set_ylabel('Empirical coverage per cell'); ax.set_ylim(-0.03,1.03); ax.set_yticks(np.arange(0,1.01,0.2))
ax.grid(axis='y',color='#e8e7e3',lw=0.5); ax.set_axisbelow(True)
h,l=ax.get_legend_handles_labels()
h+= [plt.Line2D([],[],color='#0b0b0b',lw=1.4),plt.Line2D([],[],color='#52514e',lw=0.9,ls='--')]; l+=['median','nominal 0.95']
ax.legend(h,l,frameon=False,ncol=5,loc='lower center',bbox_to_anchor=(0.5,1.0),markerscale=1.8,handlelength=1.5,columnspacing=1.0)
fig.savefig(OUT+'Figure2.png',bbox_inches='tight'); plt.close()
# Fig 3: 2 x 3 grid, legend in the sixth panel
Wd=pd.read_csv('work/ana/widths_d.csv'); Wd=Wd[Wd.limit=='URL']
fig,axes=plt.subplots(2,3,figsize=(W_IN,4.4),sharey=True,gridspec_kw={'hspace':0.55,'wspace':0.12}); axes=axes.ravel()
for ax,a in zip(axes,AN):
    for me in ['reflimR','RINet','refineR']:
        s=Wd[(Wd.analyte==a)&(Wd.method==me)].sort_values('n')
        ax.plot(s.n,s.hw,'-o',color=C[me],lw=1.6,ms=4,label=me+(' (UM scale)' if me=='refineR' else ''))
    ax.set_xscale('log'); ax.set_xticks([300,500,1000]); ax.set_xticklabels(['300','500','1000']); ax.minorticks_off()
    ax.set_yscale('log'); ax.set_yticks([5,10,20,50,100,200]); ax.set_yticklabels(['5','10','20','50','100','200']); ax.minorticks_off()
    ax.set_xlim(250,1200); ax.set_title(ANL[a]); ax.set_xlabel('Sample size (n)'); ax.grid(axis='y',color='#e8e7e3',lw=0.5); ax.set_axisbelow(True)
for ax in (axes[0],axes[3]): ax.set_ylabel('Upward half-width (%)')
h,l=axes[0].get_legend_handles_labels(); axes[5].axis('off'); axes[5].legend(h,l,frameon=False,loc='center',title='Method',title_fontsize=8)
fig.savefig(OUT+'Figure3.png',bbox_inches='tight'); plt.close()
# Fig 4: two stacked panels
S=pd.read_csv('results/tables/samawah_results_schemeD.csv')
grp=[('F','all'),('M','all'),('F','2324'),('F','2425'),('F','2526'),('M','2324'),('M','2425'),('M','2526')]
names={'all':'2023–2026','2324':'2023–24','2425':'2024–25','2526':'2025–26'}
ylab=[f"{'Women' if s=='F' else 'Men'} {names[p]} (n = {int(S[(S.sex==s)&(S.period.astype(str)==p)].n.iloc[0]):,})" for s,p in grp]
fig,axes=plt.subplots(2,1,figsize=(5.2,7.4),gridspec_kw={'hspace':0.35})
for ax,lim,kit,tag in zip(axes,['LRL','URL'],[0.35,5.10],['A','B']):
    for gi,(sx,pr) in enumerate(grp):
        for k,me in enumerate(['reflimR','RINet','refineR']):
            r=S[(S.method==me)&(S.limit==lim)&(S.sex==sx)&(S.period.astype(str)==pr)].iloc[0]; y=gi+(k-1)*0.25
            ax.plot([r.cal_lo,r.cal_hi],[y,y],color=C[me],lw=2.0,solid_capstyle='round')
            ax.plot([r.native_lo,r.native_hi],[y,y],color='#0b0b0b',lw=0.8)
            ax.plot(r.est,y,'o',color=C[me],ms=4,mec='white',mew=0.8)
    ax.axvline(kit,color='#52514e',ls='--',lw=0.9)
    ax.set_yticks(range(len(grp))); ax.set_yticklabels(ylab); ax.set_ylim(len(grp)-0.5,-0.5)
    for b in (1.5,4.5): ax.axhline(b,color='#cfcfcb',lw=0.6)
    ax.set_xscale('log'); tk=[0.1,0.2,0.3,0.5,1.0] if lim=='LRL' else [1.5,2,3,5,7,10,15,20]
    ax.set_xticks(tk); ax.set_xticklabels([f'{t:g}' for t in tk]); ax.minorticks_off()
    ax.set_xlabel('TSH (mIU/L), log scale'); ax.grid(axis='x',color='#e8e7e3',lw=0.5); ax.set_axisbelow(True)
    ax.set_title(f'{tag}. '+('Lower reference limit' if lim=='LRL' else 'Upper reference limit')+f'\n(dashed line: manufacturer limit {kit:.2f} mIU/L)',loc='left')
h=[plt.Line2D([],[],color=C[m],lw=2.0,marker='o',ms=4,mec='white') for m in ['reflimR','RINet','refineR']]+[plt.Line2D([],[],color='#0b0b0b',lw=0.8)]
fig.legend(h,['reflimR (calibrated 95%)','RINet (calibrated 95%)','refineR (calibrated 95%)','native interval'],loc='lower center',ncol=2,frameon=False,bbox_to_anchor=(0.5,-0.015))
fig.subplots_adjust(bottom=0.1)
fig.savefig(OUT+'Figure4.png',bbox_inches='tight'); plt.close()
# check limits fall inside axes
for lim in ['LRL','URL']:
    s=S[S.limit==lim]; print(lim, s[['cal_lo','native_lo','est']].min().min(), s[['cal_hi','native_hi','est']].max().max())
print('done')
