"""
04_trajectory_analysis.py
Individual-level metabolite trajectory visualization (TP1-6)
intensity = within-period log2FC: TP1=0 (early baseline), TP4=0 (late baseline)
"""

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import warnings; warnings.filterwarnings('ignore')

BASE = '/sessions/charming-youthful-pasteur/mnt/TMAO_pathway_analysis'
LFC_CSV = f'{BASE}/data/log2FC_values_production.csv'
OUT_DIR = f'{BASE}/output'

TARGETS = {
    2192:  ('Hexanoyl-carnitine (C6)',       'early_only', '#2E7D32'),
    3646:  ('Palmitoyl-carnitine (C16)',      'early_only', '#2E7D32'),
    2732:  ('Decanoyl-carnitine (C10)',       'both',       '#B45309'),
    2939:  ('Lauroyl-carnitine (C12)',        'both',       '#B45309'),
    3325:  ('Tetradecanoyl-carnitine (C14)',  'both',       '#B45309'),
    3926:  ('Oleoyl-carnitine (C18:1)',       'both',       '#B45309'),
    3912:  ('Linoleyl carnitine (C18:2)',     'both',       '#B45309'),
    1294:  ('L-Carnitine (free)',             'late_only',  '#6A1E6A'),
    3788:  ('Glycochenodeoxycholate',         'both',       '#B45309'),
    4629:  ('Val-conjugated cholate',         'both',       '#B45309'),
    4492:  ('Asp-conjugated CDCA',           'both',       '#B45309'),
    17189: ('Gln-conjugated CDCA',           'both',       '#B45309'),
    5279:  ('His-conjugated CDCA',           'early_only', '#2E7D32'),
    3593:  ('Hyodeoxycholic acid',            'late_only',  '#6A1E6A'),
    5209:  ('Met-conjugated CDCA',            'late_only',  '#6A1E6A'),
}
GROUP_COLOR = {'early_only': '#2E7D32', 'late_only': '#6A1E6A', 'both': '#B45309'}
GROUP_LABEL = {'early_only': 'Early_only', 'late_only': 'Late_only', 'both': 'Both'}

print("Loading trajectory data...")
lfc = pd.read_csv(LFC_CSV)
traj = lfc[lfc['Alignment_ID'].isin(list(TARGETS.keys()))].copy()
print(f"  {len(traj)} rows, {traj['Alignment_ID'].nunique()} metabolites")

def get_pivot(aid):
    return traj[traj['Alignment_ID']==aid].pivot(
        index='Subject', columns='Timepoint', values='intensity')

def spaghetti(ax, aid, show_xlabel=True, show_ylabel=True):
    label, group, color = TARGETS[aid]
    df = get_pivot(aid)
    if df.empty:
        ax.text(0.5,0.5,'No data',ha='center',va='center',transform=ax.transAxes); return
    tps = sorted(df.columns); x = np.array(tps)

    for subj in df.index:
        y = df.loc[subj, tps].values
        ax.plot(x[:3], y[:3], color='#AAAAAA', lw=0.55, alpha=0.35, zorder=1)
        ax.plot(x[3:], y[3:], color='#AAAAAA', lw=0.55, alpha=0.35, zorder=1)

    mean = df[tps].mean(); sd = df[tps].std()
    for px in [x[:3], x[3:]]:
        ix = list(px); m = mean[ix].values; s = sd[ix].values
        ax.fill_between(px, m-s, m+s, color=color, alpha=0.18, zorder=2)
        ax.plot(px, m, color=color, lw=2.2, marker='o', ms=4.5, zorder=3)

    ax.axvline(3.5, color='#455A64', lw=1.0, ls='--', alpha=0.6, zorder=4)
    ax.axhline(0, color='#BBBBBB', lw=0.8, zorder=0)
    ax.axvspan(0.5,3.5, alpha=0.04, color='#1565C0', zorder=0)
    ax.axvspan(3.5,6.5, alpha=0.04, color='#B71C1C', zorder=0)
    ax.set_xlim(0.5,6.5); ax.set_xticks([1,2,3,4,5,6])
    if show_xlabel:
        ax.set_xticklabels(['1\n(E)','2','3','4\n(L)','5','6'], fontsize=7.5)
    else:
        ax.set_xticklabels([])
    if show_ylabel:
        ax.set_ylabel('log2FC\n(within-period)', fontsize=8)
    ax.text(0.02,0.97, GROUP_LABEL[group], transform=ax.transAxes,
            fontsize=7, fontweight='bold', color='white', va='top', ha='left',
            bbox=dict(boxstyle='round,pad=0.2',facecolor=color,alpha=0.85,edgecolor='none'))
    ax.set_title(label, fontsize=8.5, fontweight='bold', pad=3)
    ax.tick_params(axis='both', labelsize=7.5)
    ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)

patch_e = mpatches.Patch(color=GROUP_COLOR['early_only'], label='Early_only')
patch_b = mpatches.Patch(color=GROUP_COLOR['both'],       label='Both')
patch_l = mpatches.Patch(color=GROUP_COLOR['late_only'],  label='Late_only')

# ── Fig 1: Carnitine series ────────────────────────────────
print("Fig 1: Carnitine series...")
carn_ids = [2192, 2732, 2939, 3325, 3646, 3926, 3912, 1294]
fig1, axes1 = plt.subplots(2,4,figsize=(14,6.5))
fig1.suptitle('Carnitine Metabolite Trajectories  ·  Individual subjects (grey) · Mean ± SD',
              fontsize=11, fontweight='bold', y=1.01)
for i, aid in enumerate(carn_ids):
    r, c = divmod(i, 4)
    spaghetti(axes1[r][c], aid, show_xlabel=(r==1), show_ylabel=(c==0))
fig1.legend(handles=[patch_e,patch_b,patch_l], loc='lower right',
            fontsize=8.5, framealpha=0.85, bbox_to_anchor=(1.0,0.0))
plt.tight_layout(rect=[0,0,1,0.98])
fig1.savefig(f'{OUT_DIR}/fig1_carnitine_trajectories.png', dpi=150, bbox_inches='tight', facecolor='white')
plt.close(fig1)
print(f"  -> {OUT_DIR}/fig1_carnitine_trajectories.png")

# ── Fig 2: Bile acids ──────────────────────────────────────
print("Fig 2: Bile acid series...")
bile_ids = [3788, 4629, 4492, 17189, 5279, 3593, 5209]
fig2, axes2 = plt.subplots(2,4,figsize=(14,6.5))
fig2.suptitle('Bile Acid Metabolite Trajectories  ·  Individual subjects (grey) · Mean ± SD',
              fontsize=11, fontweight='bold', y=1.01)
for i, aid in enumerate(bile_ids):
    r, c = divmod(i, 4)
    spaghetti(axes2[r][c], aid, show_xlabel=(r==1 or i>=4), show_ylabel=(c==0))
axes2.flat[7].set_visible(False)
fig2.legend(handles=[patch_e,patch_b,patch_l], loc='lower right',
            fontsize=8.5, framealpha=0.85, bbox_to_anchor=(1.0,0.0))
plt.tight_layout(rect=[0,0,1,0.98])
fig2.savefig(f'{OUT_DIR}/fig2_bileacid_trajectories.png', dpi=150, bbox_inches='tight', facecolor='white')
plt.close(fig2)
print(f"  -> {OUT_DIR}/fig2_bileacid_trajectories.png")

# ── Fig 3: Representative 4-panel ─────────────────────────
print("Fig 3: Representative comparison...")
rep = [(3646,'Palmitoyl-carnitine (C16)\nearly_only'),
       (1294,'L-Carnitine (free)\nlate_only'),
       (3912,'Linoleyl carnitine (C18:2)\nboth'),
       (3788,'Glycochenodeoxycholate\nboth')]
fig3, axes3 = plt.subplots(1,4,figsize=(14,4.5))
fig3.suptitle('Representative Trajectories by Classification Group', fontsize=11, fontweight='bold')
for ax, (aid, title) in zip(axes3, rep):
    _, group, color = TARGETS[aid]
    df = get_pivot(aid); tps = sorted(df.columns); x = np.array(tps)
    for subj in df.index:
        y = df.loc[subj, tps].values
        ax.plot(x[:3], y[:3], color='#CCCCCC', lw=0.7, alpha=0.5, zorder=1)
        ax.plot(x[3:], y[3:], color='#CCCCCC', lw=0.7, alpha=0.5, zorder=1)
    mean = df[tps].mean(); sd = df[tps].std()
    for px in [x[:3], x[3:]]:
        ix = list(px); m = mean[ix].values; s = sd[ix].values
        ax.fill_between(px, m-s, m+s, color=color, alpha=0.2, zorder=2)
        ax.plot(px, m, color=color, lw=2.5, marker='o', ms=6, zorder=3)
    ax.axvline(3.5, color='#455A64', lw=1.0, ls='--', alpha=0.7)
    ax.axhline(0, color='#AAAAAA', lw=0.8)
    ax.axvspan(0.5,3.5, alpha=0.04, color='#1565C0')
    ax.axvspan(3.5,6.5, alpha=0.04, color='#B71C1C')
    ax.set_xlim(0.5,6.5); ax.set_xticks([1,2,3,4,5,6])
    ax.set_xticklabels(['1\n(Early)','2','3','4\n(Late)','5','6'], fontsize=8)
    ax.set_title(title, fontsize=9, fontweight='bold', pad=4)
    ax.set_ylabel('log2FC (within-period)', fontsize=8)
    ax.tick_params(labelsize=8)
    ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)
    ax.text(0.97,0.03,'n=40',transform=ax.transAxes,fontsize=7.5,color='#888888',ha='right',va='bottom')
plt.tight_layout()
fig3.savefig(f'{OUT_DIR}/fig3_representative_comparison.png', dpi=150, bbox_inches='tight', facecolor='white')
plt.close(fig3)
print(f"  -> {OUT_DIR}/fig3_representative_comparison.png")

# ── Fig 4: Mean response bar summary ──────────────────────
print("Fig 4: Mean response summary...")
rows = []
for aid, (label, group, color) in TARGETS.items():
    df = get_pivot(aid)
    if df.empty: continue
    em = df[[2,3]].mean(axis=1)
    lm = df[[5,6]].mean(axis=1)
    rows.append({'label':label,'group':group,'color':color,
                 'early_mean':em.mean(),'early_se':em.std()/np.sqrt(40),
                 'late_mean':lm.mean(), 'late_se': lm.std()/np.sqrt(40)})
summ = pd.DataFrame(rows).sort_values(['group','early_mean'],ascending=[True,False])

fig4, axes4 = plt.subplots(1,2,figsize=(14,6),sharey=False)
fig4.suptitle('Mean Within-Period Response (TP2+3 vs TP1 for early; TP5+6 vs TP4 for late)  ±SEM, n=40',
              fontsize=10.5, fontweight='bold')
for ax, col, se_col, title in [
    (axes4[0],'early_mean','early_se','Early period (X)'),
    (axes4[1],'late_mean', 'late_se', 'Late period (Y)'),
]:
    ax.barh(summ['label'], summ[col], xerr=summ[se_col],
            color=[GROUP_COLOR[g] for g in summ['group']],
            alpha=0.82, edgecolor='white', height=0.65,
            error_kw={'elinewidth':1.2,'capsize':3,'ecolor':'#555555'})
    ax.axvline(0, color='#455A64', lw=0.8)
    ax.set_xlabel('Mean log2FC', fontsize=9)
    ax.set_title(title, fontsize=10, fontweight='bold')
    ax.tick_params(axis='y', labelsize=8)
    ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)
fig4.legend(handles=[patch_e,patch_b,patch_l], loc='lower right', fontsize=9, framealpha=0.85)
plt.tight_layout()
fig4.savefig(f'{OUT_DIR}/fig4_mean_response_summary.png', dpi=150, bbox_inches='tight', facecolor='white')
plt.close(fig4)
print(f"  -> {OUT_DIR}/fig4_mean_response_summary.png")

print("\nDone.")
