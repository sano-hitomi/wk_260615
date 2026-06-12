import pandas as pd
import numpy as np
from scipy.stats import fisher_exact
from scipy.stats import chi2_contingency
from statsmodels.stats.multitest import multipletests
import warnings
warnings.filterwarnings('ignore')

# ─── データ読み込み ───
cand = pd.read_csv('/sessions/charming-youthful-pasteur/mnt/TMAO_pathway_analysis/data/plsda_mechanism_candidates.csv')
ora  = pd.read_csv('/sessions/charming-youthful-pasteur/mnt/TMAO_pathway_analysis/data/plsda_metaboanalyst_ora.csv')
msea = pd.read_csv('/sessions/charming-youthful-pasteur/mnt/TMAO_pathway_analysis/data/plsda_metaboanalyst_msea.csv')

print("Data loaded:", cand.shape, ora.shape, msea.shape)
print("Classes:", cand['mechanism_class'].value_counts().to_dict())

# ─── 1. グループ定義 ───
early = cand[cand['mechanism_class'] == 'Early_specific'].copy()
late  = cand[cand['mechanism_class'] == 'Late_specific'].copy()
rev   = cand[cand['mechanism_class'].isin(['Reversed_Q2','Reversed_Q4'])].copy()
background = cand.copy()

print(f"\nEarly_specific: n={len(early)}, Late_specific: n={len(late)}, Reversed: n={len(rev)}")

# ─── 2. Ontologyエンリッチメント (Fisher's exact test) ───
def run_ora_ontology(subset, background, label, min_count=2):
    """各Ontologyクラスで Fisher's exact test"""
    results = []
    N = len(background)
    n = len(subset)
    
    all_classes = background['Ontology'].dropna()
    all_classes = all_classes[all_classes != '']
    class_counts = all_classes.value_counts()
    
    for cls, K in class_counts.items():
        if K < min_count:
            continue
        k = subset['Ontology'].value_counts().get(cls, 0)
        # contingency table: [k, n-k; K-k, N-n-(K-k)]
        table = [[k, n - k],
                 [K - k, N - n - (K - k)]]
        try:
            _, p = fisher_exact(table, alternative='greater')
        except:
            p = 1.0
        results.append({
            'Ontology': cls,
            'count_in_set': k,
            'set_size': n,
            'count_in_bg': K,
            'bg_size': N,
            'ratio_in_set': round(k/n, 4),
            'ratio_in_bg': round(K/N, 4),
            'fold_enrichment': round((k/n)/(K/N), 3) if K/N > 0 else np.nan,
            'p_value': p,
            'group': label
        })
    
    if not results:
        return pd.DataFrame()
    df = pd.DataFrame(results)
    # FDR correction
    _, pvals_adj, _, _ = multipletests(df['p_value'], method='fdr_bh')
    df['FDR'] = pvals_adj
    df = df.sort_values('p_value').reset_index(drop=True)
    return df

ora_early   = run_ora_ontology(early, background, 'Early_specific')
ora_late    = run_ora_ontology(late,  background, 'Late_specific')
ora_rev     = run_ora_ontology(rev,   background, 'Reversed')

print("\nEarly significant (p<0.05):", (ora_early['p_value'] < 0.05).sum())
print("Late significant (p<0.05):",  (ora_late['p_value']  < 0.05).sum())
print("Reversed significant (p<0.05):", (ora_rev['p_value'] < 0.05).sum())

# ─── 3. Early vs Late 直接比較 ───
def compare_early_vs_late(early, late, min_count=3):
    """Early_specific と Late_specific の間でOntologyの偏りを検定"""
    all_cls = set(early['Ontology'].dropna().tolist()) | set(late['Ontology'].dropna().tolist())
    n_e = len(early)
    n_l = len(late)
    results = []
    
    for cls in all_cls:
        if cls == '':
            continue
        k_e = (early['Ontology'] == cls).sum()
        k_l = (late['Ontology']  == cls).sum()
        if k_e + k_l < min_count:
            continue
        table = [[k_e, n_e - k_e],
                 [k_l, n_l - k_l]]
        try:
            _, p = fisher_exact(table)
        except:
            p = 1.0
        results.append({
            'Ontology': cls,
            'Early_count': k_e, 'Early_n': n_e,
            'Late_count':  k_l, 'Late_n':  n_l,
            'Early_pct': round(100*k_e/n_e, 2),
            'Late_pct':  round(100*k_l/n_l, 2),
            'Early_enriched': k_e/n_e > k_l/n_l,
            'p_value': p
        })
    
    if not results:
        return pd.DataFrame()
    df = pd.DataFrame(results)
    _, pvals_adj, _, _ = multipletests(df['p_value'], method='fdr_bh')
    df['FDR'] = pvals_adj
    df = df.sort_values('p_value').reset_index(drop=True)
    return df

ev_comparison = compare_early_vs_late(early, late)
sig_ev = ev_comparison[ev_comparison['p_value'] < 0.05]
print(f"\nEarly vs Late 有意差あり Ontologyクラス (p<0.05): {len(sig_ev)}")
print(sig_ev[['Ontology','Early_count','Late_count','Early_pct','Late_pct','p_value','FDR']].to_string())

# ─── 4. KEGG pathway対応マップ ───
# 主要OntologyクラスのKEGG pathway対応（手動マッピング）
ontology_to_kegg = {
    # 脂質代謝
    'Triacylglycerols':                          ['map00561 Glycerolipid metabolism'],
    'Phosphatidylcholines':                       ['map00564 Glycerophospholipid metabolism', 'map05231 Choline metabolism in cancer'],
    'Phosphatidylethanolamines':                  ['map00564 Glycerophospholipid metabolism'],
    'Lysophosphatidylcholines':                   ['map00564 Glycerophospholipid metabolism'],
    '1-acyl-sn-glycero-3-phosphocholines':        ['map00564 Glycerophospholipid metabolism'],
    'Lysophosphatidylethanolamines':              ['map00564 Glycerophospholipid metabolism'],
    'Sphingomyelins':                             ['map00600 Sphingolipid metabolism'],
    'Ceramides':                                  ['map00600 Sphingolipid metabolism', 'map04071 Sphingolipid signaling pathway'],
    'Hexosylceramides':                           ['map00600 Sphingolipid metabolism', 'map00601 Glycosphingolipid biosynthesis'],
    'Glycerophosphocholines':                     ['map00564 Glycerophospholipid metabolism'],
    # 胆汁酸
    'Bile acids and derivatives':                 ['map00120 Primary bile acid biosynthesis', 'map00121 Secondary bile acid biosynthesis', 'map04976 Bile secretion'],
    'Dihydroxy bile acids, alcohols and derivatives': ['map00120 Primary bile acid biosynthesis', 'map04976 Bile secretion'],
    'Trihydroxy bile acids, alcohols and derivatives': ['map00120 Primary bile acid biosynthesis'],
    'Glycinated bile acids and derivatives':      ['map00120 Primary bile acid biosynthesis', 'map00121 Secondary bile acid biosynthesis'],
    'Taurinated bile acids and derivatives':      ['map00120 Primary bile acid biosynthesis', 'map00430 Taurine and hypotaurine metabolism'],
    'Bile alcohols':                              ['map00120 Primary bile acid biosynthesis'],
    # アミノ酸・ペプチド
    'Amino acids and derivatives':                ['map00250 Alanine, aspartate and glutamate metabolism', 'map00260 Glycine, serine and threonine metabolism'],
    'Oligopeptides':                              ['map04974 Protein digestion and absorption'],
    'Polypeptides':                               ['map04974 Protein digestion and absorption'],
    # 脂肪酸
    'Fatty acids and conjugates':                 ['map00071 Fatty acid degradation', 'map01040 Biosynthesis of unsaturated fatty acids'],
    'Medium-chain fatty acids':                   ['map00071 Fatty acid degradation'],
    'Long-chain fatty acids':                     ['map00071 Fatty acid degradation', 'map01212 Fatty acid metabolism'],
    'Very long-chain fatty acids':                ['map01040 Biosynthesis of unsaturated fatty acids'],
    'Hydroxy fatty acids':                        ['map00071 Fatty acid degradation'],
    'Eicosanoids':                                ['map00590 Arachidonic acid metabolism', 'map04726 Serotonergic synapse'],
    'Prostaglandins':                             ['map00590 Arachidonic acid metabolism'],
    'Leukotrienes':                               ['map00590 Arachidonic acid metabolism'],
    # N-アシルアミン
    'N-acyl amines':                              ['map00590 Arachidonic acid metabolism', 'map04723 Retrograde endocannabinoid signaling'],
    'N-acylethanolamines':                        ['map04723 Retrograde endocannabinoid signaling'],
    # ステロイド
    'Steroids and steroid derivatives':           ['map00100 Steroid biosynthesis', 'map00140 Steroid hormone biosynthesis'],
    'Triterpenoids':                              ['map00100 Steroid biosynthesis', 'map00909 Sesquiterpenoid and triterpenoid biosynthesis'],
    'Steroidal saponins':                         ['map00100 Steroid biosynthesis'],
    # その他
    'Carboxylic acids and derivatives':           ['map00020 Citrate cycle (TCA cycle)'],
    'Tetracarboxylic acids and derivatives':      ['map00020 Citrate cycle (TCA cycle)'],
    'Indoles and derivatives':                    ['map00380 Tryptophan metabolism'],
    'Tryptophan derivatives':                     ['map00380 Tryptophan metabolism'],
    'Flavonoids':                                 ['map00941 Flavonoid biosynthesis', 'map00942 Anthocyanin biosynthesis'],
    'Coumarins and derivatives':                  ['map00940 Phenylpropanoid biosynthesis'],
    'Lipoxins':                                   ['map00590 Arachidonic acid metabolism'],
    'Hydroxyeicosatetraenoic acids':              ['map00590 Arachidonic acid metabolism'],
    'Resolvins':                                  ['map00590 Arachidonic acid metabolism'],
}

# ─── 5. 統合結果テーブル ───
all_ora = pd.concat([ora_early, ora_late, ora_rev], ignore_index=True)
all_ora['KEGG_pathways'] = all_ora['Ontology'].map(ontology_to_kegg).apply(
    lambda x: '; '.join(x) if isinstance(x, list) else 'No direct mapping')

# ─── 6. 代謝物リスト ───
# Early_specific (VIP高順, 名前で「low score」でない優先)
early_metabolites = early[['Alignment_ID','best_name','best_inchikey','Ontology',
                            'vip_max','loading_early','loading_late','mechanism_class']].copy()
early_metabolites['well_annotated'] = ~early_metabolites['best_name'].str.startswith('low score', na=False)
early_metabolites = early_metabolites.sort_values(['well_annotated','vip_max'], ascending=[False,False])

late_metabolites = late[['Alignment_ID','best_name','best_inchikey','Ontology',
                          'vip_max','loading_early','loading_late','mechanism_class']].copy()
late_metabolites['well_annotated'] = ~late_metabolites['best_name'].str.startswith('low score', na=False)
late_metabolites = late_metabolites.sort_values(['well_annotated','vip_max'], ascending=[False,False])

rev_metabolites = rev[['Alignment_ID','best_name','best_inchikey','Ontology',
                        'vip_max','loading_early','loading_late','mechanism_class','mechanism_score']].copy()
rev_metabolites['well_annotated'] = ~rev_metabolites['best_name'].str.startswith('low score', na=False)
rev_metabolites = rev_metabolites.sort_values(['well_annotated','vip_max'], ascending=[False,False])

print(f"\nEarly_specific metabolites: {len(early_metabolites)} (well-annotated: {early_metabolites['well_annotated'].sum()})")
print(f"Late_specific metabolites:  {len(late_metabolites)} (well-annotated: {late_metabolites['well_annotated'].sum()})")
print(f"Reversed metabolites:       {len(rev_metabolites)} (well-annotated: {rev_metabolites['well_annotated'].sum()})")

# ─── 7. 保存 ───
out_dir = '/sessions/charming-youthful-pasteur/mnt/TMAO_pathway_analysis/output'

# 個別CSV
all_ora.to_csv(f'{out_dir}/01_ontology_ORA_all_groups.csv', index=False)
ev_comparison.to_csv(f'{out_dir}/02_early_vs_late_ontology_comparison.csv', index=False)
early_metabolites.to_csv(f'{out_dir}/03_Early_specific_metabolites.csv', index=False)
late_metabolites.to_csv(f'{out_dir}/04_Late_specific_metabolites.csv', index=False)
rev_metabolites.to_csv(f'{out_dir}/05_Reversed_metabolites.csv', index=False)

print("\n✓ CSVs saved")
print("\n=== TOP RESULTS: Early vs Late Ontology差異 ===")
print(ev_comparison[['Ontology','Early_count','Late_count','Early_pct','Late_pct','p_value','FDR']].head(20).to_string())

