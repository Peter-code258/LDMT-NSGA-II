import pandas as pd, glob
for f in sorted(glob.glob('results/raw_pilot/*.csv')):
    df = pd.read_csv(f)
    algo = df['algorithm'].iloc[0]
    per_seed = df.groupby('seed').agg(MIGD=('IGD','mean'), MHV=('HV','mean'), mig_t=('migration_runtime_s','sum'))
    print(f"{algo:14s} MIGD={per_seed['MIGD'].mean():.4f}+-{per_seed['MIGD'].std():.4f}  MHV={per_seed['MHV'].mean():.4f}  sum_mig_t={per_seed['mig_t'].mean():.4f}s")
