cd "d:\qq\qq的文件\ldmt_dmoea"
python -m experiments.run_ablation --problems FDA1,FDA3,FDA4 --tau 25 --seeds 0-29 --out results/raw --n-gen 300 *> results/raw/_log_ablation.txt 2>&1
python -m experiments.run_sensitivity --problems FDA3 --tau 25 --seeds 0-29 --out results/raw --n-gen 300 *> results/raw/_log_sensitivity.txt 2>&1
python -m experiments.run_supplementary --problems DMOP1,DMOP2 --tau 25 --seeds 0-29 --out results/raw --n-gen 300 *> results/raw/_log_supplementary.txt 2>&1
"BATCH2_DONE" | Out-File -Append results/raw/_done.txt
