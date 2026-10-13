# v2.5 thesis analysis

- per-generation rows: 1,782,000
- distinct runs: 3,600
- files: 153
- code versions: v2.2, v2.4

## comparison summary [v24]

| problem | tau | algorithm | MIGD mean | MIGD std | MHV mean | post-change MIGD | mig_time_s | n |
|---|---|---|---|---|---|---|---|---|
| DMOP1 | 25 | pps | 0.0074 | 0.0011 | 0.5314 | 0.0117 | 0.01285 | 30 |
| DMOP1 | 25 | ar1_predictive | 0.0079 | 0.0012 | 0.5301 | 0.0170 | 0.01306 | 30 |
| DMOP1 | 25 | ldt_nsga2 (alpha=0.5) | 0.0088 | 0.0011 | 0.5277 | 0.0290 | 0.01564 | 30 |
| DMOP1 | 25 | nsga2_static | 0.0101 | 0.0012 | 0.5250 | 0.0495 | 0.00025 | 30 |
| DMOP1 | 25 | gmm_inspired | 0.0102 | 0.0011 | 0.5250 | 0.0514 | 0.18846 | 30 |
| DMOP1 | 25 | dnsga2_inspired | 0.0104 | 0.0011 | 0.5246 | 0.0511 | 0.00067 | 30 |
| DMOP2 | 25 | pps | 0.0069 | 0.0011 | 0.6352 | 0.0108 | 0.01296 | 30 |
| DMOP2 | 25 | ar1_predictive | 0.0073 | 0.0012 | 0.6333 | 0.0161 | 0.01316 | 30 |
| DMOP2 | 25 | ldt_nsga2 (alpha=0.5) | 0.0082 | 0.0011 | 0.6311 | 0.0276 | 0.01497 | 30 |
| DMOP2 | 25 | gmm_inspired | 0.0097 | 0.0011 | 0.6272 | 0.0495 | 0.13755 | 30 |
| DMOP2 | 25 | nsga2_static | 0.0097 | 0.0011 | 0.6272 | 0.0484 | 0.00022 | 30 |
| DMOP2 | 25 | dnsga2_inspired | 0.0098 | 0.0011 | 0.6267 | 0.0494 | 0.00068 | 30 |
| DMOP3 | 25 | nsga2_static | 0.1915 | 0.0002 | 1.3258 | 0.1915 | 0.00025 | 30 |
| DMOP3 | 25 | dnsga2_inspired | 0.1915 | 0.0002 | 1.3258 | 0.1468 | 0.00075 | 30 |
| DMOP3 | 25 | pps | 0.1920 | 0.0002 | 1.3257 | 0.1453 | 0.01468 | 30 |
| DMOP3 | 25 | ar1_predictive | 0.1920 | 0.0002 | 1.3257 | 0.1481 | 0.01389 | 30 |
| DMOP3 | 25 | ldt_nsga2 (alpha=0.5) | 0.1921 | 0.0010 | 1.3249 | 0.0608 | 0.01679 | 30 |
| DMOP3 | 25 | gmm_inspired | 0.1922 | 0.0006 | 1.3249 | 0.0596 | 0.13163 | 30 |
| FDA1 | 10 | pps | 0.0122 | 0.0014 | 0.8554 | 0.0165 | 0.03364 | 30 |
| FDA1 | 10 | ldt_nsga2 (alpha=0.5) | 0.0156 | 0.0010 | 0.8494 | 0.0289 | 0.04975 | 30 |
| FDA1 | 10 | ar1_predictive | 0.0173 | 0.0032 | 0.8463 | 0.0316 | 0.03396 | 30 |
| FDA1 | 10 | gmm_inspired | 0.0203 | 0.0010 | 0.8411 | 0.0418 | 0.49918 | 30 |
| FDA1 | 10 | nsga2_static | 0.0204 | 0.0012 | 0.8409 | 0.0391 | 0.00073 | 30 |
| FDA1 | 10 | dnsga2_inspired | 0.0211 | 0.0011 | 0.8396 | 0.0414 | 0.00197 | 30 |
| FDA1 | 25 | pps | 0.0062 | 0.0005 | 0.8668 | 0.0104 | 0.01275 | 30 |
| FDA1 | 25 | ar1_predictive | 0.0065 | 0.0005 | 0.8658 | 0.0129 | 0.01297 | 30 |
| FDA1 | 25 | ldt_nsga2 (alpha=0.5) | 0.0074 | 0.0005 | 0.8637 | 0.0238 | 0.01515 | 30 |
| FDA1 | 25 | gmm_inspired | 0.0084 | 0.0006 | 0.8617 | 0.0363 | 0.13383 | 30 |
| FDA1 | 25 | nsga2_static | 0.0086 | 0.0005 | 0.8612 | 0.0338 | 0.00034 | 30 |
| FDA1 | 25 | dnsga2_inspired | 0.0087 | 0.0005 | 0.8612 | 0.0347 | 0.00080 | 30 |
| FDA1 | 50 | pps | 0.0049 | 0.0001 | 0.8692 | 0.0112 | 0.00577 | 30 |
| FDA1 | 50 | ar1_predictive | 0.0050 | 0.0001 | 0.8690 | 0.0120 | 0.00580 | 30 |
| FDA1 | 50 | ldt_nsga2 (alpha=0.5) | 0.0052 | 0.0001 | 0.8682 | 0.0220 | 0.01116 | 30 |
| FDA1 | 50 | gmm_inspired | 0.0053 | 0.0001 | 0.8680 | 0.0316 | 0.04917 | 30 |
| FDA1 | 50 | nsga2_static | 0.0053 | 0.0001 | 0.8679 | 0.0290 | 0.00010 | 30 |
| FDA1 | 50 | dnsga2_inspired | 0.0054 | 0.0001 | 0.8678 | 0.0297 | 0.00031 | 30 |
| FDA3 | 10 | pps | 0.0132 | 0.0013 | 1.1160 | 0.0221 | 0.03406 | 30 |
| FDA3 | 10 | ldt_nsga2 (alpha=0.5) | 0.0177 | 0.0009 | 1.1064 | 0.0433 | 0.04406 | 30 |
| FDA3 | 10 | ar1_predictive | 0.0184 | 0.0026 | 1.1008 | 0.0403 | 0.03435 | 30 |
| FDA3 | 10 | nsga2_static | 0.0215 | 0.0013 | 1.0896 | 0.0416 | 0.00089 | 30 |
| FDA3 | 10 | gmm_inspired | 0.0220 | 0.0011 | 1.0914 | 0.0484 | 0.46476 | 30 |
| FDA3 | 10 | dnsga2_inspired | 0.0227 | 0.0011 | 1.0858 | 0.0447 | 0.00248 | 30 |
| FDA3 | 25 | pps | 0.0066 | 0.0005 | 1.5382 | 0.0129 | 0.01315 | 30 |
| FDA3 | 25 | ar1_predictive | 0.0070 | 0.0005 | 1.5362 | 0.0151 | 0.01327 | 30 |
| FDA3 | 25 | ldt_nsga2 (alpha=0.5) | 0.0079 | 0.0005 | 1.5319 | 0.0310 | 0.01570 | 30 |
| FDA3 | 25 | gmm_inspired | 0.0090 | 0.0006 | 1.5270 | 0.0395 | 0.12884 | 30 |
| FDA3 | 25 | nsga2_static | 0.0092 | 0.0005 | 1.5264 | 0.0355 | 0.00026 | 30 |
| FDA3 | 25 | dnsga2_inspired | 0.0092 | 0.0005 | 1.5257 | 0.0368 | 0.00071 | 30 |
| FDA3 | 50 | pps | 0.0051 | 0.0001 | 1.8233 | 0.0126 | 0.00578 | 30 |
| FDA3 | 50 | ar1_predictive | 0.0052 | 0.0001 | 1.8228 | 0.0135 | 0.00591 | 30 |
| FDA3 | 50 | ldt_nsga2 (alpha=0.5) | 0.0054 | 0.0001 | 1.8213 | 0.0255 | 0.00696 | 30 |
| FDA3 | 50 | gmm_inspired | 0.0055 | 0.0001 | 1.8206 | 0.0336 | 0.04918 | 30 |
| FDA3 | 50 | dnsga2_inspired | 0.0055 | 0.0001 | 1.8205 | 0.0307 | 0.00052 | 30 |
| FDA3 | 50 | nsga2_static | 0.0055 | 0.0001 | 1.8204 | 0.0296 | 0.00015 | 30 |
| FDA4 | 10 | pps | 0.1167 | 0.0075 | 0.5976 | 0.1359 | 0.04116 | 30 |
| FDA4 | 10 | ldt_nsga2 (alpha=0.5) | 0.1236 | 0.0047 | 0.5707 | 0.1606 | 0.04257 | 30 |
| FDA4 | 10 | ar1_predictive | 0.1507 | 0.0112 | 0.5210 | 0.2140 | 0.03666 | 30 |
| FDA4 | 10 | gmm_inspired | 0.1595 | 0.0090 | 0.4887 | 0.2213 | 0.46729 | 30 |
| FDA4 | 10 | dnsga2_inspired | 0.1628 | 0.0096 | 0.4855 | 0.2198 | 0.00216 | 30 |
| FDA4 | 10 | nsga2_static | 0.2101 | 0.0161 | 0.4070 | 0.2857 | 0.00076 | 30 |
| FDA4 | 25 | pps | 0.0881 | 0.0046 | 0.6488 | 0.1074 | 0.01401 | 30 |
| FDA4 | 25 | ldt_nsga2 (alpha=0.5) | 0.0904 | 0.0039 | 0.6343 | 0.1304 | 0.01611 | 30 |
| FDA4 | 25 | ar1_predictive | 0.0987 | 0.0082 | 0.6132 | 0.1721 | 0.01493 | 30 |
| FDA4 | 25 | gmm_inspired | 0.1011 | 0.0044 | 0.6004 | 0.1739 | 0.12942 | 30 |
| FDA4 | 25 | nsga2_static | 0.1020 | 0.0035 | 0.6024 | 0.1682 | 0.00028 | 30 |
| FDA4 | 25 | dnsga2_inspired | 0.1027 | 0.0044 | 0.6004 | 0.1750 | 0.00078 | 30 |
| FDA4 | 50 | pps | 0.0736 | 0.0019 | 0.6863 | 0.0855 | 0.00635 | 30 |
| FDA4 | 50 | ldt_nsga2 (alpha=0.5) | 0.0745 | 0.0017 | 0.6810 | 0.1073 | 0.00763 | 30 |
| FDA4 | 50 | ar1_predictive | 0.0746 | 0.0018 | 0.6807 | 0.1002 | 0.00631 | 30 |
| FDA4 | 50 | gmm_inspired | 0.0757 | 0.0018 | 0.6748 | 0.1240 | 0.05359 | 30 |
| FDA4 | 50 | dnsga2_inspired | 0.0759 | 0.0022 | 0.6753 | 0.1224 | 0.00039 | 30 |
| FDA4 | 50 | nsga2_static | 0.0764 | 0.0020 | 0.6746 | 0.1173 | 0.00014 | 30 |

### Wilcoxon (paired) + Holm, ref = ldt_nsga2 [v24]

| problem | tau | contrast | n | mean ref | mean other | p_raw | p_holm | better |
|---|---|---|---|---|---|---|---|---|
| DMOP1 | 25 | ldt_nsga2 vs nsga2_static | 30 | 0.0088 | 0.0101 | 1.863e-09 | 9.313e-09 | ldt_nsga2 |
| DMOP1 | 25 | ldt_nsga2 vs dnsga2_inspired | 30 | 0.0088 | 0.0104 | 1.863e-09 | 9.313e-09 | ldt_nsga2 |
| DMOP1 | 25 | ldt_nsga2 vs ar1_predictive | 30 | 0.0088 | 0.0079 | 1.863e-09 | 9.313e-09 | ar1_predictive |
| DMOP1 | 25 | ldt_nsga2 vs gmm_inspired | 30 | 0.0088 | 0.0102 | 1.863e-09 | 9.313e-09 | ldt_nsga2 |
| DMOP1 | 25 | ldt_nsga2 vs pps | 30 | 0.0088 | 0.0074 | 1.863e-09 | 9.313e-09 | pps |
| DMOP2 | 25 | ldt_nsga2 vs nsga2_static | 30 | 0.0082 | 0.0097 | 1.863e-09 | 9.313e-09 | ldt_nsga2 |
| DMOP2 | 25 | ldt_nsga2 vs dnsga2_inspired | 30 | 0.0082 | 0.0098 | 1.863e-09 | 9.313e-09 | ldt_nsga2 |
| DMOP2 | 25 | ldt_nsga2 vs ar1_predictive | 30 | 0.0082 | 0.0073 | 1.863e-09 | 9.313e-09 | ar1_predictive |
| DMOP2 | 25 | ldt_nsga2 vs gmm_inspired | 30 | 0.0082 | 0.0097 | 1.863e-09 | 9.313e-09 | ldt_nsga2 |
| DMOP2 | 25 | ldt_nsga2 vs pps | 30 | 0.0082 | 0.0069 | 1.863e-09 | 9.313e-09 | pps |
| DMOP3 | 25 | ldt_nsga2 vs nsga2_static | 30 | 0.1921 | 0.1915 | 1.863e-09 | 9.313e-09 | nsga2_static |
| DMOP3 | 25 | ldt_nsga2 vs dnsga2_inspired | 30 | 0.1921 | 0.1915 | 1.863e-09 | 9.313e-09 | dnsga2_inspired |
| DMOP3 | 25 | ldt_nsga2 vs gmm_inspired | 30 | 0.1921 | 0.1922 | 1.060e-02 | 3.179e-02 | ldt_nsga2 |
| DMOP3 | 25 | ldt_nsga2 vs ar1_predictive | 30 | 0.1921 | 0.1920 | 2.534e-01 | 5.069e-01 | ar1_predictive |
| DMOP3 | 25 | ldt_nsga2 vs pps | 30 | 0.1921 | 0.1920 | 3.085e-01 | 5.069e-01 | pps |
| FDA1 | 10 | ldt_nsga2 vs nsga2_static | 30 | 0.0156 | 0.0204 | 1.863e-09 | 9.313e-09 | ldt_nsga2 |
| FDA1 | 10 | ldt_nsga2 vs dnsga2_inspired | 30 | 0.0156 | 0.0211 | 1.863e-09 | 9.313e-09 | ldt_nsga2 |
| FDA1 | 10 | ldt_nsga2 vs gmm_inspired | 30 | 0.0156 | 0.0203 | 1.863e-09 | 9.313e-09 | ldt_nsga2 |
| FDA1 | 10 | ldt_nsga2 vs pps | 30 | 0.0156 | 0.0122 | 1.863e-09 | 9.313e-09 | pps |
| FDA1 | 10 | ldt_nsga2 vs ar1_predictive | 30 | 0.0156 | 0.0173 | 1.341e-03 | 1.341e-03 | ldt_nsga2 |
| FDA1 | 25 | ldt_nsga2 vs nsga2_static | 30 | 0.0074 | 0.0086 | 1.863e-09 | 9.313e-09 | ldt_nsga2 |
| FDA1 | 25 | ldt_nsga2 vs dnsga2_inspired | 30 | 0.0074 | 0.0087 | 1.863e-09 | 9.313e-09 | ldt_nsga2 |
| FDA1 | 25 | ldt_nsga2 vs ar1_predictive | 30 | 0.0074 | 0.0065 | 1.863e-09 | 9.313e-09 | ar1_predictive |
| FDA1 | 25 | ldt_nsga2 vs gmm_inspired | 30 | 0.0074 | 0.0084 | 1.863e-09 | 9.313e-09 | ldt_nsga2 |
| FDA1 | 25 | ldt_nsga2 vs pps | 30 | 0.0074 | 0.0062 | 1.863e-09 | 9.313e-09 | pps |
| FDA1 | 50 | ldt_nsga2 vs pps | 30 | 0.0052 | 0.0049 | 3.725e-09 | 1.863e-08 | pps |
| FDA1 | 50 | ldt_nsga2 vs ar1_predictive | 30 | 0.0052 | 0.0050 | 9.313e-09 | 3.725e-08 | ar1_predictive |
| FDA1 | 50 | ldt_nsga2 vs dnsga2_inspired | 30 | 0.0052 | 0.0054 | 8.326e-07 | 2.498e-06 | ldt_nsga2 |
| FDA1 | 50 | ldt_nsga2 vs nsga2_static | 30 | 0.0052 | 0.0053 | 1.106e-04 | 2.212e-04 | ldt_nsga2 |
| FDA1 | 50 | ldt_nsga2 vs gmm_inspired | 30 | 0.0052 | 0.0053 | 1.130e-02 | 1.130e-02 | ldt_nsga2 |
| FDA3 | 10 | ldt_nsga2 vs nsga2_static | 30 | 0.0177 | 0.0215 | 1.863e-09 | 9.313e-09 | ldt_nsga2 |
| FDA3 | 10 | ldt_nsga2 vs dnsga2_inspired | 30 | 0.0177 | 0.0227 | 1.863e-09 | 9.313e-09 | ldt_nsga2 |
| FDA3 | 10 | ldt_nsga2 vs gmm_inspired | 30 | 0.0177 | 0.0220 | 1.863e-09 | 9.313e-09 | ldt_nsga2 |
| FDA3 | 10 | ldt_nsga2 vs pps | 30 | 0.0177 | 0.0132 | 1.863e-09 | 9.313e-09 | pps |
| FDA3 | 10 | ldt_nsga2 vs ar1_predictive | 30 | 0.0177 | 0.0184 | 1.981e-01 | 1.981e-01 | ldt_nsga2 |
| FDA3 | 25 | ldt_nsga2 vs nsga2_static | 30 | 0.0079 | 0.0092 | 1.863e-09 | 9.313e-09 | ldt_nsga2 |
| FDA3 | 25 | ldt_nsga2 vs dnsga2_inspired | 30 | 0.0079 | 0.0092 | 1.863e-09 | 9.313e-09 | ldt_nsga2 |
| FDA3 | 25 | ldt_nsga2 vs ar1_predictive | 30 | 0.0079 | 0.0070 | 1.863e-09 | 9.313e-09 | ar1_predictive |
| FDA3 | 25 | ldt_nsga2 vs gmm_inspired | 30 | 0.0079 | 0.0090 | 1.863e-09 | 9.313e-09 | ldt_nsga2 |
| FDA3 | 25 | ldt_nsga2 vs pps | 30 | 0.0079 | 0.0066 | 1.863e-09 | 9.313e-09 | pps |
| FDA3 | 50 | ldt_nsga2 vs pps | 30 | 0.0054 | 0.0051 | 1.863e-09 | 9.313e-09 | pps |
| FDA3 | 50 | ldt_nsga2 vs ar1_predictive | 30 | 0.0054 | 0.0052 | 3.725e-09 | 1.490e-08 | ar1_predictive |
| FDA3 | 50 | ldt_nsga2 vs dnsga2_inspired | 30 | 0.0054 | 0.0055 | 1.684e-06 | 5.051e-06 | ldt_nsga2 |
| FDA3 | 50 | ldt_nsga2 vs nsga2_static | 30 | 0.0054 | 0.0055 | 3.239e-06 | 6.478e-06 | ldt_nsga2 |
| FDA3 | 50 | ldt_nsga2 vs gmm_inspired | 30 | 0.0054 | 0.0055 | 7.911e-05 | 7.911e-05 | ldt_nsga2 |
| FDA4 | 10 | ldt_nsga2 vs nsga2_static | 30 | 0.1236 | 0.2101 | 1.863e-09 | 9.313e-09 | ldt_nsga2 |
| FDA4 | 10 | ldt_nsga2 vs dnsga2_inspired | 30 | 0.1236 | 0.1628 | 1.863e-09 | 9.313e-09 | ldt_nsga2 |
| FDA4 | 10 | ldt_nsga2 vs ar1_predictive | 30 | 0.1236 | 0.1507 | 1.863e-09 | 9.313e-09 | ldt_nsga2 |
| FDA4 | 10 | ldt_nsga2 vs gmm_inspired | 30 | 0.1236 | 0.1595 | 1.863e-09 | 9.313e-09 | ldt_nsga2 |
| FDA4 | 10 | ldt_nsga2 vs pps | 30 | 0.1236 | 0.1167 | 3.148e-07 | 3.148e-07 | pps |
| FDA4 | 25 | ldt_nsga2 vs nsga2_static | 30 | 0.0904 | 0.1020 | 1.863e-09 | 9.313e-09 | ldt_nsga2 |
| FDA4 | 25 | ldt_nsga2 vs dnsga2_inspired | 30 | 0.0904 | 0.1027 | 1.863e-09 | 9.313e-09 | ldt_nsga2 |
| FDA4 | 25 | ldt_nsga2 vs ar1_predictive | 30 | 0.0904 | 0.0987 | 1.863e-09 | 9.313e-09 | ldt_nsga2 |
| FDA4 | 25 | ldt_nsga2 vs gmm_inspired | 30 | 0.0904 | 0.1011 | 1.863e-09 | 9.313e-09 | ldt_nsga2 |
| FDA4 | 25 | ldt_nsga2 vs pps | 30 | 0.0904 | 0.0881 | 3.239e-06 | 3.239e-06 | pps |
| FDA4 | 50 | ldt_nsga2 vs nsga2_static | 30 | 0.0745 | 0.0764 | 4.408e-05 | 2.204e-04 | ldt_nsga2 |
| FDA4 | 50 | ldt_nsga2 vs dnsga2_inspired | 30 | 0.0745 | 0.0759 | 1.374e-04 | 5.496e-04 | ldt_nsga2 |
| FDA4 | 50 | ldt_nsga2 vs gmm_inspired | 30 | 0.0745 | 0.0757 | 6.084e-04 | 1.825e-03 | ldt_nsga2 |
| FDA4 | 50 | ldt_nsga2 vs pps | 30 | 0.0745 | 0.0736 | 6.640e-03 | 1.328e-02 | pps |
| FDA4 | 50 | ldt_nsga2 vs ar1_predictive | 30 | 0.0745 | 0.0746 | 5.158e-01 | 5.158e-01 | ldt_nsga2 |

### Friedman + Nemenyi [v24]

| problem | tau | k | N | chi2 | p | CD | avg ranks |
|---|---|---|---|---|---|---|---|
| DMOP1 | 25 | 6 | 30 | 136.59 | 9.49e-28 | 1.3767 | nsga2_static=4.67, dnsga2_inspired=5.53, ar1_predictive=2.00, gmm_inspired=4.80, pps=1.00, ldt_nsga2=3.00 |
| DMOP2 | 25 | 6 | 30 | 133.79 | 3.73e-27 | 1.3767 | nsga2_static=4.90, dnsga2_inspired=5.27, ar1_predictive=2.00, gmm_inspired=4.83, pps=1.00, ldt_nsga2=3.00 |
| DMOP3 | 25 | 6 | 30 | 109.05 | 6.51e-22 | 1.3767 | nsga2_static=1.50, dnsga2_inspired=1.50, ar1_predictive=4.50, gmm_inspired=5.17, pps=4.33, ldt_nsga2=4.00 |
| FDA1 | 10 | 6 | 30 | 124.88 | 2.91e-25 | 1.3767 | nsga2_static=4.63, dnsga2_inspired=5.60, ar1_predictive=3.00, gmm_inspired=4.47, pps=1.00, ldt_nsga2=2.30 |
| FDA1 | 25 | 6 | 30 | 138.59 | 3.57e-28 | 1.3767 | nsga2_static=5.30, dnsga2_inspired=5.37, ar1_predictive=2.00, gmm_inspired=4.33, pps=1.00, ldt_nsga2=3.00 |
| FDA1 | 50 | 6 | 30 | 114.97 | 3.64e-23 | 1.3767 | nsga2_static=4.90, dnsga2_inspired=5.27, ar1_predictive=1.80, gmm_inspired=4.17, pps=1.27, ldt_nsga2=3.60 |
| FDA3 | 10 | 6 | 30 | 127.41 | 8.44e-26 | 1.3767 | nsga2_static=4.37, dnsga2_inspired=5.60, ar1_predictive=2.77, gmm_inspired=4.83, pps=1.00, ldt_nsga2=2.43 |
| FDA3 | 25 | 6 | 30 | 139.39 | 2.41e-28 | 1.3767 | nsga2_static=5.03, dnsga2_inspired=5.60, ar1_predictive=2.00, gmm_inspired=4.37, pps=1.00, ldt_nsga2=3.00 |
| FDA3 | 50 | 6 | 30 | 113.56 | 7.23e-23 | 1.3767 | nsga2_static=4.87, dnsga2_inspired=4.87, ar1_predictive=1.77, gmm_inspired=4.73, pps=1.27, ldt_nsga2=3.50 |
| FDA4 | 10 | 6 | 30 | 138.00 | 4.76e-28 | 1.3767 | nsga2_static=6.00, dnsga2_inspired=4.50, ar1_predictive=3.30, gmm_inspired=4.20, pps=1.10, ldt_nsga2=1.90 |
| FDA4 | 25 | 6 | 30 | 110.55 | 3.13e-22 | 1.3767 | nsga2_static=4.63, dnsga2_inspired=5.07, ar1_predictive=4.00, gmm_inspired=4.30, pps=1.13, ldt_nsga2=1.87 |
| FDA4 | 50 | 6 | 30 | 64.23 | 1.62e-12 | 1.3767 | nsga2_static=4.80, dnsga2_inspired=4.50, ar1_predictive=3.10, gmm_inspired=4.27, pps=1.67, ldt_nsga2=2.67 |

overall mean rank (all cells averaged): pps=1.397, ar1_predictive=2.686, ldt_nsga2=2.856, gmm_inspired=4.539, nsga2_static=4.633, dnsga2_inspired=4.889

## comparison summary [v25]

| problem | tau | algorithm | MIGD mean | MIGD std | MHV mean | post-change MIGD | mig_time_s | n |
|---|---|---|---|---|---|---|---|---|
| DMOP1A | 25 | ldt_alpha1 (alpha=1.0) | 0.0570 | 0.0171 | 0.5703 | 0.0634 | 0.01747 | 30 |
| DMOP1A | 25 | gmm_inspired | 0.0571 | 0.0166 | 0.5699 | 0.0676 | 0.22772 | 30 |
| DMOP1A | 25 | ar1_predictive | 0.0575 | 0.0184 | 0.5697 | 0.0646 | 0.01813 | 30 |
| DMOP1A | 25 | pps | 0.0575 | 0.0184 | 0.5697 | 0.0636 | 0.01763 | 30 |
| DMOP1A | 25 | dnsga2_inspired | 0.0578 | 0.0189 | 0.5693 | 0.0634 | 0.00127 | 30 |
| DMOP1A | 25 | nsga2_static | 0.0583 | 0.0200 | 0.5688 | 0.0636 | 0.00073 | 30 |
| DMOP2A | 25 | ldt_alpha1 (alpha=1.0) | 0.0265 | 0.0162 | 0.5629 | 0.0679 | 0.01972 | 30 |
| DMOP2A | 25 | pps | 0.0267 | 0.0162 | 0.5622 | 0.0615 | 0.01959 | 30 |
| DMOP2A | 25 | ar1_predictive | 0.0284 | 0.0164 | 0.5589 | 0.0823 | 0.02127 | 30 |
| DMOP2A | 25 | gmm_inspired | 0.0452 | 0.0173 | 0.5279 | 0.3448 | 0.25492 | 30 |
| DMOP2A | 25 | nsga2_static | 0.0465 | 0.0162 | 0.5251 | 0.3556 | 0.00072 | 30 |
| DMOP2A | 25 | dnsga2_inspired | 0.0477 | 0.0167 | 0.5229 | 0.3633 | 0.00139 | 30 |
| DMOP3A | 25 | pps | 0.0479 | 0.0096 | 0.8237 | 0.4221 | 0.01793 | 30 |
| DMOP3A | 25 | ldt_alpha1 (alpha=1.0) | 0.0525 | 0.0135 | 0.8198 | 0.4013 | 0.01735 | 30 |
| DMOP3A | 25 | ar1_predictive | 0.0567 | 0.0126 | 0.8129 | 0.7182 | 0.01852 | 30 |
| DMOP3A | 25 | nsga2_static | 0.0662 | 0.0113 | 0.8008 | 0.4940 | 0.00074 | 30 |
| DMOP3A | 25 | gmm_inspired | 0.0674 | 0.0114 | 0.8019 | 0.4901 | 0.24119 | 30 |
| DMOP3A | 25 | dnsga2_inspired | 0.0697 | 0.0105 | 0.7991 | 0.5025 | 0.00125 | 30 |
| FDA1 | 10 | ldt_alpha1 (alpha=1.0) | 0.0120 | 0.0009 | 0.8561 | 0.0212 | 0.07706 | 30 |
| FDA1 | 10 | pps | 0.0122 | 0.0014 | 0.8554 | 0.0165 | 0.09353 | 30 |
| FDA1 | 10 | ar1_predictive | 0.0173 | 0.0032 | 0.8463 | 0.0316 | 0.08299 | 30 |
| FDA1 | 10 | gmm_inspired | 0.0203 | 0.0010 | 0.8411 | 0.0418 | 0.83747 | 30 |
| FDA1 | 10 | nsga2_static | 0.0204 | 0.0012 | 0.8409 | 0.0391 | 0.00282 | 30 |
| FDA1 | 10 | dnsga2_inspired | 0.0211 | 0.0011 | 0.8396 | 0.0414 | 0.00620 | 30 |
| FDA1 | 25 | pps | 0.0062 | 0.0005 | 0.8668 | 0.0104 | 0.03787 | 30 |
| FDA1 | 25 | ldt_alpha1 (alpha=1.0) | 0.0062 | 0.0005 | 0.8667 | 0.0174 | 0.03014 | 30 |
| FDA1 | 25 | ar1_predictive | 0.0065 | 0.0005 | 0.8658 | 0.0129 | 0.03250 | 30 |
| FDA1 | 25 | gmm_inspired | 0.0084 | 0.0006 | 0.8617 | 0.0363 | 0.37408 | 30 |
| FDA1 | 25 | nsga2_static | 0.0086 | 0.0005 | 0.8612 | 0.0338 | 0.00103 | 30 |
| FDA1 | 25 | dnsga2_inspired | 0.0087 | 0.0005 | 0.8612 | 0.0347 | 0.00237 | 30 |
| FDA1 | 50 | pps | 0.0049 | 0.0001 | 0.8692 | 0.0112 | 0.01745 | 30 |
| FDA1 | 50 | ldt_alpha1 (alpha=1.0) | 0.0049 | 0.0001 | 0.8692 | 0.0183 | 0.01381 | 30 |
| FDA1 | 50 | ar1_predictive | 0.0050 | 0.0001 | 0.8690 | 0.0120 | 0.01380 | 30 |
| FDA1 | 50 | gmm_inspired | 0.0053 | 0.0001 | 0.8680 | 0.0316 | 0.21256 | 30 |
| FDA1 | 50 | nsga2_static | 0.0053 | 0.0001 | 0.8679 | 0.0290 | 0.00050 | 30 |
| FDA1 | 50 | dnsga2_inspired | 0.0054 | 0.0001 | 0.8678 | 0.0297 | 0.00108 | 30 |
| FDA3 | 10 | pps | 0.0132 | 0.0013 | 1.1160 | 0.0221 | 0.09872 | 30 |
| FDA3 | 10 | ldt_alpha1 (alpha=1.0) | 0.0136 | 0.0009 | 1.1184 | 0.0346 | 0.07782 | 30 |
| FDA3 | 10 | ar1_predictive | 0.0184 | 0.0026 | 1.1008 | 0.0403 | 0.08672 | 30 |
| FDA3 | 10 | nsga2_static | 0.0215 | 0.0013 | 1.0896 | 0.0416 | 0.00276 | 30 |
| FDA3 | 10 | gmm_inspired | 0.0220 | 0.0011 | 1.0914 | 0.0484 | 0.82136 | 30 |
| FDA3 | 10 | dnsga2_inspired | 0.0227 | 0.0011 | 1.0858 | 0.0447 | 0.00593 | 30 |
| FDA3 | 25 | pps | 0.0066 | 0.0005 | 1.5382 | 0.0129 | 0.03772 | 30 |
| FDA3 | 25 | ldt_alpha1 (alpha=1.0) | 0.0066 | 0.0005 | 1.5378 | 0.0240 | 0.03079 | 30 |
| FDA3 | 25 | ar1_predictive | 0.0070 | 0.0005 | 1.5362 | 0.0151 | 0.03196 | 30 |
| FDA3 | 25 | gmm_inspired | 0.0090 | 0.0006 | 1.5270 | 0.0395 | 0.37189 | 30 |
| FDA3 | 25 | nsga2_static | 0.0092 | 0.0005 | 1.5264 | 0.0355 | 0.00107 | 30 |
| FDA3 | 25 | dnsga2_inspired | 0.0092 | 0.0005 | 1.5257 | 0.0368 | 0.00228 | 30 |
| FDA3 | 50 | pps | 0.0051 | 0.0001 | 1.8233 | 0.0126 | 0.01718 | 30 |
| FDA3 | 50 | ldt_alpha1 (alpha=1.0) | 0.0052 | 0.0001 | 1.8230 | 0.0212 | 0.01367 | 30 |
| FDA3 | 50 | ar1_predictive | 0.0052 | 0.0001 | 1.8228 | 0.0135 | 0.01380 | 30 |
| FDA3 | 50 | gmm_inspired | 0.0055 | 0.0001 | 1.8206 | 0.0336 | 0.20735 | 30 |
| FDA3 | 50 | dnsga2_inspired | 0.0055 | 0.0001 | 1.8205 | 0.0307 | 0.00110 | 30 |
| FDA3 | 50 | nsga2_static | 0.0055 | 0.0001 | 1.8204 | 0.0296 | 0.00049 | 30 |
| FDA4 | 10 | ldt_alpha1 (alpha=1.0) | 0.1068 | 0.0050 | 0.6188 | 0.1302 | 0.06129 | 30 |
| FDA4 | 10 | pps | 0.1167 | 0.0075 | 0.5976 | 0.1359 | 0.08300 | 30 |
| FDA4 | 10 | ar1_predictive | 0.1507 | 0.0112 | 0.5210 | 0.2140 | 0.09277 | 30 |
| FDA4 | 10 | gmm_inspired | 0.1595 | 0.0090 | 0.4887 | 0.2213 | 0.84031 | 30 |
| FDA4 | 10 | dnsga2_inspired | 0.1628 | 0.0096 | 0.4855 | 0.2198 | 0.00695 | 30 |
| FDA4 | 10 | nsga2_static | 0.2101 | 0.0161 | 0.4070 | 0.2857 | 0.00312 | 30 |
| FDA4 | 25 | ldt_alpha1 (alpha=1.0) | 0.0858 | 0.0043 | 0.6527 | 0.1119 | 0.02235 | 30 |
| FDA4 | 25 | pps | 0.0881 | 0.0046 | 0.6488 | 0.1074 | 0.03131 | 30 |
| FDA4 | 25 | ar1_predictive | 0.0987 | 0.0082 | 0.6132 | 0.1721 | 0.03783 | 30 |
| FDA4 | 25 | gmm_inspired | 0.1011 | 0.0044 | 0.6004 | 0.1739 | 0.37924 | 30 |
| FDA4 | 25 | nsga2_static | 0.1020 | 0.0035 | 0.6024 | 0.1682 | 0.00126 | 30 |
| FDA4 | 25 | dnsga2_inspired | 0.1027 | 0.0044 | 0.6004 | 0.1750 | 0.00270 | 30 |
| FDA4 | 50 | pps | 0.0736 | 0.0019 | 0.6863 | 0.0855 | 0.01423 | 30 |
| FDA4 | 50 | ldt_alpha1 (alpha=1.0) | 0.0737 | 0.0022 | 0.6850 | 0.0990 | 0.00970 | 30 |
| FDA4 | 50 | ar1_predictive | 0.0746 | 0.0018 | 0.6807 | 0.1002 | 0.01830 | 30 |
| FDA4 | 50 | gmm_inspired | 0.0757 | 0.0018 | 0.6748 | 0.1240 | 0.21914 | 30 |
| FDA4 | 50 | dnsga2_inspired | 0.0759 | 0.0022 | 0.6753 | 0.1224 | 0.00119 | 30 |
| FDA4 | 50 | nsga2_static | 0.0764 | 0.0020 | 0.6746 | 0.1173 | 0.00059 | 30 |

### Wilcoxon (paired) + Holm, ref = ldt_alpha1 [v25]

| problem | tau | contrast | n | mean ref | mean other | p_raw | p_holm | better |
|---|---|---|---|---|---|---|---|---|
| DMOP1A | 25 | ldt_alpha1 vs dnsga2_inspired | 30 | 0.0570 | 0.0578 | 6.084e-04 | 3.042e-03 | ldt_alpha1 |
| DMOP1A | 25 | ldt_alpha1 vs gmm_inspired | 30 | 0.0570 | 0.0571 | 1.457e-03 | 5.830e-03 | ldt_alpha1 |
| DMOP1A | 25 | ldt_alpha1 vs nsga2_static | 30 | 0.0570 | 0.0583 | 2.987e-03 | 8.961e-03 | ldt_alpha1 |
| DMOP1A | 25 | ldt_alpha1 vs ar1_predictive | 30 | 0.0570 | 0.0575 | 4.491e-02 | 5.548e-02 | ldt_alpha1 |
| DMOP1A | 25 | ldt_alpha1 vs pps | 30 | 0.0570 | 0.0575 | 2.774e-02 | 5.548e-02 | ldt_alpha1 |
| DMOP2A | 25 | ldt_alpha1 vs nsga2_static | 30 | 0.0265 | 0.0465 | 1.863e-09 | 9.313e-09 | ldt_alpha1 |
| DMOP2A | 25 | ldt_alpha1 vs dnsga2_inspired | 30 | 0.0265 | 0.0477 | 1.863e-09 | 9.313e-09 | ldt_alpha1 |
| DMOP2A | 25 | ldt_alpha1 vs gmm_inspired | 30 | 0.0265 | 0.0452 | 1.863e-09 | 9.313e-09 | ldt_alpha1 |
| DMOP2A | 25 | ldt_alpha1 vs ar1_predictive | 30 | 0.0265 | 0.0284 | 3.239e-06 | 6.478e-06 | ldt_alpha1 |
| DMOP2A | 25 | ldt_alpha1 vs pps | 30 | 0.0265 | 0.0267 | 7.303e-01 | 7.303e-01 | ldt_alpha1 |
| DMOP3A | 25 | ldt_alpha1 vs dnsga2_inspired | 30 | 0.0525 | 0.0697 | 3.239e-06 | 1.620e-05 | ldt_alpha1 |
| DMOP3A | 25 | ldt_alpha1 vs nsga2_static | 30 | 0.0525 | 0.0662 | 1.218e-05 | 4.873e-05 | ldt_alpha1 |
| DMOP3A | 25 | ldt_alpha1 vs gmm_inspired | 30 | 0.0525 | 0.0674 | 2.080e-05 | 6.239e-05 | ldt_alpha1 |
| DMOP3A | 25 | ldt_alpha1 vs pps | 30 | 0.0525 | 0.0479 | 1.048e-01 | 2.097e-01 | pps |
| DMOP3A | 25 | ldt_alpha1 vs ar1_predictive | 30 | 0.0525 | 0.0567 | 2.534e-01 | 2.534e-01 | ldt_alpha1 |
| FDA1 | 10 | ldt_alpha1 vs nsga2_static | 30 | 0.0120 | 0.0204 | 1.863e-09 | 9.313e-09 | ldt_alpha1 |
| FDA1 | 10 | ldt_alpha1 vs dnsga2_inspired | 30 | 0.0120 | 0.0211 | 1.863e-09 | 9.313e-09 | ldt_alpha1 |
| FDA1 | 10 | ldt_alpha1 vs ar1_predictive | 30 | 0.0120 | 0.0173 | 1.863e-09 | 9.313e-09 | ldt_alpha1 |
| FDA1 | 10 | ldt_alpha1 vs gmm_inspired | 30 | 0.0120 | 0.0203 | 1.863e-09 | 9.313e-09 | ldt_alpha1 |
| FDA1 | 10 | ldt_alpha1 vs pps | 30 | 0.0120 | 0.0122 | 6.989e-02 | 6.989e-02 | ldt_alpha1 |
| FDA1 | 25 | ldt_alpha1 vs nsga2_static | 30 | 0.0062 | 0.0086 | 1.863e-09 | 9.313e-09 | ldt_alpha1 |
| FDA1 | 25 | ldt_alpha1 vs dnsga2_inspired | 30 | 0.0062 | 0.0087 | 1.863e-09 | 9.313e-09 | ldt_alpha1 |
| FDA1 | 25 | ldt_alpha1 vs ar1_predictive | 30 | 0.0062 | 0.0065 | 1.863e-09 | 9.313e-09 | ldt_alpha1 |
| FDA1 | 25 | ldt_alpha1 vs gmm_inspired | 30 | 0.0062 | 0.0084 | 1.863e-09 | 9.313e-09 | ldt_alpha1 |
| FDA1 | 25 | ldt_alpha1 vs pps | 30 | 0.0062 | 0.0062 | 5.158e-01 | 5.158e-01 | pps |
| FDA1 | 50 | ldt_alpha1 vs nsga2_static | 30 | 0.0049 | 0.0053 | 1.863e-09 | 9.313e-09 | ldt_alpha1 |
| FDA1 | 50 | ldt_alpha1 vs dnsga2_inspired | 30 | 0.0049 | 0.0054 | 1.863e-09 | 9.313e-09 | ldt_alpha1 |
| FDA1 | 50 | ldt_alpha1 vs gmm_inspired | 30 | 0.0049 | 0.0053 | 1.863e-09 | 9.313e-09 | ldt_alpha1 |
| FDA1 | 50 | ldt_alpha1 vs ar1_predictive | 30 | 0.0049 | 0.0050 | 5.226e-02 | 1.045e-01 | ldt_alpha1 |
| FDA1 | 50 | ldt_alpha1 vs pps | 30 | 0.0049 | 0.0049 | 6.120e-01 | 6.120e-01 | pps |
| FDA3 | 10 | ldt_alpha1 vs nsga2_static | 30 | 0.0136 | 0.0215 | 1.863e-09 | 9.313e-09 | ldt_alpha1 |
| FDA3 | 10 | ldt_alpha1 vs dnsga2_inspired | 30 | 0.0136 | 0.0227 | 1.863e-09 | 9.313e-09 | ldt_alpha1 |
| FDA3 | 10 | ldt_alpha1 vs ar1_predictive | 30 | 0.0136 | 0.0184 | 1.863e-09 | 9.313e-09 | ldt_alpha1 |
| FDA3 | 10 | ldt_alpha1 vs gmm_inspired | 30 | 0.0136 | 0.0220 | 1.863e-09 | 9.313e-09 | ldt_alpha1 |
| FDA3 | 10 | ldt_alpha1 vs pps | 30 | 0.0136 | 0.0132 | 2.563e-04 | 2.563e-04 | pps |
| FDA3 | 25 | ldt_alpha1 vs nsga2_static | 30 | 0.0066 | 0.0092 | 1.863e-09 | 9.313e-09 | ldt_alpha1 |
| FDA3 | 25 | ldt_alpha1 vs dnsga2_inspired | 30 | 0.0066 | 0.0092 | 1.863e-09 | 9.313e-09 | ldt_alpha1 |
| FDA3 | 25 | ldt_alpha1 vs ar1_predictive | 30 | 0.0066 | 0.0070 | 1.863e-09 | 9.313e-09 | ldt_alpha1 |
| FDA3 | 25 | ldt_alpha1 vs gmm_inspired | 30 | 0.0066 | 0.0090 | 1.863e-09 | 9.313e-09 | ldt_alpha1 |
| FDA3 | 25 | ldt_alpha1 vs pps | 30 | 0.0066 | 0.0066 | 2.987e-03 | 2.987e-03 | pps |
| FDA3 | 50 | ldt_alpha1 vs nsga2_static | 30 | 0.0052 | 0.0055 | 1.863e-09 | 9.313e-09 | ldt_alpha1 |
| FDA3 | 50 | ldt_alpha1 vs dnsga2_inspired | 30 | 0.0052 | 0.0055 | 1.863e-09 | 9.313e-09 | ldt_alpha1 |
| FDA3 | 50 | ldt_alpha1 vs gmm_inspired | 30 | 0.0052 | 0.0055 | 1.863e-09 | 9.313e-09 | ldt_alpha1 |
| FDA3 | 50 | ldt_alpha1 vs pps | 30 | 0.0052 | 0.0051 | 3.272e-02 | 6.544e-02 | pps |
| FDA3 | 50 | ldt_alpha1 vs ar1_predictive | 30 | 0.0052 | 0.0052 | 3.184e-01 | 3.184e-01 | ldt_alpha1 |
| FDA4 | 10 | ldt_alpha1 vs nsga2_static | 30 | 0.1068 | 0.2101 | 1.863e-09 | 9.313e-09 | ldt_alpha1 |
| FDA4 | 10 | ldt_alpha1 vs dnsga2_inspired | 30 | 0.1068 | 0.1628 | 1.863e-09 | 9.313e-09 | ldt_alpha1 |
| FDA4 | 10 | ldt_alpha1 vs ar1_predictive | 30 | 0.1068 | 0.1507 | 1.863e-09 | 9.313e-09 | ldt_alpha1 |
| FDA4 | 10 | ldt_alpha1 vs gmm_inspired | 30 | 0.1068 | 0.1595 | 1.863e-09 | 9.313e-09 | ldt_alpha1 |
| FDA4 | 10 | ldt_alpha1 vs pps | 30 | 0.1068 | 0.1167 | 1.863e-09 | 9.313e-09 | ldt_alpha1 |
| FDA4 | 25 | ldt_alpha1 vs nsga2_static | 30 | 0.0858 | 0.1020 | 1.863e-09 | 9.313e-09 | ldt_alpha1 |
| FDA4 | 25 | ldt_alpha1 vs dnsga2_inspired | 30 | 0.0858 | 0.1027 | 1.863e-09 | 9.313e-09 | ldt_alpha1 |
| FDA4 | 25 | ldt_alpha1 vs ar1_predictive | 30 | 0.0858 | 0.0987 | 1.863e-09 | 9.313e-09 | ldt_alpha1 |
| FDA4 | 25 | ldt_alpha1 vs gmm_inspired | 30 | 0.0858 | 0.1011 | 1.863e-09 | 9.313e-09 | ldt_alpha1 |
| FDA4 | 25 | ldt_alpha1 vs pps | 30 | 0.0858 | 0.0881 | 5.718e-07 | 5.718e-07 | ldt_alpha1 |
| FDA4 | 50 | ldt_alpha1 vs nsga2_static | 30 | 0.0737 | 0.0764 | 9.313e-09 | 4.657e-08 | ldt_alpha1 |
| FDA4 | 50 | ldt_alpha1 vs gmm_inspired | 30 | 0.0737 | 0.0757 | 4.657e-08 | 1.863e-07 | ldt_alpha1 |
| FDA4 | 50 | ldt_alpha1 vs dnsga2_inspired | 30 | 0.0737 | 0.0759 | 1.024e-07 | 3.073e-07 | ldt_alpha1 |
| FDA4 | 50 | ldt_alpha1 vs ar1_predictive | 30 | 0.0737 | 0.0746 | 3.643e-02 | 7.287e-02 | ldt_alpha1 |
| FDA4 | 50 | ldt_alpha1 vs pps | 30 | 0.0737 | 0.0736 | 8.078e-01 | 8.078e-01 | pps |

### Friedman + Nemenyi [v25]

| problem | tau | k | N | chi2 | p | CD | avg ranks |
|---|---|---|---|---|---|---|---|
| DMOP1A | 25 | 6 | 30 | 47.66 | 4.17e-09 | 1.3767 | nsga2_static=3.73, dnsga2_inspired=4.40, ar1_predictive=2.83, gmm_inspired=5.00, pps=2.70, ldt_alpha1=2.33 |
| DMOP2A | 25 | 6 | 30 | 133.24 | 4.89e-27 | 1.3767 | nsga2_static=5.00, dnsga2_inspired=5.57, ar1_predictive=2.97, gmm_inspired=4.43, pps=1.53, ldt_alpha1=1.50 |
| DMOP3A | 25 | 6 | 30 | 76.08 | 5.55e-15 | 1.3767 | nsga2_static=4.30, dnsga2_inspired=5.07, ar1_predictive=2.97, gmm_inspired=4.53, pps=1.67, ldt_alpha1=2.47 |
| FDA1 | 10 | 6 | 30 | 126.04 | 1.65e-25 | 1.3767 | nsga2_static=4.63, dnsga2_inspired=5.60, ar1_predictive=3.30, gmm_inspired=4.47, pps=1.63, ldt_alpha1=1.37 |
| FDA1 | 25 | 6 | 30 | 134.32 | 2.88e-27 | 1.3767 | nsga2_static=5.30, dnsga2_inspired=5.37, ar1_predictive=3.00, gmm_inspired=4.33, pps=1.47, ldt_alpha1=1.53 |
| FDA1 | 50 | 6 | 30 | 122.21 | 1.07e-24 | 1.3767 | nsga2_static=5.10, dnsga2_inspired=5.43, ar1_predictive=2.40, gmm_inspired=4.47, pps=1.67, ldt_alpha1=1.93 |
| FDA3 | 10 | 6 | 30 | 130.36 | 1.99e-26 | 1.3767 | nsga2_static=4.37, dnsga2_inspired=5.60, ar1_predictive=3.20, gmm_inspired=4.83, pps=1.20, ldt_alpha1=1.80 |
| FDA3 | 25 | 6 | 30 | 135.58 | 1.55e-27 | 1.3767 | nsga2_static=5.03, dnsga2_inspired=5.60, ar1_predictive=3.00, gmm_inspired=4.37, pps=1.33, ldt_alpha1=1.67 |
| FDA3 | 50 | 6 | 30 | 118.30 | 7.17e-24 | 1.3767 | nsga2_static=5.07, dnsga2_inspired=4.97, ar1_predictive=2.37, gmm_inspired=4.97, pps=1.60, ldt_alpha1=2.03 |
| FDA4 | 10 | 6 | 30 | 139.54 | 2.24e-28 | 1.3767 | nsga2_static=6.00, dnsga2_inspired=4.50, ar1_predictive=3.30, gmm_inspired=4.20, pps=2.00, ldt_alpha1=1.00 |
| FDA4 | 25 | 6 | 30 | 110.55 | 3.13e-22 | 1.3767 | nsga2_static=4.63, dnsga2_inspired=5.07, ar1_predictive=4.00, gmm_inspired=4.30, pps=1.87, ldt_alpha1=1.13 |
| FDA4 | 50 | 6 | 30 | 77.43 | 2.89e-15 | 1.3767 | nsga2_static=4.93, dnsga2_inspired=4.60, ar1_predictive=3.10, gmm_inspired=4.43, pps=1.83, ldt_alpha1=2.10 |

overall mean rank (all cells averaged): pps=1.708, ldt_alpha1=1.739, ar1_predictive=3.036, gmm_inspired=4.528, nsga2_static=4.842, dnsga2_inspired=5.147

## paired alpha=1.0 (v2.5) vs alpha=0.5 (v2.4)

| problem | tau | n | MIGD a=1.0 | MIGD a=0.5 | improvement | p |
|---|---|---|---|---|---|---|
| FDA1 | 10 | 30 | 0.0120 | 0.0156 | +23.2% | 1.863e-09 |
| FDA1 | 25 | 30 | 0.0062 | 0.0074 | +16.1% | 1.863e-09 |
| FDA1 | 50 | 30 | 0.0049 | 0.0052 | +4.8% | 1.863e-09 |
| FDA3 | 10 | 30 | 0.0136 | 0.0177 | +22.8% | 1.863e-09 |
| FDA3 | 25 | 30 | 0.0066 | 0.0079 | +16.4% | 1.863e-09 |
| FDA3 | 50 | 30 | 0.0052 | 0.0054 | +4.3% | 3.725e-09 |
| FDA4 | 10 | 30 | 0.1068 | 0.1236 | +13.6% | 1.863e-09 |
| FDA4 | 25 | 30 | 0.0858 | 0.0904 | +5.0% | 1.863e-09 |
| FDA4 | 50 | 30 | 0.0737 | 0.0745 | +1.1% | 1.853e-02 |

## v2.5 winner per (problem, tau)

| problem | tau | best algo | best MIGD | ldt_alpha1 MIGD | gap |
|---|---|---|---|---|---|
| DMOP1A | 25 | ldt_alpha1 | 0.0570 | 0.0570 | +0.0% |
| DMOP2A | 25 | ldt_alpha1 | 0.0265 | 0.0265 | +0.0% |
| DMOP3A | 25 | pps | 0.0479 | 0.0525 | +9.7% |
| FDA1 | 10 | ldt_alpha1 | 0.0120 | 0.0120 | +0.0% |
| FDA1 | 25 | pps | 0.0062 | 0.0062 | +0.3% |
| FDA1 | 50 | pps | 0.0049 | 0.0049 | +0.2% |
| FDA3 | 10 | pps | 0.0132 | 0.0136 | +3.4% |
| FDA3 | 25 | pps | 0.0066 | 0.0066 | +0.9% |
| FDA3 | 50 | pps | 0.0051 | 0.0052 | +0.8% |
| FDA4 | 10 | ldt_alpha1 | 0.1068 | 0.1068 | +0.0% |
| FDA4 | 25 | ldt_alpha1 | 0.0858 | 0.0858 | +0.0% |
| FDA4 | 50 | pps | 0.0736 | 0.0737 | +0.2% |

## recovery generations (generations until IGD <= previous environment's final IGD; tau_t if never)

| problem | tau | algorithm | mean | std |
|---|---|---|---|---|
| DMOP1 | 25 | pps | 11.76 | 8.37 |
| DMOP1 | 25 | ar1_predictive | 14.61 | 8.97 |
| DMOP1 | 25 | ldt_nsga2 (alpha=0.5) | 20.68 | 6.72 |
| DMOP1 | 25 | dnsga2_inspired | 21.60 | 6.56 |
| DMOP1 | 25 | gmm_inspired | 21.61 | 6.72 |
| DMOP1 | 25 | nsga2_static | 21.63 | 6.63 |
| DMOP1A | 25 | ar1_predictive | 22.34 | 7.50 |
| DMOP1A | 25 | nsga2_static | 22.35 | 7.41 |
| DMOP1A | 25 | pps | 22.39 | 7.40 |
| DMOP1A | 25 | dnsga2_inspired | 22.41 | 7.35 |
| DMOP1A | 25 | ldt_alpha1 (alpha=1.0) | 22.48 | 7.32 |
| DMOP1A | 25 | gmm_inspired | 22.57 | 7.07 |
| DMOP2 | 25 | pps | 10.56 | 8.11 |
| DMOP2 | 25 | ar1_predictive | 14.43 | 8.75 |
| DMOP2 | 25 | ldt_nsga2 (alpha=0.5) | 20.71 | 6.69 |
| DMOP2 | 25 | gmm_inspired | 21.51 | 6.70 |
| DMOP2 | 25 | nsga2_static | 21.52 | 6.61 |
| DMOP2 | 25 | dnsga2_inspired | 21.58 | 6.60 |
| DMOP2A | 25 | pps | 17.92 | 6.98 |
| DMOP2A | 25 | ar1_predictive | 17.98 | 6.88 |
| DMOP2A | 25 | ldt_alpha1 (alpha=1.0) | 18.79 | 6.70 |
| DMOP2A | 25 | gmm_inspired | 21.81 | 5.05 |
| DMOP2A | 25 | dnsga2_inspired | 22.04 | 4.92 |
| DMOP2A | 25 | nsga2_static | 22.06 | 5.12 |
| DMOP3 | 25 | gmm_inspired | 0.00 | 0.00 |
| DMOP3 | 25 | ldt_nsga2 (alpha=0.5) | 0.00 | 0.00 |
| DMOP3 | 25 | nsga2_static | 0.00 | 0.00 |
| DMOP3 | 25 | dnsga2_inspired | 0.24 | 1.61 |
| DMOP3 | 25 | ar1_predictive | 1.75 | 5.96 |
| DMOP3 | 25 | pps | 2.26 | 6.80 |
| DMOP3A | 25 | pps | 20.47 | 5.81 |
| DMOP3A | 25 | ldt_alpha1 (alpha=1.0) | 20.53 | 5.85 |
| DMOP3A | 25 | ar1_predictive | 21.06 | 4.98 |
| DMOP3A | 25 | dnsga2_inspired | 21.12 | 4.37 |
| DMOP3A | 25 | gmm_inspired | 21.25 | 4.53 |
| DMOP3A | 25 | nsga2_static | 21.41 | 4.29 |
| FDA1 | 10 | pps | 6.07 | 3.21 |
| FDA1 | 10 | ar1_predictive | 6.64 | 3.31 |
| FDA1 | 10 | nsga2_static | 7.39 | 3.12 |
| FDA1 | 10 | ldt_alpha1 (alpha=1.0) | 7.59 | 2.42 |
| FDA1 | 10 | dnsga2_inspired | 7.70 | 2.69 |
| FDA1 | 10 | ldt_nsga2 (alpha=0.5) | 7.91 | 2.33 |
| FDA1 | 10 | gmm_inspired | 8.09 | 2.24 |
| FDA1 | 25 | pps | 10.29 | 6.55 |
| FDA1 | 25 | ar1_predictive | 16.02 | 7.43 |
| FDA1 | 25 | ldt_alpha1 (alpha=1.0) | 16.57 | 6.09 |
| FDA1 | 25 | ldt_nsga2 (alpha=0.5) | 17.25 | 7.66 |
| FDA1 | 25 | LDT_AblationC | 21.30 | 4.83 |
| FDA1 | 25 | nsga2_static | 21.86 | 4.95 |
| FDA1 | 25 | dnsga2_inspired | 21.98 | 4.93 |
| FDA1 | 25 | gmm_inspired | 22.01 | 4.73 |
| FDA1 | 25 | LDT_AblationA | 22.11 | 4.74 |
| FDA1 | 25 | LDT_AblationB | 22.43 | 3.89 |
| FDA1 | 50 | pps | 17.63 | 15.58 |
| FDA1 | 50 | ldt_alpha1 (alpha=1.0) | 24.99 | 12.55 |
| FDA1 | 50 | ar1_predictive | 28.08 | 15.61 |
| FDA1 | 50 | ldt_nsga2 (alpha=0.5) | 43.12 | 6.02 |
| FDA1 | 50 | gmm_inspired | 44.62 | 5.91 |
| FDA1 | 50 | dnsga2_inspired | 45.25 | 4.85 |
| FDA1 | 50 | nsga2_static | 45.37 | 4.60 |
| FDA3 | 10 | pps | 6.82 | 2.76 |
| FDA3 | 10 | ar1_predictive | 6.90 | 3.12 |
| FDA3 | 10 | nsga2_static | 7.40 | 3.10 |
| FDA3 | 10 | dnsga2_inspired | 7.81 | 2.55 |
| FDA3 | 10 | gmm_inspired | 8.21 | 2.12 |
| FDA3 | 10 | ldt_alpha1 (alpha=1.0) | 8.29 | 2.15 |
| FDA3 | 10 | ldt_nsga2 (alpha=0.5) | 8.45 | 2.07 |
| FDA3 | 25 | pps | 13.97 | 6.79 |
| FDA3 | 25 | ar1_predictive | 17.03 | 6.60 |
| FDA3 | 25 | ldt_nsga2 (alpha=0.5) | 18.25 | 6.93 |
| FDA3 | 25 | ldt_alpha1 (alpha=1.0) | 18.32 | 5.68 |
| FDA3 | 25 | LDT_AblationC | 21.70 | 4.71 |
| FDA3 | 25 | gmm_inspired | 22.11 | 4.74 |
| FDA3 | 25 | nsga2_static | 22.13 | 4.79 |
| FDA3 | 25 | dnsga2_inspired | 22.13 | 4.80 |
| FDA3 | 25 | LDT_AblationA | 22.16 | 4.75 |
| FDA3 | 25 | LDT_AblationB | 22.53 | 3.82 |
| FDA3 | 50 | pps | 22.63 | 15.42 |
| FDA3 | 50 | ldt_alpha1 (alpha=1.0) | 28.30 | 13.30 |
| FDA3 | 50 | ar1_predictive | 32.69 | 14.82 |
| FDA3 | 50 | ldt_nsga2 (alpha=0.5) | 44.60 | 5.56 |
| FDA3 | 50 | gmm_inspired | 45.68 | 5.52 |
| FDA3 | 50 | nsga2_static | 45.90 | 4.66 |
| FDA3 | 50 | dnsga2_inspired | 46.61 | 4.29 |
| FDA4 | 10 | pps | 5.06 | 3.74 |
| FDA4 | 10 | ar1_predictive | 5.71 | 3.87 |
| FDA4 | 10 | ldt_alpha1 (alpha=1.0) | 6.06 | 3.51 |
| FDA4 | 10 | nsga2_static | 6.40 | 4.03 |
| FDA4 | 10 | dnsga2_inspired | 6.71 | 3.59 |
| FDA4 | 10 | gmm_inspired | 6.72 | 3.56 |
| FDA4 | 10 | ldt_nsga2 (alpha=0.5) | 6.74 | 3.41 |
| FDA4 | 25 | pps | 9.65 | 9.33 |
| FDA4 | 25 | ldt_alpha1 (alpha=1.0) | 13.63 | 8.91 |
| FDA4 | 25 | ldt_nsga2 (alpha=0.5) | 14.12 | 9.14 |
| FDA4 | 25 | ar1_predictive | 16.05 | 9.27 |
| FDA4 | 25 | LDT_AblationC | 17.71 | 7.64 |
| FDA4 | 25 | LDT_AblationA | 20.55 | 7.20 |
| FDA4 | 25 | nsga2_static | 20.64 | 7.16 |
| FDA4 | 25 | dnsga2_inspired | 20.67 | 7.18 |
| FDA4 | 25 | gmm_inspired | 20.69 | 7.15 |
| FDA4 | 25 | LDT_AblationB | 20.80 | 6.65 |
| FDA4 | 50 | pps | 16.99 | 18.62 |
| FDA4 | 50 | ldt_alpha1 (alpha=1.0) | 18.73 | 17.24 |
| FDA4 | 50 | ldt_nsga2 (alpha=0.5) | 26.39 | 17.32 |
| FDA4 | 50 | ar1_predictive | 27.09 | 19.88 |
| FDA4 | 50 | gmm_inspired | 34.19 | 18.32 |
| FDA4 | 50 | nsga2_static | 34.53 | 18.26 |
| FDA4 | 50 | dnsga2_inspired | 35.18 | 18.08 |
