$py = "D:\python\python.exe"
$proj = "D:\23-兴湘计科\轻量级分布迁移动态多目标优化器\LDT-NSGA-II"
Set-Location $proj
New-Item -ItemType Directory -Force -Path results\raw_v25, results\_logs_v25 | Out-Null
$jobs = @()
foreach ($p in @("FDA1","FDA3","FDA4")) {
  foreach ($tau in @(10,25,50)) {
    $log = "results\_logs_v25\cmp_${p}_tau${tau}.log"
    $jobs += Start-Job -Name "cmp_${p}_${tau}" -ScriptBlock {
      param($py,$proj,$p,$tau,$log)
      Set-Location $proj
      & $py -m experiments.run_v25_alpha --problems $p --taus $tau --seeds 0-29 --out results/raw_v25 --phase comparison *> $log
    } -ArgumentList $py,$proj,$p,$tau,$log
  }
}
Write-Output "LAUNCHED $($jobs.Count) jobs: $($jobs.Name -join ', ')"
$jobs | Wait-Job | Out-Null
$jobs | ForEach-Object { "=== $($_.Name) state=$($_.State) ==="; Receive-Job $_ | Select-Object -Last 3 }
$jobs | Remove-Job
Write-Output "COMPARISON_PHASE_DONE"
