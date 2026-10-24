# Paper R1 step 3 (robustness to the SMC model): v4 tissue with SmoothM replaced by Poh et al. 2012 human jejunal SMC (smcModel=2,
# src/PohSMC.h, GRL1 sub-steps 0.01 ms). g_c sweep 0/0.5/1/1.5/2/3/5 nS, 40 s, no input; tonic K+ 0 (the model's own rest, -73.5 mV)
# and 0.4 nS (as SmoothM). Output: out\paper_poh\gK<g>_gc<gc>
$d = Split-Path -Parent (Split-Path -Parent $MyInvocation.MyCommand.Path)
Set-Location $d
$env:OMP_NUM_THREADS = 1
New-Item -ItemType Directory -Force out\paper_poh | Out-Null
$running = @()
foreach ($gk in 0,0.4) { foreach ($gc in 0,0.5,1,1.5,2,3,5) {
    $pre = "out\paper_poh\gK$($gk)_gc$gc"
    if (Test-Path "$pre.log") { continue }
    $argl = @("active","$gc","0","0","40000","0.1","$pre.csv","m=0.3","fI=0","fP=0.1","gSP=3","pdModel=1","smcModel=2","gK0run=$gk","iccScale=0.72")
    $running += Start-Process -PassThru -WindowStyle Hidden -FilePath "$d\bin\gisim.exe" -WorkingDirectory $d -RedirectStandardOutput "$pre.log" -ArgumentList $argl
} }
$running | ForEach-Object { $_.WaitForExit() }
"all done"
