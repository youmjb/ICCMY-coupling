# Paper R1, FLC step 3 (robustness to the FLC model): v4 tissue with the modified PDGFa FLC replaced by the minimal mean-field
# PDGFRa cell (pdModel=0: Cm 25 pF, background 0.3 nS at -20 mV, SK 58 nS with Ca sigmoid Kd 0.4 uM, Hill 4, tau 10 ms; purinergic
# drive raises Ca by 1.5 uM x u with tau 150 ms; no STHs). g_c sweep 0/0.5/1/1.5/2/3 (40 s, no input) and the 0.5-s purinergic pulse
# at 8 phases + control at g_c 1.5 (30 s), as run_v4ijp.ps1. Output: out\paper_flcmf\
$d = Split-Path -Parent (Split-Path -Parent $MyInvocation.MyCommand.Path)
Set-Location $d
$env:OMP_NUM_THREADS = 1
New-Item -ItemType Directory -Force out\paper_flcmf | Out-Null
$kv = @("m=0.3","fI=0","fP=0.1","gSP=3","pdModel=0","smcModel=1","gK0run=0.4","iccScale=0.72")
$jobs = New-Object System.Collections.ArrayList
foreach ($gc in 0,0.5,1,1.5,2,3) { [void]$jobs.Add(@("out\paper_flcmf\gc$($gc)_u0", $gc, 0, 40000, @())) }
[void]$jobs.Add(@("out\paper_flcmf\ijp_ctrl", 1.5, 0, 30000, @()))
foreach ($t0 in 20000,20450,20900,21350,21800,22250,22700,23150) { [void]$jobs.Add(@("out\paper_flcmf\ijp_p$t0", 1.5, 1, 30000, @("uOn=$t0","uOff=$($t0 + 500)"))) }
$running = @()
foreach ($j in $jobs) {
    $pre = $j[0]
    if (Test-Path "$pre.log") { continue }
    $argl = @("active","$($j[1])","$($j[2])","0","$($j[3])","0.1","$pre.csv") + $kv + $j[4]
    $running += Start-Process -PassThru -WindowStyle Hidden -FilePath "$d\bin\gisim.exe" -WorkingDirectory $d -RedirectStandardOutput "$pre.log" -ArgumentList $argl
}
$running | ForEach-Object { $_.WaitForExit() }
"all done"
