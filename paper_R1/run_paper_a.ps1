# Paper R1, priority 2 (does the ICC-MY RMP defect matter?) and direct v4 R_in. v4 base keys, no model change.
#   out\paper_rmp : sustained inward current into ICC-MY (iccInj < 0, whole run) to depolarize its rest by ~0-10 mV;
#                   g_c 0/1/1.5/2/5 x iccInj 0/-10/-20/-30/-40 pA, 40 s, no neural input (slow-wave amplitudes from the log)
#   out\paper_rijp: 0.5 s purinergic pulse at 8 phases (+ control) as run_v4ijp.ps1, but with iccInj -30 pA throughout; g_c 1/1.5/2
#   out\paper_rin : +10 pA (outward) 300 ms into ICC-MY at 8 phases (+ control), m 0.3 g_c 0/1/1.5/2/5/20 and m 1 g_c 1/2/5, 25 s
$d = Split-Path -Parent (Split-Path -Parent $MyInvocation.MyCommand.Path)
Set-Location $d
$env:OMP_NUM_THREADS = 1
foreach ($f in 'paper_rmp','paper_rijp','paper_rin') { New-Item -ItemType Directory -Force "out\$f" | Out-Null }
$base = @("fI=0","fP=0.1","gSP=3","pdModel=1","smcModel=1","gK0run=0.4","iccScale=0.72")
$phases = 20000,20450,20900,21350,21800,22250,22700,23150
$jobs = @()   # each job: @(outPrefix, gc, u, T, extraArgs)
foreach ($gc in 0,1,1.5,2,5) { foreach ($i in 0,-10,-20,-30,-40) {
    $jobs += ,@("out\paper_rmp\gc$($gc)_inj$i", $gc, 0, 40000, @("m=0.3","iccInj=$i","injOn=0","injOff=40000"))
} }
foreach ($gc in 1,1.5,2) {
    $jobs += ,@("out\paper_rijp\gc$($gc)_ctrl", $gc, 0, 30000, @("m=0.3","iccInj=-30","injOn=0","injOff=30000"))
    foreach ($t0 in $phases) { $jobs += ,@("out\paper_rijp\gc$($gc)_p$t0", $gc, 1, 30000, @("m=0.3","iccInj=-30","injOn=0","injOff=30000","uOn=$t0","uOff=$($t0 + 500)")) }
}
foreach ($c in @(@(0.3,0),@(0.3,1),@(0.3,1.5),@(0.3,2),@(0.3,5),@(0.3,20),@(1,1),@(1,2),@(1,5))) {
    $m = $c[0]; $gc = $c[1]
    $jobs += ,@("out\paper_rin\m$($m)_gc$($gc)_ctrl", $gc, 0, 25000, @("m=$m"))
    foreach ($t0 in $phases) { $jobs += ,@("out\paper_rin\m$($m)_gc$($gc)_i$t0", $gc, 0, 25000, @("m=$m","iccInj=10","injOn=$t0","injOff=$($t0 + 300)")) }
}
$running = @()
foreach ($j in $jobs) {
    $pre = $j[0]
    if (Test-Path "$pre.log") { continue }    # already started or finished (delete the log to rerun)
    while (($running | Where-Object { -not $_.HasExited }).Count -ge 28) { Start-Sleep 5 }
    $argl = @("active","$($j[1])","$($j[2])","0","$($j[3])","0.1","$pre.csv") + $base + $j[4]
    $running += Start-Process -PassThru -WindowStyle Hidden -FilePath "$d\bin\gisim.exe" -WorkingDirectory $d -RedirectStandardOutput "$pre.log" -ArgumentList $argl
}
$running | ForEach-Object { $_.WaitForExit() }
"all done"
