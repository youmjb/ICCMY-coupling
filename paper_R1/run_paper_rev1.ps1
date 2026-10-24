# Paper R1, revision 1 (reviewer, 2026-10-15). No model-code change; v4 base keys.
# (1) time-step convergence at the calibrated column: g_c 1 and 1.5 nS with dt 0.05 and 0.025 ms (base runs used 0.1) -> out\paper_rev1\dt<dt>_gc<gc>
# (2) generality of the load limit: mouse-like column (L 7, m 1.3, G_load 5.9 nS) with a different ICC-MY own conductance
#     (iccGL -0.6 and +1.8 nS) and a larger ICC-MY (iccScale 1.0, 25 pF) -> out\paper_rev1\L7_<variant>_gc<gc>
$d = Split-Path -Parent (Split-Path -Parent $MyInvocation.MyCommand.Path)
Set-Location $d
$env:OMP_NUM_THREADS = 1
New-Item -ItemType Directory -Force out\paper_rev1 | Out-Null
$base = @("fI=0","fP=0.1","gSP=3","pdModel=1","smcModel=1","gK0run=0.4")
$jobs = @()
foreach ($dt in "0.05","0.025") { foreach ($gc in "1","1.5") {
    $jobs += ,@("out\paper_rev1\dt$($dt)_gc$gc", $gc, $dt, (@("m=0.3","iccScale=0.72") + $base)) } }
foreach ($v in @(@("gLm0.6", "iccGL=-0.6", "iccScale=0.72"), @("gL1.8", "iccGL=1.8", "iccScale=0.72"), @("s1.0", "iccGL=0", "iccScale=1.0"))) {
    foreach ($gc in "1","2","3","4","6","8","12") {
        $jobs += ,@("out\paper_rev1\L7_$($v[0])_gc$gc", $gc, "0.1", (@("m=1.3","L=7",$v[1],$v[2]) + $base)) } }
$running = @()
foreach ($j in $jobs) {
    if (Test-Path "$($j[0]).log") { continue }
    $argl = @("active",$j[1],"0","0","40000",$j[2],"$($j[0]).csv") + $j[3]
    $running += Start-Process -PassThru -WindowStyle Hidden -FilePath "$d\bin\gisim.exe" -WorkingDirectory $d -RedirectStandardOutput "$($j[0]).log" -ArgumentList $argl
}
"started $($running.Count) runs"
$running | ForEach-Object { $_.WaitForExit() }
"all done"
