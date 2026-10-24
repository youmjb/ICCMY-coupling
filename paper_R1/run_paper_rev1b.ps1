# Paper R1, revision 1b: time-step convergence of the g_c estimate at the calibrated column (dt 0.1 ms showed a monotonic dt trend
# in the muscle amplitude, 23.2 / 24.1 / 24.9 mV at dt 0.1 / 0.05 / 0.025 ms, g_c 1 nS). Adds g_c 0.5 and 2 nS at dt 0.025 ms
# (so the ratio-0.429 crossing is bracketed at that dt) and g_c 1 nS at dt 0.0125 ms. Output: out\paper_rev1\dt<dt>_gc<gc>
$d = Split-Path -Parent (Split-Path -Parent $MyInvocation.MyCommand.Path)
Set-Location $d
$env:OMP_NUM_THREADS = 1
$base = @("m=0.3","iccScale=0.72","fI=0","fP=0.1","gSP=3","pdModel=1","smcModel=1","gK0run=0.4")
$jobs = @(@("0.025","0.5"), @("0.025","0.75"), @("0.025","2"), @("0.0125","1"), @("0.0125","0.75"))
$running = @()
foreach ($j in $jobs) {
    $pre = "out\paper_rev1\dt$($j[0])_gc$($j[1])"
    if (Test-Path "$pre.log") { continue }
    $argl = @("active",$j[1],"0","0","40000",$j[0],"$pre.csv") + $base
    $running += Start-Process -PassThru -WindowStyle Hidden -FilePath "$d\bin\gisim.exe" -WorkingDirectory $d -RedirectStandardOutput "$pre.log" -ArgumentList $argl
}
"started $($running.Count) runs"
$running | ForEach-Object { $_.WaitForExit() }
"all done"
