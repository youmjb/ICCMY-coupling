# Paper R1, revision 1c: runs for the mouse amplitude-ratio target (Kito & Suzuki 2003, ratio 0.604).
# (1) SK blocked (apamin = skKd 1000) at g_c 2.5 and 3 nS, calibrated column, to bracket ratio 0.604 for the passive/SK/SMC
#     decomposition; g_c 1.5 repeated (check that the command reproduces out/paper_cal/apa_gc1.5.log) -> out\paper_rev1\apa_gc<gc>
# (2) m 0.45 column at g_c 5 nS (ratio 0.604 not reached by g_c 3) -> out\paper_sens\m0.45\gc5_u0
$d = Split-Path -Parent (Split-Path -Parent $MyInvocation.MyCommand.Path)
Set-Location $d
$env:OMP_NUM_THREADS = 1
$base = @("fI=0","fP=0.1","gSP=3","pdModel=1","smcModel=1","gK0run=0.4","iccScale=0.72")
$jobs = @()
foreach ($gc in "1.5","2.5","3") { $jobs += ,@("out\paper_rev1\apa_gc$gc", $gc, (@("m=0.3","skKd=1000") + $base)) }
$jobs += ,@("out\paper_sens\m0.45\gc5_u0", "5", (@("m=0.45") + $base))
$running = @()
foreach ($j in $jobs) {
    if (Test-Path "$($j[0]).log") { continue }
    $argl = @("active",$j[1],"0","0","40000","0.1","$($j[0]).csv") + $j[2]
    $running += Start-Process -PassThru -WindowStyle Hidden -FilePath "$d\bin\gisim.exe" -WorkingDirectory $d -RedirectStandardOutput "$($j[0]).log" -ArgumentList $argl
}
"started $($running.Count) runs"
$running | ForEach-Object { $_.WaitForExit() }
"all done"
