# Paper R1, revision 2 (event-aligned amplitudes): extra g_c values where the mouse ratio 0.604 is not reached by the saved sweep.
#   one-at-a-time sensitivity (Fig. 6) where g_c* > 3 nS: m 0.45 g_c 4, ICC-MY +gL 1.8 / 3.6 nS g_c 4, 5
#   g_SS 60 column: g_c 4, 5 -> out\paper_sens\gSS60\gc<gc>_u0 ;  L 7, m 1.3 column: g_c 10, 12 -> out\paper_anat\L7m1.3_gc<gc>
$d = Split-Path -Parent (Split-Path -Parent $MyInvocation.MyCommand.Path)
Set-Location $d
$env:OMP_NUM_THREADS = 1
$base = @("fI=0","fP=0.1","gSP=3","pdModel=1","smcModel=1","gK0run=0.4","iccScale=0.72")
$jobs = @()
foreach ($gc in "4","5") { $jobs += ,@("out\paper_sens\gSS60\gc$($gc)_u0", $gc, (@("m=0.3","gSS=60","L=30","iccGL=0","seed=1") + $base)) }
$jobs += ,@("out\paper_sens\m0.45\gc4_u0", "4", (@("m=0.45","gSS=30","L=30","iccGL=0","seed=1") + $base))
foreach ($gc in "4","5") { foreach ($gl in "1.8","3.6") { $jobs += ,@("out\paper_sens\iccGL$gl\gc$($gc)_u0", $gc, (@("m=0.3","gSS=30","L=30","iccGL=$gl","seed=1") + $base)) } }
foreach ($gc in "10","12") { $jobs += ,@("out\paper_anat\L7m1.3_gc$gc", $gc, (@("m=1.3","L=7","gSS=30") + $base)) }
$running = @()
foreach ($j in $jobs) {
    if (Test-Path "$($j[0]).log") { continue }
    $argl = @("active",$j[1],"0","0","40000","0.1","$($j[0]).csv") + $j[2]
    $running += Start-Process -PassThru -WindowStyle Hidden -FilePath "$d\bin\gisim.exe" -WorkingDirectory $d -RedirectStandardOutput "$($j[0]).log" -ArgumentList $argl
}
"started $($running.Count) runs"
$running | ForEach-Object { $_.WaitForExit() }
"all done"
