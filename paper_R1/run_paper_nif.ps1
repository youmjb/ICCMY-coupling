# Paper R1: nifedipine condition (Kito et al. 2014 recorded in 3 uM nifedipine throughout). SmoothM L-type Ca2+ conductance scaled by
# env SMS_CAL (0.1 = 90 % block, 0 = full block; default 1). ICC-MY unchanged (Kito: nifedipine had little effect on slow waves).
# v4 base keys. Per level: g_c sweep 0/0.5/1/1.5/2/3/5 nS (40 s) and the 0.5 s purinergic pulse at 8 phases + control at g_c 1.5 (30 s).
# Output: out\paper_nif\cal<level>\
$d = Split-Path -Parent (Split-Path -Parent $MyInvocation.MyCommand.Path)
Set-Location $d
$env:OMP_NUM_THREADS = 1
$kv = @("m=0.3","fI=0","fP=0.1","gSP=3","pdModel=1","smcModel=1","gK0run=0.4","iccScale=0.72")
$jobs = @()
foreach ($cal in "0.1","0") {
    $dir = "out\paper_nif\cal$cal"; New-Item -ItemType Directory -Force $dir | Out-Null
    foreach ($gc in 0,0.5,1,1.5,2,3,5) { $jobs += ,@($cal, "$dir\gc$($gc)_u0", $gc, 0, 40000, $kv) }
    $jobs += ,@($cal, "$dir\ijp_ctrl", 1.5, 0, 30000, $kv)
    foreach ($t0 in 20000,20450,20900,21350,21800,22250,22700,23150) { $jobs += ,@($cal, "$dir\ijp_p$t0", 1.5, 1, 30000, ($kv + @("uOn=$t0","uOff=$($t0 + 500)"))) }
}
$running = @()
foreach ($j in $jobs) {
    $pre = $j[1]
    if (Test-Path "$pre.log") { continue }
    while (($running | Where-Object { -not $_.HasExited }).Count -ge 28) { Start-Sleep 5 }
    $env:SMS_CAL = $j[0]    # inherited by the child process
    $argl = @("active","$($j[2])","$($j[3])","0","$($j[4])","0.1","$pre.csv") + $j[5]
    $running += Start-Process -PassThru -WindowStyle Hidden -FilePath "$d\bin\gisim.exe" -WorkingDirectory $d -RedirectStandardOutput "$pre.log" -ArgumentList $argl
}
$running | ForEach-Object { $_.WaitForExit() }
"all done"
