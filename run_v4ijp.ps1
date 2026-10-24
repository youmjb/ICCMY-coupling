# Result 3 with the v4 candidate (ICC-MY 18 pF, iccScale 0.72, SMC density m 0.3): 0.5 s purinergic pulse at 8 phases (+ control), g_c 1 / 1.5 / 2 nS
$d = Split-Path -Parent $MyInvocation.MyCommand.Path
Set-Location $d
$env:OMP_NUM_THREADS = 1
New-Item -ItemType Directory -Force out\v4ijp | Out-Null
$jobs = @()
foreach ($gc in 1,1.5,2) {
    $jobs += ,@($gc, 'ctrl', 0, -1, -1)
    foreach ($t0 in 20000,20450,20900,21350,21800,22250,22700,23150) { $jobs += ,@($gc, "p$t0", 1, $t0, ($t0 + 500)) }
}
$running = @()
foreach ($j in $jobs) {
    while (($running | Where-Object { -not $_.HasExited }).Count -ge 23) { Start-Sleep 5 }
    $gc=$j[0]; $tag=$j[1]; $u=$j[2]; $on=$j[3]; $off=$j[4]
    $argl = @("active","$gc","$u","0","30000","0.1","out\v4ijp\gc$($gc)_$tag.csv","m=0.3","fI=0","fP=0.1","gSP=3","pdModel=1","smcModel=1","gK0run=0.4","iccScale=0.72")
    if ($u -gt 0) { $argl += @("uOn=$on","uOff=$off") }
    $running += Start-Process -PassThru -WindowStyle Hidden -FilePath "$d\bin\gisim.exe" -WorkingDirectory $d -RedirectStandardOutput "out\v4ijp\gc$($gc)_$tag.log" -ArgumentList $argl
}
$running | ForEach-Object { $_.WaitForExit() }
"all done"
