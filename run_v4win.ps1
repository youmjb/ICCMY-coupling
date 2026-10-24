# Results 1 (coupling window) and 2 (load) with the v4 candidate (ICC-MY 18 pF, iccScale 0.72; SMS SMC, gK0 0.4; modified PDGFa FLC, gSP 3)
#   window: m 0.3, gSS 30, g_c sweep x (no input, sustained purinergic)
#   load:   m x gSS grid, g_c sweep, no input
$d = Split-Path -Parent $MyInvocation.MyCommand.Path
Set-Location $d
$env:OMP_NUM_THREADS = 1
New-Item -ItemType Directory -Force out\v4win | Out-Null
$jobs = @()
$jobs += ,@(0.3, 30, 0, 0)
foreach ($gc in 0.1,0.2,0.5,1,1.5,2,3,5,10,20,50,100,200) { foreach ($u in 0,1) { $jobs += ,@(0.3, 30, $gc, $u) } }
foreach ($m in 0.1,0.3,1,3,10) { foreach ($gSS in 10,30,100) { foreach ($gc in 0.5,1,2,3,5,10,20,50) {
    if (-not ($m -eq 0.3 -and $gSS -eq 30)) { $jobs += ,@($m, $gSS, $gc, 0) }
} } }
$jobs = $jobs | Sort-Object { -$_[2] }
$running = @()
foreach ($j in $jobs) {
    $m=$j[0]; $gSS=$j[1]; $gc=$j[2]; $u=$j[3]
    $tag = "m$($m)_gSS$($gSS)_gc$($gc)_u$u"
    if (Test-Path "out\v4win\$tag.log") { continue }    # already started or finished (delete the log to rerun)
    while (($running | Where-Object { -not $_.HasExited }).Count -ge 28) { Start-Sleep 5 }
    $argl = @("active","$gc","$u","0","40000","0.1","out\v4win\$tag.csv","m=$m","gSS=$gSS","fI=0","fP=0.1","gSP=3","pdModel=1","smcModel=1","gK0run=0.4","iccScale=0.72")
    $running += Start-Process -PassThru -WindowStyle Hidden -FilePath "$d\bin\gisim.exe" -WorkingDirectory $d -RedirectStandardOutput "out\v4win\$tag.log" -ArgumentList $argl
}
$running | ForEach-Object { $_.WaitForExit() }
"all done"
