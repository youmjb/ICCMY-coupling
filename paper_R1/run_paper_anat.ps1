# Paper R1, anatomy check: does g_c* / G_load stay ~0.6-0.9 for a short, wide column closer to mouse anatomy
# (L 8 layers, m 3 SMC per ICC-MY per layer)? g_SS 30 nS (G_load ~14 nS, tanh-corrected) and 300 nS (~21 nS); v4 base otherwise.
# Output: out\paper_anat\gSS<g>_gc<gc>
$d = Split-Path -Parent (Split-Path -Parent $MyInvocation.MyCommand.Path)
Set-Location $d
$env:OMP_NUM_THREADS = 1
New-Item -ItemType Directory -Force out\paper_anat | Out-Null
$jobs = @()
foreach ($gc in 3,6,10,15,25,40) { $jobs += ,@(30, $gc) }
foreach ($gc in 5,10,15,20,30,50) { $jobs += ,@(300, $gc) }
$running = @()
foreach ($j in $jobs) {
    $pre = "out\paper_anat\gSS$($j[0])_gc$($j[1])"
    if (Test-Path "$pre.log") { continue }
    $argl = @("active","$($j[1])","0","0","40000","0.1","$pre.csv","m=3","L=8","gSS=$($j[0])","fI=0","fP=0.1","gSP=3","pdModel=1","smcModel=1","gK0run=0.4","iccScale=0.72")
    $running += Start-Process -PassThru -WindowStyle Hidden -FilePath "$d\bin\gisim.exe" -WorkingDirectory $d -RedirectStandardOutput "$pre.log" -ArgumentList $argl
}
$running | ForEach-Object { $_.WaitForExit() }
"all done"
