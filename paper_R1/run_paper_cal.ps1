# Paper R1: SK-blocked (apamin = skKd 1000) runs of the calibrated column at g_c 0, 1, 1.5, 2 nS -> out\paper_cal\apa_gc<gc>
# (made by hand on 2026-10-14 for Fig. 1C and the SK decomposition; same command form as run_paper_rev1c.ps1). Skips existing logs.
$d = Split-Path -Parent (Split-Path -Parent $MyInvocation.MyCommand.Path)
Set-Location $d
$env:OMP_NUM_THREADS = 1
New-Item -ItemType Directory -Force "out\paper_cal" | Out-Null
$base = @("m=0.3","skKd=1000","fI=0","fP=0.1","gSP=3","pdModel=1","smcModel=1","gK0run=0.4","iccScale=0.72")
$running = @()
foreach ($gc in "0","1","1.5","2") {
    $pre = "out\paper_cal\apa_gc$gc"
    if (Test-Path "$pre.log") { continue }
    $argl = @("active",$gc,"0","0","40000","0.1","$pre.csv") + $base
    $running += Start-Process -PassThru -WindowStyle Hidden -FilePath "$d\bin\gisim.exe" -WorkingDirectory $d -RedirectStandardOutput "$pre.log" -ArgumentList $argl
}
$running | ForEach-Object { $_.WaitForExit() }
