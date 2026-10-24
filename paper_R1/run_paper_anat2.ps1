# Same CM : ICC-MY membrane ratio as the model column (m*L = 9 SMC per ICC-MY) but mouse-like thickness: L 7, m 1.3, g_SS 30.
$d = Split-Path -Parent (Split-Path -Parent $MyInvocation.MyCommand.Path)
Set-Location $d
$env:OMP_NUM_THREADS = 1
$running = @()
foreach ($gc in 1,2,3,4,6,8) {
    $pre = "out\paper_anat\L7m1.3_gc$gc"
    if (Test-Path "$pre.log") { continue }
    $argl = @("active","$gc","0","0","40000","0.1","$pre.csv","m=1.3","L=7","gSS=30","fI=0","fP=0.1","gSP=3","pdModel=1","smcModel=1","gK0run=0.4","iccScale=0.72")
    $running += Start-Process -PassThru -WindowStyle Hidden -FilePath "$d\bin\gisim.exe" -WorkingDirectory $d -RedirectStandardOutput "$pre.log" -ArgumentList $argl
}
$running | ForEach-Object { $_.WaitForExit() }
