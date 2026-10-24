# Paper R1: does the ICC-MY slow-wave duration matter? ICC-MY Vm prescribed as a trapezoid (keys wavePer/waveHW/waveRise/waveFall/
# waveRest/waveAmp, added 2026-10-14), v4 tissue (m 0.3, gSS 30, ...), g_c 0.5-5 nS, 40 s, no neural input. Output: out\paper_wave\
#   W1 model-like  (period 3.33 s, half-width 0.67 s, rise 60 ms, fall 400 ms, -70 -> -20 mV)
#   W2 period 3.33, HW 1.5 | W3 period 6, HW 0.67 | W4 period 6, HW 1.5 | W5 period 6, HW 2.4 (Kito half-width and rate)
#   W6 = W5 with rise 5 ms (Kito dV/dt ~10 V/s) | W7 = W6 with rest -65 mV and amplitude 57 mV (Kito ICC-MY)
# plus a regression run (wavePer 0) that must reproduce out\v4win\m0.3_gSS30_gc1.5_u0.csv.
$d = Split-Path -Parent (Split-Path -Parent $MyInvocation.MyCommand.Path)
Set-Location $d
$env:OMP_NUM_THREADS = 1
New-Item -ItemType Directory -Force out\paper_wave | Out-Null
$kv = @("m=0.3","fI=0","fP=0.1","gSP=3","pdModel=1","smcModel=1","gK0run=0.4","iccScale=0.72")
$W = @(
    @("W1", @("wavePer=3330","waveHW=670","waveRise=60")),  @("W2", @("wavePer=3330","waveHW=1500","waveRise=60")),
    @("W3", @("wavePer=6000","waveHW=670","waveRise=60")),  @("W4", @("wavePer=6000","waveHW=1500","waveRise=60")),
    @("W5", @("wavePer=6000","waveHW=2400","waveRise=60")), @("W6", @("wavePer=6000","waveHW=2400","waveRise=5")),
    @("W7", @("wavePer=6000","waveHW=2400","waveRise=5","waveRest=-65","waveAmp=57"))
)
$jobs = New-Object System.Collections.ArrayList
[void]$jobs.Add(@("out\paper_wave\regress_gc1.5", 1.5, @()))
foreach ($w in $W) { foreach ($gc in 0.5,1,1.5,2,3,5) { [void]$jobs.Add(@("out\paper_wave\$($w[0])_gc$gc", $gc, $w[1])) } }
"jobs: $($jobs.Count)"
$running = @()
foreach ($j in $jobs) {
    $pre = $j[0]
    if (Test-Path "$pre.log") { continue }
    while (($running | Where-Object { -not $_.HasExited }).Count -ge 28) { Start-Sleep 5 }
    $argl = @("active","$($j[1])","0","0","40000","0.1","$pre.csv") + $kv + $j[2]
    $running += Start-Process -PassThru -WindowStyle Hidden -FilePath "$d\bin\gisim.exe" -WorkingDirectory $d -RedirectStandardOutput "$pre.log" -ArgumentList $argl
}
$running | ForEach-Object { $_.WaitForExit() }
"all done"
