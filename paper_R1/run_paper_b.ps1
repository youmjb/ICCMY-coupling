# Paper R1, priority 2 follow-up: IJP divider with the ICC-MY rest moved into the recorded range (Kito -65.1, Cai -62.2 mV).
# out\paper_rijp (iccInj -30) put the rest at -33 to -52 mV (beyond the range), so repeat with iccInj -10 and -20 pA
# (rest about -67 and -62 mV in out\paper_rmp). Same protocol as run_v4ijp.ps1: 0.5 s purinergic pulse at 8 phases + control.
#   out\paper_rijp_i10, out\paper_rijp_i20 ; g_c 1/1.5/2, m 0.3, v4 base keys, 30 s
$d = Split-Path -Parent (Split-Path -Parent $MyInvocation.MyCommand.Path)
Set-Location $d
$env:OMP_NUM_THREADS = 1
$base = @("m=0.3","fI=0","fP=0.1","gSP=3","pdModel=1","smcModel=1","gK0run=0.4","iccScale=0.72")
$jobs = @()
foreach ($inj in 10,20) {
    New-Item -ItemType Directory -Force "out\paper_rijp_i$inj" | Out-Null
    foreach ($gc in 1,1.5,2) {
        $jobs += ,@("out\paper_rijp_i$inj\gc$($gc)_ctrl", $gc, 0, @("iccInj=-$inj","injOn=0","injOff=30000"))
        foreach ($t0 in 20000,20450,20900,21350,21800,22250,22700,23150) {
            $jobs += ,@("out\paper_rijp_i$inj\gc$($gc)_p$t0", $gc, 1, @("iccInj=-$inj","injOn=0","injOff=30000","uOn=$t0","uOff=$($t0 + 500)"))
        }
    }
}
$running = @()
foreach ($j in $jobs) {
    $pre = $j[0]
    if (Test-Path "$pre.log") { continue }
    while (($running | Where-Object { -not $_.HasExited }).Count -ge 28) { Start-Sleep 5 }
    $argl = @("active","$($j[1])","$($j[2])","0","30000","0.1","$pre.csv") + $base + $j[3]
    $running += Start-Process -PassThru -WindowStyle Hidden -FilePath "$d\bin\gisim.exe" -WorkingDirectory $d -RedirectStandardOutput "$pre.log" -ArgumentList $argl
}
$running | ForEach-Object { $_.WaitForExit() }
"all done"
