# Paper R1, priority 3: one-at-a-time sensitivity around the v4 base (m 0.3, gSS 30, L 30, gK0run 0.4, fP 0.1, gSP 3, iccGL 0, seed 1).
# Per condition: g_c sweep 0/0.5/1/1.5/2/3 nS (40 s, no input) -> g_c matching the Kito CSMC amplitude;
#                0.5 s purinergic pulse at 8 phases + control at g_c 1.5 (30 s) -> IJP divider.
# iccGL = extra passive conductance on ICC-MY (nS, E = -72 mV; key added 2026-10-14) = the ICC-MY own-conductance axis (base ~3.6-3.7 nS).
# Base condition itself: out\v4win (g_c sweep) and out\v4ijp (IJP) are reused.  Output: out\paper_sens\<cond>\
$d = Split-Path -Parent (Split-Path -Parent $MyInvocation.MyCommand.Path)
Set-Location $d
$env:OMP_NUM_THREADS = 1
$base = [ordered]@{ m = "0.3"; gSS = "30"; L = "30"; gK0run = "0.4"; fP = "0.1"; gSP = "3"; iccGL = "0"; seed = "1" }
$conds = @(
    @("gSS", "15"), @("gSS", "60"), @("L", "15"), @("L", "45"), @("m", "0.15"), @("m", "0.45"),
    @("gK0run", "0.2"), @("gK0run", "0.6"), @("fP", "0.05"), @("fP", "0.15"), @("gSP", "1.5"), @("gSP", "4.5"),
    @("iccGL", "-0.6"), @("iccGL", "1.8"), @("iccGL", "3.6"), @("seed", "2"), @("seed", "3")
)
$jobs = @()
foreach ($c in $conds) {
    $p = [ordered]@{}; foreach ($k in $base.Keys) { $p[$k] = $base[$k] }; $p[$c[0]] = $c[1]
    $tag = "$($c[0])$($c[1])"; $dir = "out\paper_sens\$tag"
    New-Item -ItemType Directory -Force $dir | Out-Null
    $kv = @("fI=0","pdModel=1","smcModel=1","iccScale=0.72") + ($p.Keys | ForEach-Object { "$_=$($p[$_])" })
    foreach ($gc in 0,0.5,1,1.5,2,3) { $jobs += ,@("$dir\gc$($gc)_u0", $gc, 0, 40000, $kv) }
    $jobs += ,@("$dir\ijp_ctrl", 1.5, 0, 30000, $kv)
    foreach ($t0 in 20000,20450,20900,21350,21800,22250,22700,23150) { $jobs += ,@("$dir\ijp_p$t0", 1.5, 1, 30000, ($kv + @("uOn=$t0","uOff=$($t0 + 500)"))) }
}
$running = @()
foreach ($j in $jobs) {
    $pre = $j[0]
    if (Test-Path "$pre.log") { continue }
    while (($running | Where-Object { -not $_.HasExited }).Count -ge 28) { Start-Sleep 5 }
    $argl = @("active","$($j[1])","$($j[2])","0","$($j[3])","0.1","$pre.csv") + $j[4]
    $running += Start-Process -PassThru -WindowStyle Hidden -FilePath "$d\bin\gisim.exe" -WorkingDirectory $d -RedirectStandardOutput "$pre.log" -ArgumentList $argl
}
$running | ForEach-Object { $_.WaitForExit() }
"all done"
