# Paper R1 supplement: ICC-DMP excitatory (cholinergic) input, extra seeds. Same commands as out\dmp1 (2026-09-29 session):
#   active 1.5 0 <aI> 60000 0.1 <csv> + v4 base keys + imModel=2 fI=3 [L=10|15]  (aI 0 = DMP baseline, 9 = ACh x10)
# seed 1 = out\dmp1\{L10_,L15_,}f3_a{0,9}; this script adds seeds 2 and 3. Output: out\paper_dmp\L<L>_a<aI>_s<seed>
$d = Split-Path -Parent (Split-Path -Parent $MyInvocation.MyCommand.Path)
Set-Location $d
$env:OMP_NUM_THREADS = 1
New-Item -ItemType Directory -Force out\paper_dmp | Out-Null
$kv = @("m=0.3","fP=0.1","gSP=3","pdModel=1","smcModel=1","gK0run=0.4","iccScale=0.72","imModel=2","fI=3")
$running = @()
foreach ($seed in 2,3) { foreach ($L in 10,15,30) { foreach ($a in 0,9) {
    $pre = "out\paper_dmp\L$($L)_a$($a)_s$seed"
    if (Test-Path "$pre.log") { continue }
    $argl = @("active","1.5","0","$a","60000","0.1","$pre.csv") + $kv + @("L=$L","seed=$seed")
    $running += Start-Process -PassThru -WindowStyle Hidden -FilePath "$d\bin\gisim.exe" -WorkingDirectory $d -RedirectStandardOutput "$pre.log" -ArgumentList $argl
} } }
$running | ForEach-Object { $_.WaitForExit() }
"all done"
