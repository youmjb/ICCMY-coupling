#!/bin/sh
# PRC pilot (2026-10-15): isolated ICC-MY (s 0.72), 200-ms rectangular pulses, dt 0.1 ms -> out/paper_prc/
# gL = passive stand-in for the muscle column, g_c*G_load/(g_c+G_load) toward -65.5 mV:
#   0 isolated; 0.60 rabbit calibrated column (g_c 0.96, G_load 1.61); 0.91 mouse (g_c 2.08); 2.0 a large load
cd "$(dirname "$0")/.."
mkdir -p out/paper_prc
for a in 10 20 50 -20 -50; do
  ./bin/gisim.exe iccprc $a 200 40 0.1 out/paper_prc/prc_a${a}_d200.csv trace=0.4
done
for gl in 0 0.6 0.91 2; do
  for a in 25 30 35 40 45 60 80; do
    ./bin/gisim.exe iccprc $a 200 20 0.1 out/paper_prc/load_gL${gl}_a${a}.csv gL=$gl
  done
done
