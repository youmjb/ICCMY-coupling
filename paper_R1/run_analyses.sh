#!/bin/sh
# Re-run every analysis and figure of the paper from the saved runs (revision 2: event-aligned amplitudes via amp_event.py).
cd "$(dirname "$0")"
export PYTHONIOENCODING=utf-8
python analyze_criteria.py > /dev/null
for f in analyze_loadrule analyze_paper_a analyze_paper_anat analyze_paper_cal analyze_paper_flc analyze_paper_nif analyze_paper_poh \
         analyze_paper_wave analyze_paper_sens ratio_uncertainty mouse_target analyze_rev1; do python $f.py > $f.out.txt 2>&1 || echo "FAILED $f"; done
python null_model.py > /dev/null; python event_amplitude.py > /dev/null; python icc_loss.py > /dev/null
python fig2_assemble.py > /dev/null && python figs_assemble.py > /dev/null && python graphical_abstract.py > /dev/null && python fig_s3_prc.py > /dev/null
echo done
