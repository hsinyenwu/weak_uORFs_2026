#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""rerun_all_modules.py — runs the main script, the figure driver and every revision module in the order that their inputs need
(2026-11-22; first used for V53, after the correction of the CAGE positions by correct_CAGE_positions.R).
Each step starts when the steps it depends on have ended without error; up to six steps run at the same time.
Logs: output/rerun_logs/<step>.log; status: output/rerun_logs/<step>.json.
Run from the repository root (the folder that contains 20260827_weak_uORFs.R):
    python3 revision_modules/rerun_all_modules.py            # every step that has not ended without error yet
    python3 revision_modules/rerun_all_modules.py A31 A42    # only these steps (the steps they depend on must be done)
    python3 revision_modules/rerun_all_modules.py --fresh    # forget earlier runs and run everything
Order of the modules (a module after the arrow needs an output of the one before it):
    main script -> A2 -> A3 -> A20 -> A11 -> A12 -> A31 -> A42;  A12 -> A13, A38, A39, A40, A41;  A20 -> A29
    A2 -> A14, A19, A25 -> A25_tables, A25_figure, A26 -> Figure 3C panel (CAGE-seq, then TSS-seq)
    main script -> A33 -> A34 -> A18;  A23 -> A24;  A35 -> A37
    main script -> A5, A6, A7, A8, A9, A15, A16, A17, A21, A22, A30, A36, Supplemental Fig. S8 (FigS2_TSSseq_panels.R)
    figure driver (built from the main script by build_revised_figure_driver.py) after the main script."""
import os, subprocess, time, json, sys, threading
DRV = '_to_delete/driver_compact_20261028'; LOG = 'output/rerun_logs'
assert os.path.exists('20260827_weak_uORFs.R') and os.path.isdir('revision_modules'), 'run from the repository root'
os.makedirs(LOG, exist_ok=True); os.makedirs(DRV, exist_ok=True)
M = lambda name: ['Rscript', f'revision_modules/{name}.R']
STEPS = {   # name: (command, steps it depends on, environment)
 'main':   (['Rscript', '20260827_weak_uORFs.R'], [], {}),
 'driver_build': (['python3', 'revision_modules/build_revised_figure_driver.py', '20260827_weak_uORFs.R', f'{DRV}/revised_figure_driver.R', DRV], [], {}),
 'driver': (['Rscript', f'{DRV}/revised_figure_driver.R'], ['driver_build', 'main'], {}),
 'A2': (M('A2_leader_splicing_audit'), ['main'], {}), 'A3': (M('A3_uIF_abundance_and_independence'), ['A2'], {}),
 'A20': (M('A20_flag_neighbouring_gene_CDS_and_table_legends'), ['A3'], {}),
 'A5': (M('A5_riboseq_QC_figure'), ['main'], {}), 'A6': (M('A6_nulls_and_uncertainty'), ['main'], {}), 'A7': (M('A7_ATG_composition_controls'), ['main'], {}),
 'A8': (M('A8_CPuORF_matched_controls'), ['main'], {}), 'A9': (M('A9_GO_figures_FDR'), ['main'], {}),
 'A11': (M('A11_uncapped_vs_uIF'), ['A20'], {}), 'A12': (M('A12_PARE_uncapped_vs_uIF'), ['A11'], {}), 'A13': (M('A13_union_leader_window_sensitivity'), ['A12'], {}),
 'A14': (M('A14_junction_supported_splicing'), ['A2'], {}), 'A15': (M('A15_intronic_TSS_in_leaders'), ['main'], {}), 'A16': (M('A16_predicted_only_uORF_TE'), ['main'], {}),
 'A17': (M('A17_Figure3C_same_gene_mapping'), ['main'], {}), 'A19': (M('A19_TuORF_TSS_vs_splicing_counts'), ['A2'], {}), 'A21': (M('A21_mORF_window_check'), ['main'], {}),
 'A22': (M('A22_uIF_TE_abundance_null'), ['main'], {}), 'A23': (M('A23_GO_pooled_mid_high'), ['main'], {}), 'A24': (M('A24_GO_by_max_pUP_figures_table'), ['A23'], {}),
 'A25': (M('A25_ouORF_leader_splicing_audit'), ['A2'], {}), 'A25_tables': (M('A25_tables'), ['A25', 'A20'], {}), 'A25_figure': (M('A25_figure'), ['A25'], {}),
 'A26': (M('A26_uORF_pool_fraction_TSS_splicing'), ['A25'], {}), 'A29': (M('A29_pUP_TE_abundance_adjusted'), ['A20'], {}), 'A30': (M('A30_same_stop_TuORF_merge_sensitivity'), ['main'], {}),
 'A31': (M('A31_decay_by_inclusion_and_splicing'), ['A12', 'A2'], {}), 'A33': (M('A33_predicted_uORFs_with_overlapping'), ['main'], {}), 'A34': (M('A34_predicted_uORFs_TSS_splicing'), ['A33', 'A2'], {}),
 'A18': (M('A18_mRNA_weighted_TSS_splicing'), ['A34'], {}), 'A35': (M('A35_housekeeping_tau'), ['main'], {}), 'A36': (M('A36_GO_umbrella_terms_low_group'), ['main'], {}),
 'A37': (M('A37_low_group_expression_level'), ['A35', 'A6'], {}), 'A38': (M('A38_matched_decay_comparison'), ['A12', 'A31'], {}), 'A39': (M('A39_uncapped_TE_halflife_partial'), ['A12'], {}),
 'A40': (M('A40_uncapped_halflife_TE_figure'), ['A12'], {}), 'A41': (M('A41_sorenson_halflife_check'), ['A12'], {}), 'A42': (M('A42_FigureS22_Sorenson_uncapped'), ['A31'], {}),
 'FigS2': (M('FigS2_TSSseq_panels'), ['main'], {}), 'Fig3C_CAGE': (M('Fig3C_prevalence_ladder_panel'), ['A26'], {}),
 'Fig3C_TSSseq': (M('Fig3C_prevalence_ladder_panel'), ['Fig3C_CAGE'], {'LADDER_SOURCE': 'TSSseq'}),
}
def status(s):
    p = f'{LOG}/{s}.json'; return json.load(open(p))['status'] if os.path.exists(p) else None
def run(s):
    cmd, deps, env = STEPS[s]; t0 = time.time()
    with open(f'{LOG}/{s}.log', 'w') as lg: p = subprocess.run(cmd, stdout=lg, stderr=subprocess.STDOUT, env=dict(os.environ, **env))
    st = f'exit {p.returncode}'; json.dump({'step': s, 'status': st, 'seconds': round(time.time() - t0), 'ended': time.strftime('%Y-%m-%d %H:%M:%S')}, open(f'{LOG}/{s}.json', 'w'))
    print(f'{time.strftime("%H:%M:%S")} {s:14s} {st:8s} {round(time.time() - t0):5d} s', flush=True)
args = sys.argv[1:]
if '--fresh' in args:
    for s in STEPS:
        if os.path.exists(f'{LOG}/{s}.json'): os.remove(f'{LOG}/{s}.json')
    args = [a for a in args if a != '--fresh']
todo = args or [s for s in STEPS if status(s) != 'exit 0']
for s in todo:
    assert s in STEPS, s
    if os.path.exists(f'{LOG}/{s}.json'): os.remove(f'{LOG}/{s}.json')
lock = threading.Lock(); running = set(); done = set(s for s in STEPS if status(s) is not None); sem = threading.Semaphore(6)
def worker(s):
    with sem: run(s)
    with lock: running.discard(s); done.add(s)
while True:
    with lock:
        left = [s for s in todo if s not in done and s not in running]
        if not left and not running: break
        for s in left:
            deps = STEPS[s][1]
            if any(d in done and status(d) != 'exit 0' for d in deps):
                json.dump({'step': s, 'status': 'not run (a step before it failed)', 'seconds': 0}, open(f'{LOG}/{s}.json', 'w')); done.add(s); print(f'{s}: not run (a step before it failed)', flush=True); continue
            if all(status(d) == 'exit 0' for d in deps): running.add(s); threading.Thread(target=worker, args=(s,), daemon=True).start()
    time.sleep(2)
bad = [s for s in STEPS if status(s) != 'exit 0']; print('ALL DONE' if not bad else 'NOT DONE: ' + ', '.join(bad), flush=True); sys.exit(1 if bad else 0)
