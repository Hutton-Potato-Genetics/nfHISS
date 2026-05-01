include { agrenseq } from './modules/agrenseq.nf'
include { drenseq } from './modules/drenseq.nf'
include { smrtrenseq } from './modules/smrtrenseq.nf'

workflow {
    main:
    if (params.workflow == "agrenseq") {
        agrenseq()
    } else if (params.workflow == "drenseq") {
        drenseq()
    } else if (params.workflow == "smrtrenseq") {
        smrtrenseq()
    } else {
        error("Unknown workflow: ${params.workflow}")
    }

    publish:
    contigs = assembly
    rep = report
    stat = stats
    ann_txt = annotator_text
    ann_fa = annotator_fa
    nlr_sum = nlr_summary
    in_stat = input_stats
    nlr_sort_bed = sorted_bed
    cov_parse = parsed_coverage
    passed_genes = passed
    missed_genes = missed
    cov = transposed_coverage
    association_txt = association
    association_plot = ag_plot
    bl_plot = blast_plot
    contigs_out = filtered_contigs
    cand_fa = candidates_fa
    cand_bed = candidates_bed
    cand_nlr_pos = nlr_candidates
}
