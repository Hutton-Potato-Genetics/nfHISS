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
    association_txt = association
    association_plot = ag_plot
    bl_plot = blast_plot
    contigs = filtered_contigs
    cand_fa = candidates_fa
    cand_bed = candidates_bed
    cand_nlr_pos = nlr_candidates
    passed_genes = passed
    missed_genes = missed
    cov = transposed_coverage
    contigs_out = assembly
    rep = report
    stat = stats
    ann_txt = annotator_text
    ann_fa = annotator_fa
    nlr_sum = nlr_summary
    in_stat = input_stats
    nlr_sort_bed = sorted_bed
    cov_parse = parsed_coverage
}

output {
    contigs {
        path { sample, assembled_contigs -> "${sample}" }
    }

    rep {
        path { sample, assembly_report -> "${sample}" }
    }

    stat {
        path { sample, seqkit_out -> "${sample}" }
    }

    ann_txt {
        path { sample, nlr_annotator_txt -> "${sample}" }
    }

    ann_fa {
        path { sample, nlr_annotator_fa -> "${sample}" }
    }

    nlr_sum {
        path { sample, summary_of_nlrs -> "${sample}" }
    }

    in_stat {
        path { sample, stats_input -> "${sample}" }
    }

    nlr_sort_bed {
        path { sample, sorted_nlr_bed -> "${sample}" }
    }

    cov_parse {
        path { sample, parsed_nlr_coverage -> "${sample}" }
    }

    passed_genes {
        path 'diagnostics'
    }

    missed_genes {
        path 'diagnostics'
    }

    cov {
    }

    association_txt {
    }

    association_plot {
    }

    bl_plot {
    }

    contigs_out {
    }

    cand_fa {
    }

    cand_bed {
    }

    cand_nlr_pos {
    }
}
