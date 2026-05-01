include { agrenseq } from './modules/agrenseq.nf'
include { drenseq } from './modules/drenseq.nf'
include { smrtrenseq } from './modules/smrtrenseq.nf'

workflow {
    main:
    if (params.workflow == "agrenseq") {
        ag_results = agrenseq()
    } else if (params.workflow == "drenseq") {
        dren_results = drenseq()
    } else if (params.workflow == "smrtrenseq") {
        smrt_results = smrtrenseq()
    } else {
        error("Unknown workflow: ${params.workflow}")
    }

    publish:
    if (params.workflow == "agrenseq") {
        association = ag_results.association_txt
        ag_plot = ag_results.association_plot
        blast_plot = ag_results.bl_plot
        filtered_contigs = ag_results.contigs
        candidates_fa = ag_results.cand_fa
        candidates_bed = ag_results.cand_bed
        nlr_candidates = ag_results.cand_nlr_pos
    } else if (params.workflow == "drenseq") {
        passed = dren_results.passed_genes
        missed = dren_results.missed_genes
        transposed_coverage = dren_results.cov
    } else if (params.workflow == "smrtrenseq") {
        assembly = smrt_results.contigs_out
        report = smrt_results.rep
        stats = smrt_results.stat
        annotator_text = smrt_results.ann_txt
        annotator_fa = smrt_results.ann_fa
        nlr_summary = smrt_results.nlr_sum
        input_stats = smrt_results.in_stat
        sorted_bed = smrt_results.nlr_sort_bed
        parsed_coverage = smrt_results.cov_parse
    }
}

output {
    if (params.workflow == "agrenseq") {
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
    } else if (params.workflow == "drenseq") {
        passed_genes {
            path 'diagnostics'
        }

        missed_genes {
            path 'diagnostics'
        }

        cov {
        }
    } else if (params.workflow == "smrtrenseq") {
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
    }
}
