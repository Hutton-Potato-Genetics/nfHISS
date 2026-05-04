include { agrenseq } from './modules/agrenseq.nf'
include { drenseq } from './modules/drenseq.nf'
include { smrtrenseq } from './modules/smrtrenseq.nf'

workflow {
    main:
    if (params.workflow == "agrenseq") {
        association = agrenseq().association_txt
        ag_plot = agrenseq().association_plot
        blast_plot = agrenseq().bl_plot
        filtered_contigs = agrenseq().contigs
        candidates_fa = agrenseq().cand_fa
        candidates_bed = agrenseq().cand_bed
        nlr_candidates = agrenseq().cand_nlr_pos
        passed = channel.empty()
        missed = channel.empty()
        transposed_coverage = channel.empty()
        assembly = channel.empty()
        report = channel.empty()
        stats = channel.empty()
        annotator_text = channel.empty()
        annotator_fa = channel.empty()
        nlr_summary = channel.empty()
        input_stats = channel.empty()
        sorted_bed = channel.empty()
        parsed_coverage = channel.empty()
    } else if (params.workflow == "drenseq") {
        association = channel.empty()
        ag_plot = channel.empty()
        blast_plot = channel.empty()
        filtered_contigs = channel.empty()
        candidates_fa = channel.empty()
        candidates_bed = channel.empty()
        nlr_candidates = channel.empty()
        passed = drenseq().passed_genes
        missed = drenseq().missed_genes
        transposed_coverage = drenseq().cov
        assembly = channel.empty()
        report = channel.empty()
        stats = channel.empty()
        annotator_text = channel.empty()
        annotator_fa = channel.empty()
        nlr_summary = channel.empty()
        input_stats = channel.empty()
        sorted_bed = channel.empty()
        parsed_coverage = channel.empty()
    } else if (params.workflow == "smrtrenseq") {
        association = channel.empty()
        ag_plot = channel.empty()
        blast_plot = channel.empty()
        filtered_contigs = channel.empty()
        candidates_fa = channel.empty()
        candidates_bed = channel.empty()
        nlr_candidates = channel.empty()
        passed = channel.empty()
        missed = channel.empty()
        transposed_coverage = channel.empty()
        assembly = smrtrenseq().contigs_out
        report = smrtrenseq().rep
        stats = smrtrenseq().stat
        annotator_text = smrtrenseq().ann_txt
        annotator_fa = smrtrenseq().ann_fa
        nlr_summary = smrtrenseq().nlr_sum
        input_stats = smrtrenseq().in_stat
        sorted_bed = smrtrenseq().nlr_sort_bed
        parsed_coverage = smrtrenseq().cov_parse
    } else {
        error("Unknown workflow: ${params.workflow}")
    }

    publish:
    association_txt = association
    association_plot = ag_plot
    bl_plot = blast_plot
    contigs_out = filtered_contigs
    cand_fa = candidates_fa
    cand_bed = candidates_bed
    cand_nlr_pos = nlr_candidates
    passed_genes = passed
    missed_genes = missed
    cov = transposed_coverage
    contigs = assembly
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
    passed_genes {
        path 'diagnostics'
    }

    missed_genes {
        path 'diagnostics'
    }

    cov {
    }
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
