include { agrenseq } from './modules/agrenseq.nf'
include { drenseq } from './modules/drenseq.nf'
include { smrtrenseq } from './modules/smrtrenseq.nf'

workflow {
    main:
    if (params.workflow == "agrenseq") {
        association = agrenseq().out.association_txt
        ag_plot = agrenseq().out.association_plot
        blast_plot = agrenseq().out.bl_plot
        filtered_contigs = agrenseq().out.contigs
        candidates_fa = agrenseq().out.cand_fa
        candidates_bed = agrenseq().out.cand_bed
        nlr_candidates = agrenseq().out.cand_nlr_pos
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
        passed = drenseq().out.passed_genes
        missed = drenseq().out.missed_genes
        transposed_coverage = drenseq().out.cov
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
        assembly = smrtrenseq().out.contigs_out
        report = smrtrenseq().out.rep
        stats = smrtrenseq().out.stat
        annotator_text = smrtrenseq().out.ann_txt
        annotator_fa = smrtrenseq().out.ann_fa
        nlr_summary = smrtrenseq().out.nlr_sum
        input_stats = smrtrenseq().out.in_stat
        sorted_bed = smrtrenseq().out.nlr_sort_bed
        parsed_coverage = smrtrenseq().out.cov_parse
    } else {
        error("Unknown workflow: ${params.workflow}")
    }

    publish:
    association = association
    ag_plot = ag_plot
    blast_plot = blast_plot
    filtered_contigs = filtered_contigs
    candidates_fa = candidates_fa
    candidates_bed = candidates_bed
    nlr_candidates = nlr_candidates
    passed = passed
    missed = missed
    transposed_coverage = transposed_coverage
    assembly = assembly
    report = report
    stats = stats
    annotator_text = annotator_text
    annotator_fa = annotator_fa
    nlr_summary = nlr_summary
    input_stats = input_stats
    sorted_bed = sorted_bed
    parsed_coverage = parsed_coverage
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
