include { agrenseq } from './modules/agrenseq.nf'
include { drenseq } from './modules/drenseq.nf'
include { smrtrenseq } from './modules/smrtrenseq.nf'

workflow {
    main:
    if (params.workflow == "agrenseq") {
        (association, ag_plot, blast_plot, filtered_contigs, candidates_fa, candidates_bed, nlr_candidates) = agrenseq()
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
        (passed, missed, transposed_coverage) = drenseq()
        association = channel.empty()
        ag_plot = channel.empty()
        blast_plot = channel.empty()
        filtered_contigs = channel.empty()
        candidates_fa = channel.empty()
        candidates_bed = channel.empty()
        nlr_candidates = channel.empty()
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
        (assembly, report, stats, annotator_text, annotator_fa, nlr_summary, input_stats, sorted_bed, parsed_coverage) = smrtrenseq()
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
    } else {
        error("Unknown workflow: ${params.workflow}")
    }

    publish:
    association_text = association
    association_plot = ag_plot
    blast_location_plot = blast_plot
    contigs_out = filtered_contigs
    candidates_fasta = candidates_fa
    candidates_location_bed = candidates_bed
    candidates_nlr_positive_fasta = nlr_candidates
    passed_genes = passed
    missed_genes = missed
    coverage = transposed_coverage
    assembled_contigs = assembly
    assembly_report = report
    asssembly_stats = stats
    nlr_annotator_text = annotator_text
    nlr_annotator_fasta = annotator_fa
    nlr_annotator_summary = nlr_summary
    input_reads_stats = input_stats
    nlr_annotator_sorted_bed = sorted_bed
    nlr_coverage = parsed_coverage
}

output {
    association_text {
    }

    association_plot {
    }

    blast_location_plot {
    }

    contigs_out {
    }

    candidates_fasta {
    }

    candidates_location_bed {
    }

    candidates_nlr_positive_fasta {
    }
    passed_genes {
        path 'diagnostics'
    }

    missed_genes {
        path 'diagnostics'
    }

    coverage {
    }

    assembled_contigs {
        path { sample, assembled_contigs -> "${sample}" }
    }

    assembly_report {
        path { sample, assembly_report -> "${sample}" }
    }

    asssembly_stats {
        path { sample, seqkit_out -> "${sample}" }
    }

    nlr_annotator_text {
        path { sample, nlr_annotator_txt -> "${sample}" }
    }

    nlr_annotator_fasta {
        path { sample, nlr_annotator_fa -> "${sample}" }
    }

    nlr_annotator_summary {
        path { sample, summary_of_nlrs -> "${sample}" }
    }

    input_reads_stats {
        path { sample, stats_input -> "${sample}" }
    }

    nlr_annotator_sorted_bed {
        path { sample, sorted_nlr_bed -> "${sample}" }
    }

    nlr_coverage {
        path { sample, parsed_nlr_coverage -> "${sample}" }
    }
}
