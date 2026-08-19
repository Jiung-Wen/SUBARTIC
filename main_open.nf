nextflow.enable.dsl=2

workflow {
    // 1. Inputs
    manifest_ch = Channel.fromPath(params.manifest, checkIfExists: true)
    ref_seqs_ch = Channel.fromPath(params.ref_seqs, checkIfExists: true)

    // 2. Workflow Logic
    IMPORT_DEMUX(manifest_ch)
    MERGE_PAIRS(IMPORT_DEMUX.out.demux)
    QUALITY_FILTER(MERGE_PAIRS.out.merged)
    DEREPLICATE(QUALITY_FILTER.out.filtered)
    FILTER_FEATURES(DEREPLICATE.out.table, DEREPLICATE.out.seqs)
    CLUSTER_FEATURES(
        FILTER_FEATURES.out.filt_seqs, 
        FILTER_FEATURES.out.filt_table, 
        ref_seqs_ch
    )
}

/* * Process 1: Import 
 * Output Prefix: 01_
 */
process IMPORT_DEMUX {
    publishDir "${params.outdir}/01_demux", mode: 'copy'

    input:
    path manifest

    output:
    path "01_paired-end-demux.qza", emit: demux
    path "01_demux.qzv"

    script:
    """
    qiime tools import \
      --type 'SampleData[PairedEndSequencesWithQuality]' \
      --input-path ${manifest} \
      --output-path 01_paired-end-demux.qza \
      --input-format PairedEndFastqManifestPhred33V2

    qiime demux summarize \
      --i-data 01_paired-end-demux.qza \
      --p-n 10000 \
      --o-visualization 01_demux.qzv
    """
}

/* * Process 2: Merge 
 * Output Prefix: 02_
 */
process MERGE_PAIRS {
    publishDir "${params.outdir}/02_merged", mode: 'copy'
    cpus params.threads

    input:
    path demux

    output:
    path "02_joined-reads.qza", emit: merged
    path "02_unmerged-reads.qza"

    script:
    """
    qiime vsearch merge-pairs \
      --i-demultiplexed-seqs ${demux} \
      --p-minovlen 20 \
      --p-minmergelen ${params.min_mergelen} \
      --p-maxee 2 \
      --p-threads ${task.cpus} \
      --o-merged-sequences 02_joined-reads.qza \
      --o-unmerged-sequences 02_unmerged-reads.qza
    """
}

/* * Process 3: Quality Filter 
 * Output Prefix: 03_
 */
process QUALITY_FILTER {
    publishDir "${params.outdir}/03_quality", mode: 'copy'

    input:
    path joined_reads

    output:
    path "03_joined-filtered.qza", emit: filtered
    path "03_joined-filter-stats.qza"
    path "03_joined-filter-stats.qzv"

    script:
    """
    qiime quality-filter q-score \
      --i-demux ${joined_reads} \
      --o-filtered-sequences 03_joined-filtered.qza \
      --o-filter-stats 03_joined-filter-stats.qza

    qiime metadata tabulate \
      --m-input-file 03_joined-filter-stats.qza \
      --o-visualization 03_joined-filter-stats.qzv
    """
}

/* * Process 4: Dereplicate 
 * Output Prefix: 04_
 */
process DEREPLICATE {
    publishDir "${params.outdir}/04_derep", mode: 'copy'

    input:
    path filtered_reads

    output:
    path "04_table.qza", emit: table
    path "04_rep-seqs.qza", emit: seqs
    path "04_table.qzv"
    path "04_rep-seqs.qzv"

    script:
    """
    qiime vsearch dereplicate-sequences \
      --i-sequences ${filtered_reads} \
      --p-min-seq-length ${params.min_mergelen} \
      --o-dereplicated-table 04_table.qza \
      --o-dereplicated-sequences 04_rep-seqs.qza

    qiime feature-table summarize \
      --i-table 04_table.qza \
      --o-visualization 04_table.qzv

    qiime feature-table tabulate-seqs \
      --i-data 04_rep-seqs.qza \
      --o-visualization 04_rep-seqs.qzv
    """
}

/* * Process 5: Filter Features 
 * Output Prefix: 05_
 */
process FILTER_FEATURES {
    publishDir "${params.outdir}/05_filtered_features", mode: 'copy'

    input:
    path table
    path seqs

    output:
    path "05_filtered-table.qza", emit: filt_table
    path "05_filtered-rep-seqs.qza", emit: filt_seqs
    path "05_filtered-table.qzv"
    path "05_filtered-rep-seqs.qzv"

    script:
    """
    # Filter Table
    qiime feature-table filter-features \
      --i-table ${table} \
      --p-min-frequency ${params.min_frequency} \
      --o-filtered-table 05_filtered-table.qza

    # Filter Seqs based on Table
    qiime feature-table filter-seqs \
      --i-data ${seqs} \
      --i-table 05_filtered-table.qza \
      --o-filtered-data 05_filtered-rep-seqs.qza

    qiime feature-table summarize \
      --i-table 05_filtered-table.qza \
      --o-visualization 05_filtered-table.qzv

    qiime feature-table tabulate-seqs \
      --i-data 05_filtered-rep-seqs.qza \
      --o-visualization  05_filtered-rep-seqs.qzv
    """
}

/* * Process 6: Cluster 
 * Output Prefix: 06_
 */
process CLUSTER_FEATURES {
    publishDir "${params.outdir}/06_clustering", mode: 'copy'
    cpus params.threads

    input:
    path seqs
    path table
    path reference_db

    output:
    path "06_open_ref_table_*.qza"
    path "06_open_ref_seqs_*.qza"
    path "06_open_ref_new_reference_*.qza"
    path "*.qzv"

    script:
    def perc = params.cluster_perc
    """
    qiime vsearch cluster-features-open-reference \
      --i-sequences ${seqs} \
      --i-table ${table} \
      --i-reference-sequences ${reference_db} \
      --p-perc-identity ${perc} \
      --p-strand both \
      --p-threads ${task.cpus} \
      --o-clustered-table 06_open_ref_table_${perc}.qza \
      --o-clustered-sequences 06_open_ref_seqs_${perc}.qza \
      --o-new-reference-sequences 06_open_ref_new_reference_${perc}.qza

    # Generate Visualizations
    qiime feature-table summarize \
      --i-table 06_open_ref_table_${perc}.qza \
      --o-visualization 06_open_ref_table_${perc}.qzv

    qiime feature-table tabulate-seqs \
      --i-data 06_open_ref_seqs_${perc}.qza \
      --o-visualization 06_open_ref_seqs_${perc}.qzv
      
    qiime feature-table tabulate-seqs \
      --i-data 06_open_ref_new_reference_${perc}.qza \
      --o-visualization 06_open_ref_new_reference_${perc}.qzv
    """
}
