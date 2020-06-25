#!/usr/bin/env nextflow
/*
 ========================================================================================
 MARVEL
 ========================================================================================
 MARVEL Analysis Pipeline.
 #### Homepage / Documentation
 https://github.com/fuxialexander/marvel
 ----------------------------------------------------------------------------------------
 */
nextflow.preview.dsl = 2
include {cli_banner; checkHostname; create_workflow_summary; get_software_versions} from './modules/workflow_helper'
// TODO Clear iGenome in readme and usage
def helpMessage() {
    // TODO nf-core: Add to this help message with new command line parameters
    log.info cli_banner()
    log.info"""

    Usage:

    The typical command for running the pipeline is as follows:

    nextflow run marvel --reads '*_R{1,2}.fastq.gz' -profile docker

    Mandatory arguments:
      -profile                      Configuration profile to use. Can use multiple (comma separated)
                                    Available: conda, docker, singularity, awsbatch, test and more.
      --enhancer                    Bed file of enhancers
      --promoter                    Bed file of promoters
      --vcf                         VCF file for genotyped variants
      --vcf_index                   .tbi index of VCF file 
      --phenocov                    Sample phenotype and covariates file
      --fasta                       Reference FASTA
      --chunk_size                  Number of regions in each chunk
      --permutation_multiplier      Number of permutations to do
      --motif_scan_threshold        P-value threshold for motif scanning
    Options:
      --genome                      Name of genome reference

    References                      If not specified in the configuration file or you wish to overwrite any of the references.
      --fasta                       Path to Fasta reference

    Other options:
      --outdir                      The output directory where the results will be saved
      --email                       Set this parameter to your e-mail address to get a summary e-mail with details of the run sent to you when the workflow exits
      --email_on_fail               Same as --email, except only send mail if the workflow is not successful
      -name                         Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic.

    AWSBatch options:
      --awsqueue                    The AWSBatch JobQueue that needs to be set when running on AWSBatch
      --awsregion                   The AWS Region for your AWS Batch job to run on
    """.stripIndent()
}

// Show help message
if (params.help) {
    helpMessage()
    exit 0
}

/*
 * SET UP CONFIGURATION VARIABLES
 */

// Check if genome exists in the config file
if (params.genomes && params.genome && !params.genomes.containsKey(params.genome)) {
    exit 1, "The provided genome '${params.genome}' is not available in the genomes file. Currently the available genomes are ${params.genomes.keySet().join(", ")}"
}
params.fasta = params.genome ? params.genomes[ params.genome ].fasta ?: false : false

// Has the run name been specified by the user?
//  this has the bonus effect of catching both -name and --name
custom_runName = params.name
if (!(workflow.runName ==~ /[a-z]+_[a-z]+/)) {
    custom_runName = workflow.runName
}

if ( workflow.profile == 'awsbatch') {
    // AWSBatch sanity checking
    if (!params.awsqueue || !params.awsregion) exit 1, "Specify correct --awsqueue and --awsregion parameters on AWSBatch!"
    // Check outdir paths to be S3 buckets if running on AWSBatch
    // related: https://github.com/nextflow-io/nextflow/issues/813
    if (!params.outdir.startsWith('s3:')) exit 1, "Outdir not on S3 - specify S3 Bucket to run on AWSBatch!"
    // Prevent trace files to be stored on S3 since S3 does not support rolling files.
    if (workflow.tracedir.startsWith('s3:')) exit 1, "Specify a local tracedir or run without trace! S3 cannot be used for tracefiles."
}

// Stage config files
// ch_output_docs = file("$baseDir/docs/output.md", checkIfExists: true)

ch_output_docs = file("$baseDir/docs/output.md", checkIfExists: true)

// Header log info
log.info cli_banner()
def summary = [:]
if (workflow.revision) summary['Pipeline Release'] = workflow.revision
summary['Run Name']         = custom_runName ?: workflow.runName
// TODO: Report custom parameters here
summary['Enhancer'] = params.enhancer
summary['Promoter'] = params.promoter
summary['Variant'] = params.vcf
summary['Variant index'] = params.vcf_index
summary['Phenotype'] = params.phenocov
summary['Region split chunk size'] = params.chunk_size
summary['Permuation times'] = params.permutation_multiplier
summary['Motif scanning threshold'] = params.motif_scan_threshold
summary['Fasta Ref']        = params.fasta
summary['Max Resources']    = "$params.max_memory memory, $params.max_cpus cpus, $params.max_time time per job"
if (workflow.containerEngine) summary['Container'] = "$workflow.containerEngine - $workflow.container"
summary['Output dir']       = params.outdir
summary['Launch dir']       = workflow.launchDir
summary['Working dir']      = workflow.workDir
summary['Script dir']       = workflow.projectDir
summary['User']             = workflow.userName
if (workflow.profile == 'awsbatch') {
    summary['AWS Region']     = params.awsregion
    summary['AWS Queue']      = params.awsqueue
}
summary['Config Profile'] = workflow.profile
if (params.config_profile_description) summary['Config Description'] = params.config_profile_description
if (params.config_profile_contact)     summary['Config Contact']     = params.config_profile_contact
if (params.config_profile_url)         summary['Config URL']         = params.config_profile_url
if (params.email || params.email_on_fail) {
    summary['E-mail Address']    = params.email
    summary['E-mail on failure'] = params.email_on_fail
}
log.info summary.collect { k,v -> "${k.padRight(18)}: $v" }.join("\n")
log.info "-\033[2m--------------------------------------------------\033[0m-"

// Check the hostnames against configured profiles
checkHostname()


/*
 * Parse software version numbers
 */
process get_software_versions {
    publishDir "${params.outdir}/pipeline_info", mode: 'copy',
        saveAs: { filename ->
        if (filename.indexOf(".csv") > 0) filename
        else null
    }

    output:
    file 'software_versions_mqc.yaml' into software_versions_yaml
    file "software_versions.csv"

    script:
    // TODO nf-core: Get all tools to print their version number here
    """
    echo $workflow.manifest.version > v_pipeline.txt
    echo $workflow.nextflow.version > v_nextflow.txt
    samtools --version | head -n 1 > v_samtools.txt
    bcftools --version | head -n 1 > v_bcftools.txt
    scrape_software_versions.py &> software_versions_mqc.yaml
    """
}


/*
 * STEP 1 - Get Region Fasta
 */
process conda_setup_dep {
    input: 
    file glmpath
    
    output:
    val(1)

    script:
    """
    R -e "install.packages(c('glmpath'), dependencies=TRUE, repos='http://cran.rstudio.com/')" \
      && R CMD INSTALL $glmpath
    pip install https://github.com/khramts/assocplots/archive/master.zip
    """
}

process get_region_ref_fasta {
    tag "${name}"
    label 'process_light'
    publishDir "${params.outdir}/sequences"

    input:
    tuple val(ch_name), file(genome_fa), val(name), file(bed)

    output:
    tuple val(name), file("${name}.fa")

    script:
    """
    awk '{print \$1":"(\$2+1)"-"\$3}' $bed | xargs samtools faidx ${genome_fa} > ${name}.fa \
    """
}

process index_fasta {
    tag "$name"
    label 'process_light'
    publishDir "${params.outdir}/sequences"

    input:
    tuple val(name), file(fa)

    output:
    tuple val(name), file(fa), file("${fa}.fai")

    script:
    """
    samtools faidx ${fa} > ${fa}.fai
    """
}

process gunzip {
    label 'process_light'
    publishDir "${params.outdir}/sequences"

    input:
    file gz

    output:
    file "$gz.baseName"

    script:
    """
    gunzip -k --verbose --stdout --force $gz > $gz.baseName
    """
}

process count_nucleotide_frequency {
    tag "${name}"
    label 'process_light'
    publishDir "${params.outdir}/sequences"

    input:
    tuple val(name), file(fasta)

    output:
    tuple val(name), file("${name}.nuc_freq.txt")

    script:
    """
    count_fasta_base.py -o ACGT <(grep -v ">" $fasta | cat <(echo ">") -) > ${name}.nuc_freq.txt
    """
}

process get_hocomoco_motif_in_jaspar_format {
    tag "${name}"
    label 'process_light'
    publishDir "${params.outdir}/motifs"

    output:
    file "*.pwm"
    // file "motif_index.txt"
    // file "motif_symbol.txt"
    // TODO: Consider add option to use HOCOMOCO Core, or even other motif databases

    script:
    """
    get_hocomoco_motif.sh
    """
}

process chunkify_regions {
    tag "${name}"
    label 'process_light'
    publishDir "${params.outdir}/regions/chunks"

    input:
    tuple val(name), file(regions)

    output:
    file "${name}_chunk_*.bed" 
    // TODO: potentially add options to randomize the region file first
    script:
    """
    split --additional-suffix=.bed -l ${params.chunk_size} $regions ${name}_chunk_
    """
}

process get_region_chunk_vcf_indexed {
    tag "${name}"
    label 'process_light'
    publishDir "${params.outdir}/variants/chunks"

    input:
    tuple val(name), file(region_chunk), file(vcf), file(vcf_index)

    output:
    file("${name}.vcf.gz")
    file("${name}.vcf.gz.tbi")

    script:
    """
    cat <(bcftools view -h $vcf) \
        <(bcftools view -H -Ov -R $region_chunk $vcf | sort -k1,1d -k2,2n ) | \
          bgzip > ${name}.vcf.gz && \
          tabix -f -p vcf ${name}.vcf.gz
    """
}

process get_sample_fasta {
    tag "${name}"
    label 'process_high'
    cache 'lenient'
    publishDir "${params.outdir}/sequences/chunks"

    input:
    tuple val(name), file(bed), file(vcf), file(vcf_index), val(genome_name), file(fasta), file(fasta_index), file(phenocov)

    output:
    file("${name}.fa")

    script:
    """
    for sample in `tail -n +2 $phenocov | cut -f1 | sort | tr '\n' '\t' `; do \
        bcftools view --trim-alt-alleles -s \$sample $vcf \
        --min-ac=1 --no-update -Ou \
        | bcftools norm -Ou -m-any -f $fasta \
        | bcftools norm -d both -Ov \
        | awk '\$5!~"<*:"{print}' \
        | bgzip > ${name}.norm.vcf.gz;

        tabix -f -p vcf ${name}.norm.vcf.gz;

        awk '{print \$1":"(\$2)+1"-"\$3}' $bed \
        | xargs samtools faidx $fasta \
        | bcftools consensus -H A ${name}.norm.vcf.gz \
        | sed "s/^>/>\$sample /" >> ${name}.fa ;
    done
    """
}

process get_promoter_enhancer_pair {
    tag "get_promoter_enhancer_pair"
    label 'process_medium'
    publishDir "${params.outdir}/regions/"
    input:
    tuple(pname, promoter)
    tuple(ename, enhancer)
    file(chrom_size)

    output:
    file "promoter_enhancer_pair.txt"

    script:
    """
    awk '{OFS="\t"; print \$1,\$2,\$3,\$1":"\$2"-"\$3}' $promoter \
    | bedtools slop -b $params.promoter_enhancer_distance_threshold -g $chrom_size -i - \
    | sort -k1,1V -k2,2n \
    | bedtools intersect -a - -b $enhancer -wa -wb \
    | cut -f4,8 \
    | awk -F "[\t:-]" '{print \$0, \$5-\$2}' > promoter_enhancer_pair.txt
    """
}

process get_promoter_gene_pair {
    tag "get_promoter_gene_pair"
    label 'process_medium'
    publishDir "${params.outdir}/regions/"
    input:
    tuple(val(name), file(promoter))

    output:
    file "promoter_gene_pair.txt"

    script:
    """
    awk '{OFS="\t"; print \$1":"\$2"-"\$3, \$4}' $promoter > promoter_gene_pair.txt
    """
}



process scan_test {
    // first scan motif then test, useful in region-level test
    tag "${name}"
    label 'process_heavy'
    publishDir "${params.outdir}/test/${test_name}_chunks"
    input:
    tuple val(name), file(fa), file(phenocov), file(motifs), val(test_name), file(nuc_freq), val(prepared)

    output:
    file '*.npz'

    script:
    """
    scan_test.py -m $motifs/*.pwm -s $fa \
    -P $phenocov -o $name -p $params.motif_scan_threshold --batch \
    --permutation-size-multiplier $params.permutation_multiplier \
    --bg \$(cat $nuc_freq | tr '\\n' ' ')
    """
}

process collect_chunk {
    tag "${name}"
    label 'process_heavy'
    publishDir "${params.outdir}/test/"
    input:
    tuple val(name), val(chunk_path)

    output:
    file '*.npz'

    script:
    """
    collect_results.py $name $chunk_path ./
    """
}

process summarize {
    // Summarize result and plot qq-plots, AUROC, etc
    tag "${name}"
    label 'process_medium'
    publishDir "${params.outdir}/summary/"
    input:
    tuple val(name), file(results)
    file phenocov

    output:
    file "*.pdf"
    file "*.xlsx"

    script:
    """
    summarize.py ${name}*real.npz ${name}*profile.npz
    """
}

// workflow summary {
//     get: result_path
//     main:
//     results = Channel.fromPath(result_path)

// }

process profile_test {
    // test with out scan, useful for gene-level analysis
    tag "${name}"
    label 'process_heavy'
    publishDir "${params.outdir}/test/chunks"
    input:
    tuple val(name), file(profiles)
    file phenocov
    tuple val(test_name), file(nuc_freq)

    output:
    file "*.npz"

    script:
    """
    test.py -P $phenocov -o $name --permutation-size-multiplier $params.permutation_multiplier 
    """
}

// workflow gene_test {
//   get: tuple enhancer, promoter
//   main:
//   pe_pair = get_promoter_enhancer_pair(enhancer, promoter, chrom_size)
//   pg_pair = get_promoter_gene_pair(promoter)
  
// }

workflow {
    main:
    // Region-based test
    if ("$workflow.profile" =~ /conda/) {
        glmpath = Channel.fromPath( "glmpath_0.98.tar.gz" )
        prepared = conda_setup_dep(glmpath)
    } else {
        prepared = Channel.value(1)
    }

    vcf = Channel
        .fromPath( params.vcf, checkIfExists: true )
        .ifEmpty { exit 1, "VCF file not found: ${params.vcf}" }

    vcf_index = Channel
        .fromPath( params.vcf_index, checkIfExists: true )
        .ifEmpty { exit 1, "VCF file not found: ${params.vcf_index}" }

    phenocov = Channel
        .fromPath( params.phenocov, checkIfExists: true )
        .ifEmpty { exit 1, "Phenotype & covariates file not found: ${params.phenocov}" }

    if (params.fasta) { 
        ch_fasta = Channel
            .fromPath(params.fasta, checkIfExists: true)
            .ifEmpty { exit 1, "Reference FASTA file not found: ${params.fasta}" }
            .map { file -> tuple(file.baseName, file) }
    }

    if (params.enhancer) {
        enhancer = Channel
            .fromPath( params.enhancer, checkIfExists: true )
            .ifEmpty { exit 1, "Enhancer BED file not found: ${params.enhancer}" }
            .map { file -> tuple(file.baseName, file) }
        regions = enhancer
    }

    if (params.promoter) {
        promoter = Channel
            .fromPath( params.promoter, checkIfExists: true )
            .ifEmpty { exit 1, "Promoter BED file not found: ${params.promoter}" }
            .map { file -> tuple(file.baseName, file) }
        regions = promoter
    }

    if (params.other_regions) {
        other_regions = Channel
            .fromPath( params.other_regions, checkIfExists: true )
            .ifEmpty { exit 1, "Other regions BED files not found: ${params.regions}" }
            .map { file -> tuple(file.baseName, file) }
        
        regions = other_regions.concat(enhancer).concat(promoter).unique()
    } else {
        regions = enhancer.concat(promoter)
    }

    chrom_size = Channel
        .fromPath( "test_input/hg19.chrom.sizes" )

    promoter_enhancer_pair = get_promoter_enhancer_pair(promoter, enhancer, chrom_size)
    promoter_gene_pair = get_promoter_gene_pair(promoter)

    motifs = get_hocomoco_motif_in_jaspar_format()
        .flatten()
        .first()
        .map { file -> file.getParent() }
    regions_fasta = get_region_ref_fasta(ch_fasta.combine(regions))
    fasta_indexed = index_fasta(ch_fasta)
    nuc_freq = count_nucleotide_frequency(regions_fasta)

    chunks_bed = chunkify_regions(regions)
        .flatten()
        .map { file -> tuple(file.baseName, file) }
    chunks_vcf = get_region_chunk_vcf_indexed( chunks_bed.combine( vcf ).combine( vcf_index ) )
    chunks_bed_vcf = chunks_bed.join (
        chunks_vcf[0]
            .map { file -> tuple(file.simpleName, file) }
            .merge(chunks_vcf[1])
    ).combine(fasta_indexed).combine(phenocov) 
    chunks_sample_fasta = get_sample_fasta(chunks_bed_vcf)
    results = scan_test {
        chunks_sample_fasta
            .map { file -> tuple(file.baseName, file) }
            .combine(phenocov)
            .combine(motifs)
            .combine(nuc_freq)
            .combine(prepared)
            .filter { it[0].startsWith(it[4]) }
    }.unique()

    chunk_results = regions.map{
        x -> tuple(x[0], "$baseDir/${params.outdir}/test/" + x[0] + "_chunks/")
    }

    collect_chunk(chunk_results)

    // Genes-based test
    // if (enhancer_profiles && promoter_profiles) {
    //   gene_test(enhancers, promoters)
    // }
}

/*
 * Completion e-mail notification
 */
workflow.onComplete {

    // Set up the e-mail variables
    def subject = "[MARVEL] Successful: $workflow.runName"
    if (!workflow.success) {
        subject = "[MARVEL] FAILED: $workflow.runName"
    }
    def email_fields = [:]
    email_fields['version'] = workflow.manifest.version
    email_fields['runName'] = custom_runName ?: workflow.runName
    email_fields['success'] = workflow.success
    email_fields['dateComplete'] = workflow.complete
    email_fields['duration'] = workflow.duration
    email_fields['exitStatus'] = workflow.exitStatus
    email_fields['errorMessage'] = (workflow.errorMessage ?: 'None')
    email_fields['errorReport'] = (workflow.errorReport ?: 'None')
    email_fields['commandLine'] = workflow.commandLine
    email_fields['projectDir'] = workflow.projectDir
    email_fields['summary'] = summary
    email_fields['summary']['Date Started'] = workflow.start
    email_fields['summary']['Date Completed'] = workflow.complete
    email_fields['summary']['Pipeline script file path'] = workflow.scriptFile
    email_fields['summary']['Pipeline script hash ID'] = workflow.scriptId
    if (workflow.repository) email_fields['summary']['Pipeline repository Git URL'] = workflow.repository
    if (workflow.commitId) email_fields['summary']['Pipeline repository Git Commit'] = workflow.commitId
    if (workflow.revision) email_fields['summary']['Pipeline Git branch/tag'] = workflow.revision
    if (workflow.container) email_fields['summary']['Docker image'] = workflow.container
    email_fields['summary']['Nextflow Version'] = workflow.nextflow.version
    email_fields['summary']['Nextflow Build'] = workflow.nextflow.build
    email_fields['summary']['Nextflow Compile Timestamp'] = workflow.nextflow.timestamp

    // Check if we are only sending emails on failure
    email_address = params.email
    if (!params.email && params.email_on_fail && !workflow.success) {
        email_address = params.email_on_fail
    }

    // Render the TXT template
    def engine = new groovy.text.GStringTemplateEngine()
    def tf = new File("$baseDir/assets/email_template.txt")
    def txt_template = engine.createTemplate(tf).make(email_fields)
    def email_txt = txt_template.toString()

    // Render the HTML template
    def hf = new File("$baseDir/assets/email_template.html")
    def html_template = engine.createTemplate(hf).make(email_fields)
    def email_html = html_template.toString()

    // Render the sendmail template
    def smail_fields = [ email: email_address, subject: subject, email_txt: email_txt, email_html: email_html, baseDir: "$baseDir" ]
    def sf = new File("$baseDir/assets/sendmail_template.txt")
    def sendmail_template = engine.createTemplate(sf).make(smail_fields)
    def sendmail_html = sendmail_template.toString()

    // Send the HTML e-mail
    if (email_address) {
        try {
            if ( params.plaintext_email ){ throw GroovyException('Send plaintext e-mail, not HTML') }
            // Try to send HTML e-mail using sendmail
            sendMail( to: email_address,
                      subject: subject,
                      body: email_html )
            // [ 'sendmail', '-t' ].execute() << sendmail_html
            log.info "[MARVEL] Sent summary e-mail to $email_address (sendmail)"
        } catch (all) {
            // Catch failures and try with plaintext
            sendMail( to: email_address,
                      subject: subject,
                      body: email_txt )
            // [ 'mail', '-s', subject, email_address ].execute() << email_txt
            log.info "[MARVEL] Sent summary e-mail to $email_address (mail)"
        }
    }

    // Write summary e-mail HTML to a file
    def output_d = new File( "${params.outdir}/pipeline_info/" )
    if (!output_d.exists()) {
        output_d.mkdirs()
    }
    def output_hf = new File( output_d, "pipeline_report.html" )
    output_hf.withWriter { w -> w << email_html }
    def output_tf = new File( output_d, "pipeline_report.txt" )
    output_tf.withWriter { w -> w << email_txt }

    c_reset = params.monochrome_logs ? '' : "\033[0m";
    c_purple = params.monochrome_logs ? '' : "\033[0;35m";
    c_green = params.monochrome_logs ? '' : "\033[0;32m";
    c_red = params.monochrome_logs ? '' : "\033[0;31m";

    if (workflow.stats.ignoredCount > 0 && workflow.success) {
        log.info "${c_purple}Warning, pipeline completed, but with errored process(es) ${c_reset}"
        log.info "${c_red}Number of ignored errored process(es) : ${workflow.stats.ignoredCount} ${c_reset}"
        log.info "${c_green}Number of successfully ran process(es) : ${workflow.stats.succeedCount} ${c_reset}"
    }

    if (workflow.success) {
        log.info "^_^ ${c_purple}[MARVEL]${c_green} Pipeline completed successfully${c_reset}"
    } else {
        checkHostname()
        log.info "Q_Q ${c_purple}[MARVEL]${c_red} Pipeline completed with errors${c_reset}"
    }

}


