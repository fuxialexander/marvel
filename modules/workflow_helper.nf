def hasExtension(it, extension) {
    it.toString().toLowerCase().endsWith(extension.toLowerCase())
}

def cli_banner(){
    // Log colors ANSI codes
    c_reset = params.monochrome_logs ? '' : "\033[0m";
    c_dim = params.monochrome_logs ? '' : "\033[2m";
    c_white = params.monochrome_logs ? '' : "\033[0;37m";

    return """_${c_dim}_______________________________________________${c_reset}_
    ${c_white}_   _    __     ____   _    _   _____    _${c_reset}
    ${c_white}/  /|    / |    /    ) |   /    /    '   /${c_reset}
${c_dim}---${c_reset}${c_white}/| /${c_reset}${c_dim}-${c_reset}${c_white}|${c_reset}${c_dim}---${c_reset}${c_white}/__|${c_reset}${c_dim}---${c_reset}${c_white}/___ /${c_reset}${c_dim}--${c_reset}${c_white}|${c_reset}${c_dim}--${c_reset}${c_white}/${c_reset}${c_dim}----${c_reset}${c_white}/__${c_reset}${c_dim}------${c_reset}${c_white}/${c_reset}${c_dim}----${c_reset}
  ${c_white}/ |/  |  /   |  /    |   | /    /        /${c_reset}
_${c_white}/${c_reset}${c_dim}__${c_reset}${c_white}/${c_reset}${c_dim}___${c_reset}${c_white}|${c_reset}${c_dim}_${c_reset}${c_white}/${c_reset}${c_dim}____${c_reset}${c_white}|${c_reset}${c_dim}_${c_reset}${c_white}/${c_reset}${c_dim}_____${c_reset}${c_white}|${c_reset}${c_dim}___${c_reset}${c_white}|/${c_reset}${c_dim}____${c_reset}${c_white}/____,${c_reset}${c_dim}___${c_reset}${c_white}/____/${c_reset}${c_dim}_${c_reset}

    """.stripIndent()
}

def checkHostname(){
    def c_reset = params.monochrome_logs ? '' : "\033[0m"
    def c_white = params.monochrome_logs ? '' : "\033[0;37m"
    def c_red = params.monochrome_logs ? '' : "\033[1;91m"
    def c_yellow_bold = params.monochrome_logs ? '' : "\033[1;93m"
    if (params.hostnames) {
        def hostname = "hostname".execute().text.trim()
        params.hostnames.each { prof, hnames ->
            hnames.each { hname ->
                if (hostname.contains(hname) && !workflow.profile.contains(prof)) {
                    log.error "====================================================\n" +
                        "  ${c_red}WARNING!${c_reset} You are running with `-profile $workflow.profile`\n" +
                        "  but your machine hostname is ${c_white}'$hostname'${c_reset}\n" +
                        "  ${c_yellow_bold}It's highly recommended that you use `-profile $prof${c_reset}`\n" +
                        "============================================================"
                }
            }
        }
    }
}

def create_workflow_summary(summary) {
    def yaml_file = workDir.resolve('workflow_summary_mqc.yaml')
    yaml_file.text  = """
    id: 'fuxialexander-marvel-summary'
    description: " - this information is collected when the pipeline is started."
    section_name: 'MARVEL Workflow Summary'
    section_href: 'https://github.com/fuxialexander/marvel'
    plot_type: 'html'
    data: |
        <dl class=\"dl-horizontal\">
${summary.collect { k,v -> "            <dt>$k</dt><dd><samp>${v ?: '<span style=\"color:#999999;\">N/A</a>'}</samp></dd>" }.join("\n")}
        </dl>
    """.stripIndent()

    return yaml_file
}

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
