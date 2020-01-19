params.pubdir

process index_fasta {
  tag "$name"
  label 'process_light'
  publishDir "$params.pubdir"

  input:
  tuple val(name), file(fa)

  output:
  tuple file(fa), file("${fa}*")

  script:
  """
  samtools faidx ${fa} > ${fa}.fai
  """
}