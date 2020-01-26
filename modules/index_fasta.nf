params.pubdir

process index_fasta {
  tag "$name"
  label 'process_light'
  publishDir "$params.pubdir"

  input:
  tuple val(name), file(fa)

  output:
  tuple val(name), file(fa), file("${fa}.fai")

  script:
  """
  samtools faidx ${fa} > ${fa}.fai
  """
}