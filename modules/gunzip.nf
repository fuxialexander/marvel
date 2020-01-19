params.pubdir

process gunzip {
  label 'process_light'
  publishDir "$params.pubdir"

  input:
  file gz

  output:
  file "$gz.baseName"

  script:
  """
  gunzip -k --verbose --stdout --force $gz > $gz.baseName
  """
}