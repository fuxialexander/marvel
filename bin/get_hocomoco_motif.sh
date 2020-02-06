#!/usr/bin/env bash
# somehow MOODS only perform well with jaspar format converted PWMs
# Download the JASPAR format motifs
wget http://hocomoco11.autosome.ru/final_bundle/hocomoco11/full/HUMAN/mono/HOCOMOCOv11_full_HUMAN_mono_jaspar_format.txt
echo "\n" >> HOCOMOCOv11_full_HUMAN_mono_jaspar_format.txt
# Split the motifs into multiple files
while read line
do
    if [[ ${line:0:1} == '>' ]]
    then
        outfile=${line#>}.pwm
    else
        echo $line | tr ' ' '\t' >> $outfile
    fi
done < HOCOMOCOv11_full_HUMAN_mono_jaspar_format.txt

# Compile an index of motifs
grep ">" HOCOMOCOv11_full_HUMAN_mono_jaspar_format.txt | awk -F">" '{print NR, $2}' | tr ' ' '\t' > motif_index.txt
# Download motif annotation
wget http://hocomoco11.autosome.ru/final_bundle/hocomoco11/full/HUMAN/mono/HOCOMOCOv11_full_annotation_HUMAN_mono.tsv
# Get a dict of motif and corresponding symbol
cut -f1,2 HOCOMOCOv11_full_annotation_HUMAN_mono.tsv | sort -k2,2 | awk -F"[_\t]" '{print $1"\t"$3}' > motif_symbol.txt