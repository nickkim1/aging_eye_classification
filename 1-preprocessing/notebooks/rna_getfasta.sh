# This is a script that runs bedtools getfasta
# Reference genomes should be provided from the command line 

FASTA=$1
DIR=$2
for file in ../../0-data/1-experimental/bed/${DIR}/*; do 
    # if [ "$file" -eq *.bed ]; then
    filename=$(basename -- "$file")
    filename_wo_extension=${filename%.bed}
    bedtools getfasta -fi $FASTA -bed $file -fo ../../0-data/2-final/0-getfasta/${DIR}/$filename_wo_extension.out.fa -tab
    # fi 
done 