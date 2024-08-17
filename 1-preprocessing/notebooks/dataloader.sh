# this is a script that loads in data 

GTF_LINK=https://ftp.flybase.org/genomes/dmel/current/gtf/dmel-all-r6.56.gtf.gz

wget -P ../../0-data/0-refs/agingeye/ ${GTF_LINK}

# maybe have some sort of a more extnsible format (e.g. loop for downloading in BED files from encode or something)

COUNTS_LINK=...
wget -P ../../0-data/1-experimental/counts/agingeye/ ${COUNTS_LINK}