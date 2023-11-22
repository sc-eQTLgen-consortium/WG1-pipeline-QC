#!/bin/bash

which wget >/dev/null || { echo "wget: command not found"; exit 1; }
which md5sum >/dev/null || { echo "md5sum: command not found"; exit 1; }
which tar >/dev/null || { echo "tar: command not found"; exit 1; }
which gunzip >/dev/null || { echo "gunzip: command not found"; exit 1; }

mkdir -p data
cd data || exit

echo "Downloading eQTLGenImpRef.tar.gz"
wget  \
  && md5sum -c - <<<"0c6cddad981e597852eaa679a00516e1  eQTLGenImpRef.tar.gz" \
  && tar -xzf eQTLGenImpRef.tar.gz
rm eQTLGenImpRef.tar.gz

echo "Downloading 1000G.tar.gz"
wget  \
  && md5sum -c - <<<"22fa624102af2bced1a71ff98a6f39a9 1000G.tar.gz" \
  && tar -xzf 1000G.tar.gz
rm 1000G.tar.gz

echo "Downloading hg38exonsUCSC.bed"
wget https://www.dropbox.com/s/fvd4pl8no3ngg0l/hg38exonsUCSC.bed \
  && md5sum -c - <<<"c323d4fb5da690ca6352a7ee8a14e022 hg38exonsUCSC.bed" \
  || rm hg38exonsUCSC.bed

echo "Downloading GRCh37_to_GRCh38.chain.gz"
wget https://ftp.ensembl.org/pub/assembly_mapping/homo_sapiens/GRCh37_to_GRCh38.chain.gz \
  && md5sum -c - <<<"2a862e70fa1afc08f1f3e7a8cc5cda14  GRCh37_to_GRCh38.chain.gz" \
  && gunzip GRCh37_to_GRCh38.chain.gz \
  || rm GRCh37_to_GRCh38.chain.gz

echo "Downloading hg18ToHg38.over.chain.gz"
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg18/liftOver/hg18ToHg38.over.chain.gz \
  && md5sum -c - <<<"cb4ba08996a10b0ded40fae00987d1c8 hg18ToHg38.over.chain.gz" \
  && gunzip hg18ToHg38.over.chain.gz \
  || rm hg18ToHg38.over.chain.gz

echo "Ribosomal_genes.txt"
wget https://raw.githubusercontent.com/sc-eQTLgen-consortium/WG1-pipeline-QC/master/Demultiplexing/Ribosomal_genes.txt \
  && md5sum -c - <<<"1f232f4815abe534768a32e3203600a3  Ribosomal_genes.txt"

echo "Mitochondrial_genes.txt"
wget https://raw.githubusercontent.com/sc-eQTLgen-consortium/WG1-pipeline-QC/master/Demultiplexing/Mitochondrial_genes.txt \
  && md5sum -c - <<<"c5f65ae3bce28b9340615d605b2bb792  Mitochondrial_genes.txt"