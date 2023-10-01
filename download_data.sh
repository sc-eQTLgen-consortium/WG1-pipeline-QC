#!/bin/bash

which wget >/dev/null || { echo "wget: command not found"; exit 1; }
which md5sum >/dev/null || { echo "md5sum: command not found"; exit 1; }
which tar >/dev/null || { echo "tar: command not found"; exit 1; }
which gunzip >/dev/null || { echo "gunzip: command not found"; exit 1; }

mkdir -p data
cd data || exit

echo "Downloading eQTLGenImpRef.tar.gz"
wget -q https://www.dropbox.com/s/l60a2r3e4vo78mn/eQTLGenImpRef.tar.gz \
  && md5sum -c - <<<"88e3603933a21712a7403023c1f2e9df  eQTLGenImpRef.tar.gz" \
  && tar -xzf eQTLGenImpRef.tar.gz
rm eQTLGenImpRef.tar.gz

echo "Downloading 1000G.tar.gz"
wget -q https://www.dropbox.com/s/xso2vt3p9h2rh8m/1000G.tar.gz \
  && md5sum -c - <<<"77e1141441b2b443519249d331dab5ae 1000G.tar.gz" \
  && tar -xzf 1000G.tar.gz
rm 1000G.tar.gz

echo "Downloading hg38exonsUCSC.bed"
wget -q https://www.dropbox.com/s/fvd4pl8no3ngg0l/hg38exonsUCSC.bed \
  && md5sum -c - <<<"c323d4fb5da690ca6352a7ee8a14e022 hg38exonsUCSC.bed" \
  || rm hg38exonsUCSC.bed

echo "Downloading GRCh37_to_GRCh38.chain.gz"
wget -q https://ftp.ensembl.org/pub/assembly_mapping/homo_sapiens/GRCh37_to_GRCh38.chain.gz \
  && md5sum -c - <<<"2a862e70fa1afc08f1f3e7a8cc5cda14  GRCh37_to_GRCh38.chain.gz" \
  && gunzip GRCh37_to_GRCh38.chain.gz \
  || rm GRCh37_to_GRCh38.chain.gz

echo "Downloading hg18ToHg38.over.chain.gz"
wget -q http://hgdownload.cse.ucsc.edu/goldenPath/hg18/liftOver/hg18ToHg38.over.chain.gz \
  && md5sum -c - <<<"cb4ba08996a10b0ded40fae00987d1c8 hg18ToHg38.over.chain.gz" \
  && gunzip hg18ToHg38.over.chain.gz \
  || rm hg18ToHg38.over.chain.gz