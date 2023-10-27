#!/bin/bash

which wget >/dev/null || { echo "wget: command not found"; exit 1; }
which md5sum >/dev/null || { echo "md5sum: command not found"; exit 1; }
which tar >/dev/null || { echo "tar: command not found"; exit 1; }

mkdir -p TestData/Imputation
mkdir -p TestData/Demultiplexing
cd TestData/Imputation || exit

echo "Downloading ImputationTestDataset_plink.tar.gz"
wget https://www.dropbox.com/s/uy9828g1r1jt5xy/ImputationTestDataset_plink.tar.gz \
  && md5sum -c - <<<"9f40b159c2bbad789ba4a73c39147aaa  ImputationTestDataset_plink.tar.gz" \
  && tar -xzf ImputationTestDataset_plink.tar.gz
rm ImputationTestDataset_plink.tar.gz

cd ../Demultiplexing || exit

echo "Downloading TestData4PipelineSmall.tar.gz"
wget https://www.dropbox.com/s/m8u61jn4i1mcktp/TestData4PipelineSmall.tar.gz \
  && md5sum -c - <<<"d2d5f1c2d405588de8e88b2419eebffa  TestData4PipelineSmall.tar.gz" \
  && tar -xzf TestData4PipelineSmall.tar.gz
rm TestData4PipelineSmall.tar.gz

#echo "Downloading TestData4PipelineFull.tar.gz"
#wget https://www.dropbox.com/s/3oujqq98y400rzz/TestData4PipelineFull.tar.gz \
#  && md5sum -c - <<<"137afb046237d6b205023b1205287d78  TestData4PipelineFull.tar.gz" \
#  && tar -xzf TestData4PipelineSmall.tar.gz
#rm TestData4PipelineSmall.tar.gz