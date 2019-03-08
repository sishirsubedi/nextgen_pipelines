#!/bin/bash

start_vep()
{
  /opt/perl/bin/perl /opt/vep/ensembl-tools-release-83/scripts/variant_effect_predictor/variant_effect_predictor.pl \
  -i $1 \
  -o $2 \
  --offline \
  --dir_cache /opt/vep/ensembl-tools-release-83/cache/ \
  --sift p \
  --polyphen p \
  --hgvs \
  --symbol \
  --vcf \
  --pubmed \
  --fasta /home/doc/ref/ref_genome/Homo_sapiens.GRCh37.73.dna_sm.primary_assembly.fa \
  --pick_allele \
  --check_alleles \
  --force_overwrite \
  --gmaf \
  --maf_1kg
}
