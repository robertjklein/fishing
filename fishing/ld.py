"""
fishing.ld module

This module includes routines computing pairwise LD

Copyright (c) 2019 Robert J. Klein, Icahn School of Medicine at Mount Sinai
"""

from pyspark.sql.functions import abs, corr
from pyspark.sql.types import *

def create_ld_plan (first_variant_set, second_variant_set, max_dist):
  """
  Given two sets of (potentially overlapping) variants, creates the cross product of them filtered on the max_dist
  If max_dist < 0, then take all pairs
  Only select the VAR_IDX column when done
  """

  first_set = first_variant_set.withColumnRenamed("VAR_IDX","VAR_IDX1").withColumnRenamed("CHR","CHR1").withColumnRenamed("POS","POS1").withColumnRenamed("filename","filename1")
  second_set = second_variant_set.withColumnRenamed("VAR_IDX","VAR_IDX2").withColumnRenamed("CHR","CHR2").withColumnRenamed("POS","POS2").withColumnRenamed("filename","filename2")

  all_tests = first_set.crossJoin(second_set)
  if max_dist >= 0:
    tests = all_tests.filter(all_tests.CHR1 == all_tests.CHR2).filter(abs(all_tests.POS1 - all_tests.POS2)<=max_dist)
  else:
    tests = all_tests

  return(tests.select("filename1","VAR_IDX1","filename2","VAR_IDX2"))

def add_genotypes_to_ld_plan (ld_plan, gts):
  """
  Given an LD plan and a processed genotype array with GT_ADD, GT_HAPA, and GT_HAPB fields, adds these GT array fields to the plan, named appropriately
  for the first or second SNP in the plan
  """

  # Drop columns we don't need to save space
  dropped_gts = gts.drop('N_ALT','GT_STR')

  with_first_snp = ld_plan.join(dropped_gts, (ld_plan.VAR_IDX1 == dropped_gts.VAR_IDX) & (ld_plan.filename1 == dropped_gts.filename), "inner").drop("filename","VAR_IDX").withColumnRenamed("SAMPLE_IDX","SAMPLE_IDX1").withColumnRenamed("GT_ADD","GT_ADD1").withColumnRenamed("GT_HAPA","GT_HAPA1").withColumnRenamed("GT_HAPB","GT_HAPB1")
  with_second_snp = with_first_snp.join(dropped_gts, (with_first_snp.filename2 == dropped_gts.filename) & (with_first_snp.VAR_IDX2 == dropped_gts.VAR_IDX) & (with_first_snp.SAMPLE_IDX1 == dropped_gts.SAMPLE_IDX) ,"inner").drop("filename","VAR_IDX","SAMPLE_IDX1").withColumnRenamed("GT_ADD","GT_ADD2").withColumnRenamed("GT_HAPA","GT_HAPA2").withColumnRenamed("GT_HAPB","GT_HAPB2")

  return(with_second_snp).repartition(with_second_snp.VAR_IDX1, with_second_snp.filename1, with_second_snp.VAR_IDX2, with_second_snp.filename2)

def add_hap_r2_to_gt_ld_plan (ld_plan):
  """
  Given an LD plan to which haplotypes have been added with add_genotypes_to_ld_plan, and which we assume has been repartitioned by VAR_IDX and filename 1&2, 
  aggregates by variant to compute hap_r2 statistics
  Does so as follows:
  Standard r2 calculation requires 4 numbers: n00, n01, n10, n11
  We can easily get the following: n_1 by summing GT_HAPA2 with GT_HAPB2
    n1_ by summing GT_HAPA1 with GT_HAPB1
  n by taking 2*count
  n11 is sum of (GT_HAPA2*GT_HAPA1)+(GT_HAPB2*GT_HAPB1)
  Then, can easily compute r2 = (n00*n11 - n01*n10)^2/(n1_ * n_1 * n0_ * n_0)
  """
  
  # First, remove NA
  data = ld_plan.select("filename1","VAR_IDX1","filename2","VAR_IDX2","SAMPLE_IDX","GT_HAPA1","GT_HAPB1","GT_HAPA2","GT_HAPB2").dropna()

  # Next, compute column for n11
  data_with_n11 = data.withColumn("GT_HAPAB_11", (data.GT_HAPA1 * data.GT_HAPA2)  + (data.GT_HAPB1 * data.GT_HAPB2))

  # Now aggregate
  x = data_with_n11.groupby(['filename1','VAR_IDX1','filename2','VAR_IDX2']).agg({'GT_HAPA1':'sum', 'GT_HAPA2':'sum','GT_HAPB1':'sum','GT_HAPB2':'sum','GT_HAPAB_11':'sum', 'SAMPLE_IDX':'count'}) \
    .withColumnRenamed('sum(GT_HAPA1)','n1_a').withColumnRenamed('sum(GT_HAPB1)','n1_b').withColumnRenamed('sum(GT_HAPA2)','n_1a').withColumnRenamed('sum(GT_HAPB2)','n_1b') \
    .withColumnRenamed('sum(GT_HAPAB_11)','n11').withColumnRenamed('count(SAMPLE_IDX)','half_n') 

  # Add n1_ and n_1 columns
  x1 = x.withColumn("n1_", x.n1_a + x.n1_b)
  d = x1.withColumn("n_1",x.n_1a+x.n_1b)

  # And do the computation in one fell swoop, as annotated here
  hap_r2_result = d.withColumn('hap_r2',   \
    ((((2.0 * d.half_n) - d.n1_ - d.n_1 + d.n11) * d.n11) - ((d.n1_ - d.n11)*(d.n_1 - d.n11))) *  \
    ((((2.0 * d.half_n) - d.n1_ - d.n_1 + d.n11) * d.n11) - ((d.n1_ - d.n11)*(d.n_1 - d.n11))) /  \
    (d.n1_ * d.n_1 * (2.0*d.half_n - d.n_1) * (2.0*d.half_n - d.n1_)))                         

  return(hap_r2_result.select("filename1","VAR_IDX1","filename2","VAR_IDX2","hap_r2"))

def add_geno_r2_to_gt_ld_plan (ld_plan):
  """
  Given an LD plan to which additive genotypes have been added with add_genotypes_to_ld_plan, 
  aggregates by variant to compute geno_r2 statistics
  r2 calculation comes from Rogers and Huff 2008 -- it's just Pearson's r2
  for the additive genotype
  """
  
  # First, remove NA
  data = ld_plan.select("filename1","VAR_IDX1","filename2","VAR_IDX2","SAMPLE_IDX","GT_ADD1","GT_ADD2").dropna()

  # Now aggregate and compute r
  x = data.groupby(['filename1','VAR_IDX1','filename2','VAR_IDX2']).agg(corr("GT_ADD1","GT_ADD2").alias("r"))

  # Now, square it to get r2
  geno_r2_result = x.withColumn('geno_r2',x.r * x.r)
  
  return(geno_r2_result.select("filename1","VAR_IDX1","filename2","VAR_IDX2","geno_r2"))

  

