"""
fishing.gts module

This module includes routines processing the GT data

Copyright (c) 2019 Robert J. Klein, Icahn School of Medicine at Mount Sinai
"""

from pyspark.sql.functions import size, split, element_at, array_position, posexplode, lit, array
from pyspark.sql.types import *

def process_GT_field (variants, gts):
  """
  Given a variant table to get the "FORMAT" field and a gts file already filtered by subject,
  Extracts GT information from the gts field into the following fields:
    GT_HAPA -- 0/1 for haplotype 1, ND if no data
    GT_HAPB -- same for Hap 2
    GT_STR -- simple string of the GT field
    GT_ADD -- additive dosage (0,1,2) of the ALT allele, or ND for missing (includes case of multi-allelic)
  
    Returns GT df with these added fields
  """

  # First, add a column containing N_ALT alleles to the variants data frame, while removing all but VAR_IDX and FORMAT
  var_with_nalt = variants.select("filename","VAR_IDX","FORMAT", size(split(variants.ALT,",")).alias("N_ALT"))

  # Next, add GT_pos string with position of GT in format
  var_with_nalt_gtpos = var_with_nalt.select("filename","VAR_IDX","N_ALT",array_position(split(variants.FORMAT, ":"),"GT").alias("GT_pos"))
  var_with_nalt_gtpos1 = var_with_nalt_gtpos.filter(var_with_nalt_gtpos.GT_pos > 0)
 
  # Now, merge with gt data, and select fields of interest
  gts_varann_join = gts.join(var_with_nalt_gtpos1, ["VAR_IDX","filename"], "inner")
  
  gts_varann = gts_varann_join. \
    select("filename","VAR_IDX","SAMPLE_IDX","N_ALT","GT_pos","GT").repartition("VAR_IDX","filename")
 
  # Split GT data by :
  gts_varann_split = gts_varann.select("*",split(gts_varann.GT,":").alias("splitdata"))
  
  # And select the column with the GT field
  gts_varann_split_explode = gts_varann_split.select("*",posexplode(gts_varann_split.splitdata))
  gts_varann_gtonly = gts_varann_split_explode.filter((gts_varann_split_explode.pos+1) == gts_varann_split_explode.GT_pos).select("filename","VAR_IDX","SAMPLE_IDX","N_ALT","col").withColumnRenamed("col","GT_STR")

  # Add missing data for default for GT_ADD, GT_HAPA, and GT_HAPB
  gts_with_all_columns = gts_varann_gtonly.select("*", lit(None).alias('GT_ADD').cast(IntegerType()), lit(None).alias('GT_HAPA').cast(IntegerType()), lit(None).alias('GT_HAPB').cast(IntegerType()))

  # Now split into N_ALT == 1 and N_ALT > 1
  one_alt = gts_with_all_columns.filter(gts_with_all_columns.N_ALT == 1)
  poly_alt = gts_with_all_columns.exceptAll(one_alt)

  # Now, compute additive GT for the one_alt group
  one_alt_add_gt = one_alt.drop('GT_ADD').withColumn("GT_ADD", one_alt.GT_STR.substr(1,1).cast(ByteType()) + one_alt.GT_STR.substr(3,1).cast(ByteType()))

  # Subset to valid haplotype or not
  one_alt_for_hap = one_alt_add_gt.filter(one_alt_add_gt.GT_STR.substr(2,1) == '|')
  one_alt_not_hap = one_alt_add_gt.exceptAll(one_alt_for_hap)

  one_alt_with_haps = one_alt_for_hap.drop('GT_HAPA','GT_HAPB').withColumn("GT_HAPA",one_alt_for_hap.GT_STR.substr(1,1).cast(ByteType())).withColumn("GT_HAPB",one_alt_for_hap.GT_STR.substr(3,1).cast(ByteType()))

  retval = one_alt_with_haps.unionByName(one_alt_not_hap).unionByName(poly_alt)

  return(retval)

"""
simple_additive_gt
Given a dataframe with "RAWGT" column consisting of plink-like genotype (two 
alleles separated by a space), a REF column, and an ALT column, creates the
additive GT GT_ADD with values 0,1,2
"""
def simple_additive_gt (data):
  df = data.select("*", element_at(split(data.RAWGT, "\s+"), 1).alias("a1"), element_at(split(data.RAWGT, "\s+"), 2).alias("a2"))
  df = df.select("*", (df.a1 == df.REF).alias("a1R"), (df.a2 == df.REF).alias("a2R"), (df.a1 == df.ALT).alias("a1A"), (df.a2 == df.ALT).alias("a2A"))

  # Now, filter out those for which a1 isn't REF or ALT or a2 isn't REF or ALT
  df = df.drop("a1","a2").withColumn("a1", df.a1A | df.a1R).withColumn("a2",df.a2A | df.a2R)

  df.filter((df.a1 == False) | (df.a2 == False)).select("VAR_IDX","RAWGT","REF","ALT").distinct().show()

  # And make the additive
  df = df.filter(df.a1 == True).filter(df.a2== True).filter("RAWGT != '0 0'").withColumn("GT_ADD", df.a1A.cast(IntegerType()) + df.a2A.cast(IntegerType())).drop("a1A","a2A","a1R","a2R","a1","a2","REF","ALT")
  return(df)

  
