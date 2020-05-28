"""
fishing.variants module
adds a
This module includes routines processing the variant table

Copyright (c) 2019 Robert J. Klein, Icahn School of Medicine at Mount Sinai
"""

import pyspark.sql.functions as f
from pyspark.sql.types import *

def filter_variants_by_name (spark, variant_df, snp_list):
  """
  Given a spark SQL session, a list of SNP ids in snp_list, returns the subset of rows from variant_df that have an ID matching the SNP 
  name.  
  """

  # First, create a single column data-frame from the list, find unique elements, and mark it as being broadcastable
  snp_list_df = spark.createDataFrame(snp_list, schema=StringType()).toDF("ID")
  f.broadcast(snp_list_df)

  # Now, do the join with the variant_df, and return after dropping SELECTSNP column
  return(variant_df.join(snp_list_df, "ID", 'inner'))

def standardize_alleles (variants):
  """
  Given a variant data frame containing columns of REF, ALT, and alleles_present,
  adds a column specifiying how to convert alleles to REF/ALT.  The options are
  'ambiguous '(cannot do it -- probably toss), 'reference' -- the alleles match, 
  and 'revcomp' -- take the reverse complement.  Also return 'mismatch' when there
  is a mismatch.
  """

  def determine_type (x):
    """
    Given an array of REF, ALT, and alleles, determines the type of allele change and returns
    a string
    """
    ambig = frozenset(['AT','CG','GC','TA'])
    nts = frozenset(['A','C','G','T'])
    revcomp = {'A':'T', 'C':'G','G':'C', 'T':'A'}
    # Hard code accepting A/C/G/T SNVs only -- not sure how to do build conversion, etc. for indels
    for i in x:
      if i not in nts:
        return("mismatch")
    if (x[0]+x[1]) in ambig:
      return("ambiguous")
    if len(x) == 3:
      if (x[0] == x[2]) | (x[1] == x[2]):
        return('reference')
      else:
        return('revcomp')    
    if len(x) == 4:
      ref = x[0]
      alt = x[1]
      a = x[2]
      b = x[3]
      if ((ref == a) & (alt == b)) | ((ref == b) & (alt == a)):
        return('reference')
      if ((ref == revcomp[a]) & (alt == revcomp[b])) | ((ref == revcomp[b]) & (alt == revcomp[a])):
        return('revcomp')
    return('mismatch')
    
  # Set the UDF for determining type
  determine_type_udf = f.pandas_udf(lambda x: x.apply(determine_type),StringType())

  # First, create a column "allele_array" that is an array of REF, ALT, and the
  # alleles_present array
  with_allele_array = variants.withColumn("allele_array", f.concat(f.array(variants.REF, variants.ALT),variants.alleles_present))

  # Now apply the UDF and return it
  
  return(with_allele_array.withColumn("allele_conversion_type", determine_type_udf(with_allele_array.allele_array)).drop("allele_array"))

"""
Given a GT table, with GT_ADD field, computes variant allele frequency
"""
def get_vaf (gts):
  data = gts.dropna(how='any',subset='GT_ADD').groupBy("filename","VAR_IDX").agg(f.sum('GT_ADD'),f.count('GT_ADD')).withColumnRenamed('sum(GT_ADD)','num').withColumnRenamed('count(GT_ADD)','denom')
  return(data.select(data.filename, data.VAR_IDX, ((data.num * 0.5)/(data.denom)).alias("VAF")))
