"""
fishing.qc module

QC routines 

Copyright (c) 2019 Robert J. Klein, Icahn School of Medicine at Mount Sinai
"""

from pyspark.sql.functions import pandas_udf, array
from pyspark.sql.types import FloatType, IntegerType

"""
Given a data frame witha  column called "missing", a Bool for missigness, returns missingness on a per-subject basis
"""
def subject_level_missingness (data):

  totals = data.groupBy('SAMPLE_IDX').count().withColumnRenamed('count','total')
  missing = data.filter(data.missing==True).groupBy('SAMPLE_IDX').count().withColumnRenamed('count','samp_missing').withColumnRenamed('SAMPLE_IDX','SAMPLE_IDX2')
  j = totals.join(missing, totals.SAMPLE_IDX == missing.SAMPLE_IDX2, "left_outer").drop('SAMPLE_IDX2')
  retval = j.select('SAMPLE_IDX','samp_missing','total', (j.samp_missing / j.total).alias("frac")).fillna(0,['samp_missing',"var_missing","frac"])
  return(retval)

"""
Given a data frame witha  column called "missing", a Bool for missigness, returns missingness on a per-variant basis
"""
def var_level_missingness (data):
  totals = data.groupBy('filename','VAR_IDX').count().withColumnRenamed('count','total')
  missing = data.filter(data.missing==True).groupBy('filename','VAR_IDX').count().withColumnRenamed('count','var_missing')
  j = totals.join(missing, ['filename','VAR_IDX'], "left_outer")
  retval = j.select('filename', 'VAR_IDX','total', (j.var_missing / j.total).alias("frac")).fillna(0,["var_missing","frac"])
  return(retval)

"""
Given a split df, a column name, a number of cols to skip, a stride, and a missing function, computes the fraction missing
on that row
"""
def get_missingness_fraction (col):
  n_missing = 0
  for i in col:
    if i[0] == '.':
      n_missing += 1
  tot = len(col)
  return(n_missing * 1.0 / tot)
  
variant_level_missingness_udf = pandas_udf (lambda x: x.apply(get_missingness_fraction), FloatType())

"""
Given a genotype data frame, compute HWE for each autosomal SNP.
Returns a df with SNP identifiers and chi2 statistic for each SNP
tested.  I don't compute p here because for simple filtering I can 
just compute the chi2 cutoff rather than compute p-value for each SNP
Assumes all genotypes passed in are from variants that are at a diploid site in the individual.
Remove X chromosome from males and Y chromosome in calling code, not here
"""
def get_hwe(data):

  def compute_hwe (obs):
    if len(obs) != 3:
      return(-1.0)
    n = obs[0] + obs[1] + obs[2]
    p = (obs[0] + (0.5 * obs[1]))/n
    chi = 0.0

    # If monomorph, return 0.0 -- 
    if ((p == 0.0) or (p == 1.0)):
      return(0.0)
    
    # AA -- e = p^2
    e = p * p * n
    chi += (((obs[0]-e)**2)/e) 

    # AB == 2p(1-p)
    e = 2.0 *p * (1.0 - p) * n
    chi += (((obs[1] - e)**2)/e)

    #BB -- (1-p)^2
    e = ((1.0 - p)** 2) * n
    chi += (((obs[2] - e)**2)/e)

    return(chi)
  
  compute_hwe_udf = pandas_udf(lambda x: x.apply(compute_hwe), FloatType())

  # First, perform substitutions on the GT data frame to account for potential phased calls as well as 1/0 instead
  # of 0/1
  replaced_data = data.replace({'0|0':'0/0', '0|1': '0/1', '1|0':'0/1', '1/0':'0/1', '1|1':'1/1'}, subset='GT')
  
  # Remove all rows that aren't 0/1 0/1 or 1/1
  filtered_data = replaced_data.filter((replaced_data.GT == '0/0') | (replaced_data.GT == '0/1') | (replaced_data.GT == '1/1'))
  
  # Now get the counts for each of the three genotypes via a pivot
  counts = filtered_data.groupBy("filename","VAR_IDX").pivot("GT", ["0/0","0/1","1/1"]).count()

  # And get the chi2 via a call to the pandas UDF
  return(counts.select("filename","VAR_IDX", "0/0", "0/1", "1/1", compute_hwe_udf(array('0/0','0/1','1/1')).alias("hwe_chi2")))
  

  
