"""
fishing.variants module

This module includes routines processing the samples

Copyright (c) 2019 Robert J. Klein, Icahn School of Medicine at Mount Sinai
"""

from pyspark.sql.functions import broadcast as can_broadcast
from pyspark.sql.types import *

def filter_samples_by_name (spark, sample_df, sample_list):
  """
  Given a spark SQL session, alist of sample ids in sample list, returns the subset of rows from sample_df that have an ID matching the SNP name
  """

  # First, create a single column data-frame from the list, find unique elements, and mark it as being broadcastable
  sample_list_df = spark.createDataFrame(sample_list,schema=StringType()).toDF("SELECTSAMPLE")
  can_broadcast(sample_list_df)

  # Now, do the join with the variant_df, and return after dropping SELECTSNP column
  return(sample_df.join(sample_list_df, sample_df.SAMPLE_NAME == sample_list_df.SELECTSAMPLE, 'inner').drop("SELECTSAMPLE"))



  
