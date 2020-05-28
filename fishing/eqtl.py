"""
fishing.eqtl module

This module includes routines for eqtl calculation

Copyright (c) 2020 Robert J. Klein, Icahn School of Medicine at Mount Sinai
"""


import pandas as pd
import statsmodels.api as sm
from scipy.stats import t, norm

from pyspark.sql.functions import array, pandas_udf, PandasUDFType, lit, covar_samp, var_samp, count, sqrt, abs
from pyspark.sql.types import *

"""
covariate_adjustment

Given a dataframe, grouping columns, an input outcome column name, an output output column name, and a list of 
covariates, adjusts the input outcome for the covariates and returns a dataframe with the residuals added as the 
output column.
"""
def covariate_adjustment (df, grouping_columns, input_outcome_name, output_outcome_name, covariate_list):
  schema = df.withColumn(output_outcome_name, lit(1.0)).schema

  @pandas_udf(schema, PandasUDFType.GROUPED_MAP)
  def adjust (pdf):
    y = pdf[[input_outcome_name]]
    x = pdf[covariate_list]
    rlm_result = sm.RLM(y, x, m=sm.robust.norms.TukeyBiweight())
    results = rlm_result.fit()
    resid = results.resid.to_frame(name=output_outcome_name)
    return(pd.merge(pdf, resid, right_index=True, left_index=True))

  return(df.groupBy(grouping_columns).apply(adjust))


"""
simple_regression
X data frame: genetic data  GT_SAMPLE_NAME|    STUDY|CHR|      POS|REF|ALT|     SNP|          GT_dosage|    TN
Y data frame: expression ('GT_SAMPLE_NAME', 'GENE','STUDY','TN', 'ADJ_EXP')
Given an x data frame and a y data frame, and link data frame, computes the eQTL regressions
Based on https://mathworld.wolfram.com/LeastSquaresFitting.html and Hogg and Tanis 9th edition section 6.5 and 9.6
"""
def simple_regression (x, y, link):
  # First, join everything together

  joined_data = x.join(link, 'SNP', 'inner').join(y, ['GENE', 'GT_SAMPLE_NAME', 'STUDY','TN'], 'inner')

  df = joined_data.groupBy('SNP','GENE','STUDY','TN').agg(var_samp('GT_dosage'), var_samp('ADJ_EXP'), covar_samp('GT_dosage','ADJ_EXP'), count('GT_dosage')).withColumnRenamed('var_samp(GT_dosage)', 'ss_xx').withColumnRenamed('var_samp(ADJ_EXP)','ss_yy').withColumnRenamed('covar_samp(GT_dosage, ADJ_EXP)','ss_xy').withColumnRenamed('count(GT_dosage)','n')
  return(df.select('SNP','GENE','STUDY','TN', (df.ss_xy / df.ss_xx).alias('BETA'), (sqrt((df.ss_yy - (df.ss_xy * df.ss_xy / df.ss_xx))/(df.n - 2.0)) / sqrt(df.ss_xx)).alias('SE_BETA'), 'n').na.drop())

"""
p_from_beta
Given a df with BETA, SE_BETA, and n, adds column P with the p-value from the t-distribution
"""
def p_from_beta (df):
  def compute_p_udf (x):
    return(2*t.sf(x[0]/x[1], x[2] - 2.0))
  compute_p = pandas_udf(lambda x : x.apply(compute_p_udf), 'float')
   
  return(df.withColumn('P', compute_p(array(abs(df.BETA), df.SE_BETA, df.n.cast('float')))))

"""
inverse_rank_transform
Given a data frame, group by grouping_cols and inverse rank transform the value_column
"""
def inverse_rank_transform (df, grouping_cols, value_column):
  schema = df.schema

  @pandas_udf(schema, PandasUDFType.GROUPED_MAP)
  def irt (pdf):
    pdf[value_column] = norm.ppf((pdf[value_column].rank()-0.5)/(pdf.shape[0]))
    return(pdf)

  return(df.groupBy(grouping_cols).apply(irt))

 
