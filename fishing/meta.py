"""
fishing.meta module

This module includes routines for meta-analysis

Copyright (c) 2020 Robert J. Klein, Icahn School of Medicine at Mount Sinai
"""

import pyspark.sql.functions as f

from scipy.stats import norm

"""
fixed_effects_beta
Given a data frame with grouping columns, groups on the grouping columns and uses Willer and Abecasis meta-analysis
method from METAL to do meta-anlysis using beta and se_beta
"""
def fixed_effects_beta (df, grouping_columns):
  temp = df.withColumn('w_i', 1.0 / (df.SE_BETA * df.SE_BETA))
  temp2 = temp.withColumn('w_i_beta_i', temp.BETA * temp.w_i)
  grouped = temp2.withColumn('studies', f.lit(1)).groupBy(grouping_columns).agg(f.sum('n'), f.sum('w_i'), f.sum('w_i_beta_i'), f.sum('studies')).withColumnRenamed('sum(w_i)', 'sum_w_i').withColumnRenamed('sum(w_i_beta_i)', 'sum_w_i_beta_i')
  final = grouped.withColumn('SE', f.sqrt(1.0 / grouped.sum_w_i)).withColumn('BETA', grouped.sum_w_i_beta_i /grouped.sum_w_i)
  return(final.withColumn('Z', final.BETA/final.SE).select(grouping_columns + ['BETA','SE','Z','sum(studies)','sum(n)']))

"""
fixed_effects_p
Given a data frame with grouping columns, groups on the grouping columns and uses Willer and Abecasis Sample Size based
meta-analysis method from METAL using P, n, and sign of BETA
"""
def fixed_effects_p (df, grouping_columns):
  inverse_normal_udf = f.pandas_udf (lambda x: x.apply(norm.ppf), 'float')
  temp = df.withColumn('inverse_normal', inverse_normal_udf (df.P)).withColumn('sign', df.BETA / f.sqrt(df.BETA*df.BETA))
  temp1 = temp.withColumn('Z_i', f.abs(temp.inverse_normal) * temp.sign).withColumn('w_i', f.sqrt(temp.n))
  temp2 = temp1.withColumn('Z_i_w_i', temp1.Z_i * temp1.w_i).withColumn('w_i_sq', temp1.w_i * temp1.w_i)
  grouped = temp2.withColumn('studies', f.lit(1)).groupBy(grouping_columns).agg(f.sum('n'), f.sum('Z_i_w_i'), f.sum('w_i_sq'), f.sum('studies')).withColumnRenamed('sum(Z_i_w_i)','sum_Z_i_w_i').withColumnRenamed('sum(w_i_sq)', 'sum_w_i_sq')
  final = grouped.withColumn('Z', grouped.sum_Z_i_w_i / f.sqrt(grouped.sum_w_i_sq))
  return(final.select(grouping_columns + ['sum(n)','sum(studies)','Z']))

def add_p_to_z (df):
  normal_udf = f.pandas_udf(lambda x: x.apply(norm.sf), 'float')
  return(df.withColumn('P', 2.0 * normal_udf(f.abs(df.Z))))
