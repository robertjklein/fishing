"""
fishing.prs module

Routine for computing PRS

Copyright (c) 2019 Robert J. Klein, Icahn School of Medicine at Mount Sinai
"""

from fishing.variants import get_vaf

import pyspark.sql.functions as f


"""
compute_prs (vars, gts, prs)
Given a variant table, a genotype table, and a PRS table, computes the PRS and returns
a data frame containing subject identifiers (with filename), PRS, and number of vars in
PRS and number in the var table
"""
def compute_prs (vars, gts, prs):
  # First, extract the variant sites from the var table and merge with PRS
  prs.cache()

  num_prs_snps = prs.count()
  vars_to_use = vars.join(prs, ["CHR","POS","REF","ALT"], "inner").distinct().cache()
  num_vars = vars_to_use.count()

  gt_subset = gts.join(vars_to_use, ["filename","VAR_IDX"], "inner").withColumn("havegt",f.lit(True)).distinct()
  subject_by_var = gts.select("SAMPLE_IDX").distinct().crossJoin(vars_to_use.select("VAR_IDX").distinct()).distinct().cache()
  gt_subset_with_vaf = gt_subset.join(subject_by_var, ["VAR_IDX","SAMPLE_IDX"], "right_outer").join(get_vaf(gt_subset).select("filename","VAR_IDX","VAF"), ["filename","VAR_IDX"], "left_outer").replace(to_replace=["REF","ALT"],value=["-1.0","1.0"],subset="effect_allele").na.fill(False, "havegt").cache()



  with_have_gt = gt_subset_with_vaf
#  print("There are " + str(with_have_gt.select("SAMPLE_IDX").distinct().count()) + " subjects and " + str(with_have_gt.select("VAR_IDX").distinct().count()) + " variants")
#  print(with_have_gt.columns)
# 'VAR_IDX', 'SAMPLE_IDX', 'RAWGT', 'alleles', 'GT_ADD', 'CHR', 'POS', 'REF', 'ALT',
  gt1 = with_have_gt.groupBy('SAMPLE_IDX','VAR_IDX').count().toDF('SAMPLE_IDX','VAR_IDX','c').filter("c>1")
  gt1.join(with_have_gt, ["SAMPLE_IDX","VAR_IDX"], "inner").show()
  with_have_gt.groupBy("havegt").count().show()

  have_gt = with_have_gt.filter(with_have_gt.havegt == True).drop("VAF").cache()
  no_have_gt = with_have_gt.filter(with_have_gt.havegt == False).drop("GT_ADD").withColumnRenamed("VAF","GT_ADD").cache()
#  print("Using " + str(have_gt.count()) + " origs and " + str(no_have_gt.count()) + " imputed from VAF")
  gt_imputed = no_have_gt.unionByName(have_gt).select("filename","SAMPLE_IDX","effect_allele","GT_ADD","BETA")
#  print(gt_imputed.select("SAMPLE_IDX").groupBy("SAMPLE_IDX").count().collect())
  result = gt_imputed.withColumn("score",gt_imputed.BETA*((gt_imputed.effect_allele.cast("float") * (gt_imputed.GT_ADD - 1.0)) + 1.0)).groupBy("filename","SAMPLE_IDX").sum("score").withColumnRenamed("sum(score)","prs_result")
  return(result)

