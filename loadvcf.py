"""
loadvcf.py

Given a VCF file, loads it to a parquet store in the same place after parsing
"""

import sys


from pyspark.sql import SparkSession
spark = SparkSession.builder.appName("LoadVcf").config("spark.sql.shuffle.partitions","50").config("spark.rdd.compress","true").getOrCreate()

from pyspark import StorageLevel

infile=sys.argv[1]
outdir=sys.argv[2]

from fishing.io import *

sparkfile = read_file(spark, infile)

print("RJK: Doing metadata")
metadata = metadata_from_vcf(sparkfile)
metadata.write.parquet(outdir + "/meta",mode="overwrite", compression="gzip")
print("RJK: Meta data done.  Doing samples")
samples = samples_from_vcf(spark, sparkfile)
samples.write.parquet(outdir + "/subjects",mode="overwrite", compression="gzip")
print("RJK: Samples done.  Doing variants")
variants = variants_from_vcf(sparkfile).drop("FULLGTDATA")
variants.write.parquet(outdir + "/vars",mode="overwrite", compression="gzip")
print ("RJK: Variants done.  Doing gts")
gts = gts_from_vcf(sparkfile)
gts.write.parquet(outdir + "/gts",mode="overwrite", compression="gzip")
print ("Done")
