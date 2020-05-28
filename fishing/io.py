"""
fishing.io module

This module reads in various genetic file formats to an internal, standard Spark SQL format.
The internal format consists of 4 tables:
  metadata: one-column table of strings to store meta-data for future parsing
  variant: ten-column table consisting of an int INDEX and the first 9 columns from VCF format
  samples: minimum two-column table consisting of index and sample name.  
  gt: three-column table with n_variant x n_samples rows consisting of sample index, variant index, and genotype data string.  Genotype data is in a VCF-like
      format
To help with memory usage, it first reads into a VCF (or VCF-like) format, and then processes to 
the internal format so we can do filtering by variant or subject ID first
Copyright (c) 2019 Robert J. Klein, Icahn School of Medicine at Mount Sinai
"""



import pyspark.sql.functions as f
from pyspark.sql.types import *
from pyspark.storagelevel import StorageLevel

import fishing.rksparkutil as rkutil
from fishing.qc import variant_level_missingness_udf

def read_file (spark, txt_file, *read_on_driver):
  """
  Reads a file into memory without parsing.

  Inputs: a SparkSession, and the name of the file in Spark-compatible format (e.g. hdfs://, etc. ok).  Optional True if it should be
   read on driver node singly to avoid out of memory errors for unsplittable formats (esp. .gz files)
  Returns: A data frame with three columns: filename, lineid, and data
"""
#  if read_on_driver == () ||  read_on_driver == (False):
  read_file = spark.read.text(txt_file).toDF('data').select('data', f.monotonically_increasing_id().alias("lineid"),f.lit(txt_file).alias('filename'))
#  else:
#    if read_on_driver == (True):
#      raise "Cannot read on driver yet"
#    else:
#      raise "Invalid call to read_file with " + str(read_on_driver)
  return(read_file)

def metadata_from_vcf (vcf):
  """
  Given a dataframe consisting of filename and data, returns only those lines where data starts with ##
  
  """
  return(vcf.filter(vcf.data[0:2] == "##"))

def samples_from_vcf (vcf):
  """
  Given a dataframe consisting of filename and data, in vcf-like format, returns a list of the subjects
  with column number

  Samples come from line that starts with one #
  Split the line by "\t" to get an array of column names
  """
  temp = vcf.filter(vcf.data[0:2] == '#C').select('filename',f.posexplode(f.split('data','\t'))).filter('pos > 9')
  return(temp.select((temp.pos - 1).alias('SAMPLE_IDX'), temp.col.alias('SAMPLE_NAME'), 'filename'))

def variants_from_vcf (vcf):
  """
  Given a VCF file in a data frame, extract the first 9 variant columns and give them unique identifiers.  Include genotype columns as an array parsed out
  with a pandas udf
  """

  # Get the main data and put a unique index on each variant
  maindata = vcf.filter(vcf.data.startswith('#') == False)
  splitdata = maindata.select("filename",f.split(f.substring_index('data',"[\t ]+",9),"[\t ]+").alias("split_data"),maindata.lineid.alias("VAR_IDX"))
    
  # Now pull out the columns one at a time, casting non-strings to appropriate type.  Split out INFO and FORMAT here
  variant = splitdata.select("filename","VAR_IDX",\
    f.element_at(splitdata.split_data,1).alias("CHR"),\
    f.element_at(splitdata.split_data,2).cast(IntegerType()).alias("POS"),\
    f.element_at(splitdata.split_data,3).alias("ID"),\
    f.element_at(splitdata.split_data,4).alias("REF"),\
    f.element_at(splitdata.split_data,5).alias("ALT"),\
    f.element_at(splitdata.split_data,6).cast(FloatType()).alias("QUAL"),\
    f.element_at(splitdata.split_data,7).alias("FILTER"),\
    f.split(f.element_at(splitdata.split_data,8), ";").alias("INFO"),\
    f.split(f.element_at(splitdata.split_data,9), ":").alias("FORMAT"))
  return(variant)

# TODO: Change to store data as arrays for each field of interest
def gts_from_vcf (vcf, *variant):

  """
  Given a VCF file, and an optional variant file to filter on, returns the gt matrix
  """
  
  if (len(variant)==0):
     maindata = vcf.filter(vcf.data[0:1] != "#")
  else:
    var = variant[0].withColumnRenamed('filename','filename2')
    maindata = vcf.join(var,(var.filename2 == vcf.filename) & (var.VAR_IDX == vcf.lineid), "inner")

  splitdata = maindata.select(maindata.filename,f.split(maindata.data,"[\t ]+").alias("split_data"),maindata.lineid.alias("VAR_IDX"))
    
  # Explode out the genotype data, but leave out the first nine columns
  gt = splitdata.select('filename','VAR_IDX', f.posexplode(splitdata.split_data)).toDF("filename","VAR_IDX","SAMPLE_IDX","GT").filter("SAMPLE_IDX >= 9")
  return(gt)  


def variants_from_tped (tped):
  """
  Given a tped file in a data frame, extracts the variant information to match the VCF format above with applicable columns (CHR, ID, POS), along with MAP and
  an array of available alleles.
  Uses pandas UDFs to convert the splitdata array to get the non-'0' alleles present
  """

  # Define the UDFs locally here.  Since I'm working on arrays, the best way I found to do this is to use the Pandas.Series apply function
  def pandas_get_alleles_present (x):
    return(list(frozenset(x[4:]) - frozenset(['0'])))
  tped_get_alleles_present = f.pandas_udf(lambda x: x.apply(pandas_get_alleles_present), ArrayType(StringType()))

  # Compute missingness directly from the allele array
  def pandas_get_missing (x):
    count = 0
    for i in range(4, len(x), 2):
      if (x[i] == '0') | (x[i+1] == '0'):
        count += 1
    return(count)
  tped_get_missing = f.pandas_udf(lambda x: x.apply(pandas_get_missing), IntegerType())

  # Split the tped file
  splitdata = tped.select("filename",f.split(tped.data,"[\t ]+").alias("split_data"),tped.lineid.alias("VAR_IDX"))

  # Pull out the first four columns with appropriate casts, and get frac_missing and alleles_present from the UDFs above
  with_alleles = splitdata.select("filename","VAR_IDX", \
    f.element_at(splitdata.split_data,1).alias("CHR"), \
    f.element_at(splitdata.split_data,2).alias("ID"), \
    f.element_at(splitdata.split_data,3).cast(FloatType()).alias("MAP"), \
    f.element_at(splitdata.split_data,4).cast(IntegerType()).alias("POS"), \
    ((f.size(splitdata.split_data) - 4) / 2).alias("n_samples"), \
    tped_get_alleles_present("split_data").alias("alleles_present"), \
    tped_get_missing("split_data").cast(FloatType()).alias("missingcnt"))

  return(with_alleles)


def gts_from_tped (tped):
  """
  Given a tped file in a data frame, along with optional REF and ALT alleles and a boolean on whether to take the reverse complement, eextracts genotype data and returns it 
  either as a VCF genotype (if we have REF and ALT) or as a pseudo-VCF genotype (the actual two alleles separated by a /
  """

  def pandas_vcf (x, full_vcf):
    revcomp = {'A':'T', 'C':'G','G':'C', 'T':'A'}
    retval = []
    if full_vcf:
      ref = x[0]
      alt = x[1]
      revcomp_bool = x[2]
      data = x[7:]
    else:
      data = x[4:]
      revcomp_bool = 'false'
    for i in range(0, len(data), 2):
      if data[i] == '0' or data[i+1] == '0':
        retval.append('./.')
      else:
        a = data[i]
        b = data[i+1]
        if full_vcf & (revcomp_bool == 'true'):
          a = revcomp[a]
          b = revcomp[b]
        if full_vcf:
          if a == ref:
            a = '0'
          else:
            if a == alt:
              a = '1'
          if b == ref:
            b = '0'
          else:
            if b == alt:
              b = '1'
        retval.append(a + '/' + b)
    return(retval)
  tped_full_vcf = f.pandas_udf(lambda x: x.apply(pandas_vcf, args=(True,)), ArrayType(StringType()))
  tped_pseudo_vcf = f.pandas_udf(lambda x: x.apply(pandas_vcf, args=(False,)), ArrayType(StringType()))

  # Split the tped file
  colnames = frozenset(tped.columns)
  if ('REF' in colnames) & ('ALT' in colnames) & ('revcomp' in colnames):
    full_vcf = True
  else:
    full_vcf = False
    
  splitdata = tped.select("*", f.split(tped.data,"[\t ]+").alias("split_data")).drop("data")
  if full_vcf:
    splitdata1 = splitdata.withColumn("revcomp_str", splitdata.revcomp.cast(StringType()))
    splitdata = splitdata1.withColumn("new_split_data", f.concat(f.array(splitdata1.REF, splitdata1.ALT, splitdata1.revcomp_str), splitdata1.split_data)).drop("split_data").withColumnRenamed("new_split_data","split_data")
    vcf_array = splitdata.select("filename","VAR_IDX", tped_full_vcf("split_data").alias("vcf")).drop("split_data")
  else:
    vcf_array = splitdata.select("filename", "VAR_IDX", tped_pseudo_vcf("split_data").alias("vcf")).drop("split_data")
  
  # Posexplode the vcf array 
  gt = vcf_array.select("filename", "VAR_IDX", f.posexplode(vcf_array.vcf)).withColumnRenamed("pos","SAMPLE_IDX").withColumnRenamed("col","GT").withColumn("truevcf", f.lit(full_vcf))
  gt_with_missing_status = gt.select("*",gt.GT.contains(".").alias('missing'))
 
  return(gt_with_missing_status)

  
"""
variants_from_impute
Given an impute2 file that is read in, extracts chromosome from column name, the first five columns, 
and returns variant table including these plus the full fileline
"""
def variants_from_impute (infile):
  
  # First map the filename to CHR
  chrs = dict()
  for filenamerow in infile.select("filename").distinct().collect():
    s = filenamerow.filename
    i = s.find("chr")
    if (i>0):
      start = i+3
      stop = (i+3)+s[start:].find(".")
      cur_chr = s[start:stop]
      if cur_chr == 'X':
        cur_chr = '23'
    else:
      cur_chr = 'ND'
    chrs[s] = cur_chr
  
      
  
  # Get the main data and put a unique index on each variant.  Add in the CHR here.
  maindata = infile.filter(infile.data[0:1] != "#")
  splitdata = maindata.select(maindata.filename, \
                              maindata.lineid.alias("VAR_IDX"), \
                              maindata.data, \
                              f.split(maindata.data,"[\t ]+").alias("split_data"), \
                              maindata.filename.alias("CHR")).replace(chrs, subset="CHR")
                              
  # Now pull out the first five columns one at a time, casting non-strings to appropriate type.  
  variant = splitdata.select("filename","VAR_IDX","data","CHR",\
    f.element_at(splitdata.split_data,1).alias("COL1ID"),\
    f.element_at(splitdata.split_data,2).alias("RAWID"),\
    f.element_at(splitdata.split_data,3).cast(IntegerType()).alias("POS"),\
    f.element_at(splitdata.split_data,4).alias("ALLELE1"),\
    f.element_at(splitdata.split_data,5).alias("ALLELE2"))

  # Next, get the rsID, if present in the RAWID column
  variant2 = variant.select("*", f.split(variant.RAWID, ":").alias("split_id"))
  variant3 = variant2.select("*", f.element_at(variant2.split_id, 1).alias("EXTRACTID"))
  
  return(variant3)
  
"""
gts_from_impute
Given an input file data frame, returns the genotype probabilities for each
subject in columns P11, P12, P22
"""
def gts_from_impute (infile):
  # Get the main data and put a unique index on each variant
  maindata = infile.filter(infile.data[0:1] != "#")
  splitdata = maindata.select("filename",f.split(maindata.data,"[\t ]+").alias("split_data"),maindata.lineid.alias("VAR_IDX"))

  gtdata1 = splitdata.select("filename", "VAR_IDX", f.posexplode(splitdata.split_data)).toDF("filename","VAR_IDX","COLUMN_IDX","GTPROB").filter("COLUMN_IDX > 4")
  # Now, get subject ID and which GT
  gtdata2 = gtdata1.select("filename", "VAR_IDX", "GTPROB", "COLUMN_IDX", f.floor((gtdata1.COLUMN_IDX - 5) / 3).alias("SAMPLE_IDX"), ((gtdata1.COLUMN_IDX - 5) % 3).cast(StringType()).alias("GT_IDX"))
  gtdata3 = rkutil.withColumnsRenamed(gtdata2.groupBy("filename","VAR_IDX","SAMPLE_IDX").pivot("GT_IDX",["0","1","2"]).agg(f.collect_list("GTPROB")), ["0","1","2"],["c0","c1","c2"])
  gtdata4 = gtdata3.select("filename","VAR_IDX","SAMPLE_IDX", f.element_at(gtdata3.c0, 1).cast(FloatType()).alias("P11"), f.element_at(gtdata3.c1, 1).cast(FloatType()).alias("P12"), f.element_at(gtdata3.c2, 1).cast(FloatType()).alias("P22"))
  return(gtdata4)
