"""
fishing.rksparkutil

Simple routines to aid in using Spark

Copyright (c) 2019 Robert J. Klein, Icahn School of Medicine at Mount Sinai
"""

import pyspark.sql.functions as f

from pyspark.sql.types import *

import pandas as pd

"""
Given an array x, removes the first n elements
Called by pandas_udf to enable simple removal of first meta-columns from genotype data, etc.
"""
def slice_array_except_first_n (x, n):
  return(x[n:])

"""
Given an array column of Strings and a split string, splits each element of the array with the split string
and makes  map with the first element as key and the second as val (both strings).
Uses two calls to a python UDF as I don't see a better way to do it
"""
def map_from_array (theArray, theDelim):
  def pull_key_val (x, d, kind):
    retval = []
    index = -1
    if (kind == "key"):
      index = 0
    if (kind == "val"):
      index = 1
    if index == -1:
      raise "Bad input"
    for i in x:
      retval.append(i.split(d)[index])
    return(retval)

  pull_key_udf = f.pandas_udf(lambda x: x.apply(pull_key_val, args=(theDelim,"key")),ArrayType(StringType()))
  pull_val_udf = f.pandas_udf(lambda x: x.apply(pull_key_val, args=(theDelim,"val")),ArrayType(StringType()))

  return(f.map_from_arrays(pull_key_udf(theArray), pull_val_udf(theArray)))

"""
Given two delimited string columns of equal length, returns an array of three Integers representing the # of matches, the # of
mismatches, and the total # of entries
Use the delimited string since pyarrow doesn't allow passing in complex structures
"""
def compare_two_delimited_strings (s1, s2, delim="\t"):

  def compare_strings (x, the_delim):
    # x is array of the two strings
    match = 0
    mismatch = 0
    total = 0
    a = x[0].split(the_delim)
    b = x[1].split(the_delim)
    if (len(a) == len(b)):
      for i in range(0,len(a)):
        total += 1
        if (a[i] == b[i]):
          match += 1
        else:
          mismatch += 1
    
    return (match, mismatch, total)

  compare_strings_udf = f.pandas_udf(lambda x: x.apply(compare_strings, args=(delim)), ArrayType(IntegerType()))

  return(compare_array_udf(f.array(s1, s2)))


"""
parse_array_of_maps
Takes a list of dict as input, and returns a delim-separated String for the 
list of the values sorted by the keys
"""
def parse_array_of_maps (data, delim):
  retval = []
  temp = dict()
  for i in data:
    for j in i.keys():
      temp[j] = i[j]
  k = temp.keys()
  k.sort()
  for i in k:
    retval.append(temp[i])
  return(delim.join(retval))

def withColumnsRenamed(data, orig, renamed):
  retval = data
  for i in range(0,len(orig)):
    retval = retval.withColumnRenamed(orig[i], renamed[i])
  return(retval)

"""
Ordered array from pairs
Given an array of pairs index=value, fills in the values to the appropriate place by index
Returned array starts with index 0 and ends with index max(index)
"""
def ordered_array_from_pairs (x):
  input_map = dict()
  for i in x:
    s = i.split('=')
    input_map[int(s[0])] = s[1]
  retval = [ None for q in range(max(input_map.keys())+1) ]
  for i in input_map:
    retval[i] = input_map[i]
  return(retval)
