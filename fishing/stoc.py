# Simple Tool of Compression (STOC)
# Given a delimited text file, keeps the first N columns, compresses next M 
# singly, and then compresses the rest together.
# Compressed data stored as base64 strings
# This is both a module with callable functions and a command line tool
#
# Robert J. Klein
# October 19, 2018

import base64
try:
  import lzma
  have_lzma = True
except:
  have_lzma = False
import gzip

"""
Function string_encode.  Internal function that takes a string and encodes
it with lzma if available
"""
def string_encode(string_to_encode):
  if (have_lzma == False):
    raise "Need lzma to encode string"
  return((base64.b64encode(lzma.compress(string_to_encode.encode('utf-8'), \
                        lzma.FORMAT_XZ,lzma.CHECK_SHA256) \
                        )).decode('utf-8')
def encode (fullline, the_delim, ncols_raw, ncols_sep_compress):
    splitcols = fullline.rstrip().split(the_delim)
    t = splitcols[:ncols_raw]
    if t[0][:1] == '#':
        string_to_encode = fullline.rstrip()
        for i in range(0,ncols_raw):
            if (i == len(t)):
                t.append('#')
            else:
                t[i] = '#'
        for i in range(0,n_cols_sep_compress):
          t.append('#')
    else:
        for i in range(ncols_raw, ncols_raw + ncols_sep_compress):
          t[i] = string_encode(splitcols[i])
        string_to_encode = the_delim.join(splitcols[(ncols_raw + ncols_sep_compress):])
    t.append(string_encode(string_to_encode))
    return(t)

import argparse
import sys



parser = argparse.ArgumentParser()
parser.add_argument("file",help="The file to compress with STOC")
parser.add_argument("-n",help="number of columns not to compress",type=int,default=1)
parser.add_argument("-d",help="column delimiter",default = "\t")
parser.add_argument("-t",help="number of processing threads",default = '1')

args = parser.parse_args()

n_data_threads = int(args.t)

data_queue = queue.Queue(maxsize = 2 * n_data_threads)
output_queue = queue.PriorityQueue()

input_done = threading.Event()
output_done = threading.Event()



lineno = 0

def output_function ():
    if (n_data_threads == 1):
        (q_line_no, output_string) = output_queue.get()
        print(output_string)
    else:
        last_written = 0
        while (not input_done.is_set()) or (not output_queue.empty()) or (lineno != last_written):
            try:
                (q_line_no, output_string) = output_queue.get()
                if (q_line_no == last_written+1):
                    print(output_string)
                    last_written = q_line_no
                else:
                    output_queue.put((q_line_no, output_string))
            except:
                excepted = 1
        sys.stdout.flush()
        sys.stderr.flush()
        output_done.set()
    return


def data_function ():
    if (n_data_threads == 1):
        (line, delim, cols, lineno) = data_queue.get()
        t = encode_function(line, delim, cols)
        output_queue.put((lineno, delim.join(t)))
    else:
        while (not input_done.is_set()) or (not data_queue.empty()):
            try:
                (line, delim, cols, lineno) = data_queue.get()
                t = encode_function(line, delim, cols)
                output_queue.put((lineno, delim.join(t)))
            except:
                excepted = 1
    return

if args.file[-3:] == '.gz':
    file = gzip.open(args.file,"rt")
else:
    file = open(args.file,"r")

if (n_data_threads > 1):
    data_threads = []
    for i in (range(n_data_threads)):
        data_threads.append(threading.Thread(target=data_function,name="Data" + str(i)))
    for i in data_threads:
        i.start()
    output_thread = threading.Thread(target=output_function,name="Output")
    output_thread.start()

for line in file:
    lineno += 1
    added = 0
    while added == 0:
        try:
            data_queue.put((line, args.d, args.n, lineno))
            added = 1
        except queue.Full: 
            added = 0
    if (n_data_threads == 1):
        data_function()
        output_function()


if (n_data_threads > 1):
    input_done.set()
    output_done.wait()

sys.exit()







