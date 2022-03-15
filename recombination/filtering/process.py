#!/bin/python3

import time
import random
import sys
import subprocess
import datetime

version = sys.argv[1]
mat = sys.argv[2]
start_range = sys.argv[3]
end_range = sys.argv[4]
BUCKET_ID = sys.argv[5]
# Remote GCP Storage Bucket location to put end results
RESULTS = sys.argv[6]

def convert(n):
    return str(datetime.timedelta(seconds = n))

def parse_ripples_command(version, MAT, start, end):
    # Expecting recombination.txt and descendents.txt in recombination/filtering to start pipeline
    command = "{} -i {} -n 2 -S {} -E {} -d recombination/filtering/".format(version, MAT, start, end)
    print(command)
    return command

# Start timing test
#start_time = time.time()

# Copy over protobuf from GCP Storage into /data directory inside container
#TODO: Fix destination directory here
subprocess.run(["gsutil", "cp", "gs://{}/{}".format(BUCKET_ID, mat), "/home/usher"])

cmd = parse_ripples_command(version, mat, start_range, end_range)
# Run ripples on current GCP instance
subprocess.run(cmd)

#TODO: Verify correct top level path inside Docker container on remote machine
filtration = ["run_ripples_filtration.sh"]
# Run putative recombinants through post-processing filtration pipeline
subprocess.run(filtration)

# Specify output directory on GCP Storage and copy over list of true-positive recombinants
output_remote = "gs://ripples-testing/ripples_out"
RESULTS = sys.argv[6]
#TODO: again path check
#TODO: You need unique output naming 
local_results = "recombination/filtering/results/"
copy_cmd = ["gsutil", "cp", "-r", local_results, output_remote]
subprocess.run(copy_cmd)

# End timing
#end_time = time.time()
#delta = end_time - start_time
#f = open("/data/time10.txt", "w+")
#f.write(str(convert(delta)))
#f.close()
#copy_time = ["gsutil", "cp", "/data/time10.txt", "gs://ripples-testing/logging/"]
#subprocess.run(copy_time)

