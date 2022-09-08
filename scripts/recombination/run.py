#!/bin/python3
#
# Launch script to run parallel ripples jobs on GCP
import subprocess
from multiprocessing import Process
from google.cloud import storage
import sys
import time
import os
import os.path
import pathlib
import re
import json
import time
import timeit
import datetime
from datetime import timedelta
import yaml


def get_config():
    with open('ripples.yaml') as f:
      data = yaml.load(f, Loader=yaml.FullLoader)
      return data 

def auth():
    cmd1 = ["gcloud", "config", "set", "project", project_id]
    cmd2 = ["gcloud", "auth", "activate-service-account", "--key-file", key_file]
    print("Setting up GCP access through gcloud.")
    subprocess.run(cmd1)
    subprocess.run(cmd2)

def get_partitions(long_branches, instances):
    partitions = []
    per_instance = long_branches // instances
    k = 0
    for i in range(1, instances+1):
        # Last partition gets extra 
        if i == instances:
            partitions.append((k, long_branches))
            break
        partitions.append((k, k + per_instance))
        k += per_instance + 1
    return partitions

def create_bucket_folder(project_id, bucket_id, path_from_bucket_root):
    # Check to make sure folder ends with backslash, add if not
    if not path_from_bucket_root.endswith("/"):
        path_from_bucket_root += "/"
    client = storage.Client(project=project_id)
    bucket = client.get_bucket(bucket_id)
    blob = bucket.blob(path_from_bucket_root)
    result = blob.upload_from_string('')

def convert(n):
    return str(datetime.timedelta(seconds = n))

def parse_command(mat, start, end, out, logging):
    command = "python3 process.py {} {} {} {} {} {} {} {} {} {} {}".format(version,
            mat, start, end, out, bucket_id, results, reference, date, logging, num_descendants)
    return command

def run_chronumental(newick, metadata):
    reference_node = "DP0803|LC571037.1|2020-02-17"
    steps = "100"
    log_file = open("chronumental_stdout.log", "w")
    command = ["chronumental", "--tree", newick, "--dates", metadata, "--reference_node", reference_node, "--steps", steps]
    p = subprocess.Popen(command, stdout=log_file)
    return p

def newick_from_mat(mat, out_newick):
    command = ["matUtils", "extract", "-i", mat, "-t", out_newick]
    p = subprocess.run(command)
    return p

# Generate run command for final recombination list and ranking
def post_process(mat, filtration_results_file, chron_dates_file, date, recomb_output_file):
    cmd = ["post_filtration", "-i", mat, "-f", filtration_results_file, "-c", chron_dates_file, "-d", date, "-r", recomb_output_file]
    return cmd

def gcloud_run(command, log):
    cmd = ["gcloud", "beta", "lifesciences", "pipelines", "run",
        "--location", "us-central1",
        "--regions", "us-central1",
        "--machine-type", machine_type, 
        "--boot-disk-size", boot_disk_size,
        "--logging", "gs://{}/{}".format(bucket_id, log),
        "--docker-image", docker_image,
        "--command-line", command,
        #"--outputs", results,
        "--format=json"]

    print(" ".join(cmd))
    out = subprocess.check_output(cmd)
    result = json.loads(out)
    name = result['name']
    id = re.search("operations\/(\d+)", name).groups()[0]
    return {'operation_id': id, 'result': result}

def gcloud_describe(operation_id, key_file):
     # Refresh service account credentials
     auth = ["gcloud", "auth", "activate-service-account", "--quiet", "--key-file", key_file]
     print("Refreshing Service Account credentials.")
     subprocess.run(auth)
     print("Checking job status")
     cmd = ["gcloud", "beta", "lifesciences",
            "operations", "describe", "--format=json", operation_id]
     out = subprocess.check_output(cmd)
     result = json.loads(out)
     done = 'done' in result and result['done'] == True
     return {'done': done, 'result': result}

# Configs and credentials from yaml
config = get_config()

# Authenticate GCP account
bucket_id = config["bucket_id"]
project_id = config["project_id"]  
key_file = config["key_file"]

# Activate credentials for GCP Console and gcloud util
auth()

# Set ripples job config options
#docker_image = "mrkylesmith/ripples_pipeline:latest"
#TODO: Using dev image at the moment
docker_image = "mrkylesmith/ripples_pipeline_dev:latest"
boot_disk_size = str(config["boot_disk_size"])
instances = config["instances"] # Number of remote machines to parallelize ripples across 
machine_type = config["machine_type"]
logging = config["logging"]
version = config["version"]
mat = config["mat"]
newick = config["newick"]
metadata = config["metadata"]
date = config["date"]
reference = config["reference"]
num_descendants = config["num_descendants"]

# Remote location in GCP Storge Bucket to copy filtered recombinants
results = "gs://{}/{}".format(bucket_id, config["results"])

# Check logging file created on GCP bucket
if not logging.endswith("/"):
    logging += "/"

# Create remote logging folder
create_bucket_folder(project_id, bucket_id, logging)
print("Created empty GCP storage bucket folder for logging: {}".format(config["logging"]))

# Create remote results folder
create_bucket_folder(project_id, bucket_id, config["results"])
print("Created empty GCP storage bucket folder for results: {}".format(config["results"]))

current = str(os.getcwd())

# Copy over protobuf and metadata file from GCP storage bucket to local container
if not os.path.isfile("{}/{}".format(current,mat)):
  print("Copying input MAT: {} from GCP Storage bucket into local directory in container.".format(mat))
  subprocess.run(["gsutil", "cp", "gs://{}/{}".format(bucket_id, mat), current])
else:
    print("Input MAT found in local directory.")

if not os.path.isfile("{}/{}".format(current,metadata)):
    print("Copying input MAT metadata: {} from GCP Storage bucket into local directory in container.".format(metadata))
    subprocess.run(["gsutil", "cp", "gs://{}/{}".format(bucket_id, metadata), current])
else:
    print("Input MAT metadata found in local directory.")

# Generate newick tree input for Chronumental, if doesn't exist already
if not os.path.isfile("{}/{}".format(current,newick)):
    print("Generating newick file from MAT using matUtils extract")
    newick_from_mat(mat, newick)

# Launch Chronumental job locally
p1 = Process(target=run_chronumental, args=(newick, metadata)
p1.start()

print("Finding the number of long branches.")
if num_descendants == None:
  init = "ripplesInit -i {}".format(mat)
elif isinstance(num_descendants, int):
  init = "ripplesInit -i {} -n {}".format(mat, num_descendants)
else:
    print("Check ripples.yaml file configuration for num_descendants. Provide int number of descendents to consider or leave field blank to use default value 2.")
    exit(1)

# Run ripples init scripts 
try:
  # Get number of long branches to search
  long_branches = int(subprocess.getoutput(init))
except:
  print("Empty tree input or error with input metadata. Check 'ripples.yaml' config file")
  exit(1)

print("Found {} long branches in {} to split across {} instances (number of instances set in 'ripples.yaml').".format(long_branches, mat, instances))

# Num of long branches to search per instance
branches_per_instance = long_branches//instances        

partitions = get_partitions(long_branches, instances)
print("long branches: {}, instances: {}, branches_per_instance: {}".format(long_branches, instances, branches_per_instance))
print("partitions: {}".format(partitions))


processes = []
completed = []
for partition in partitions:

    start_range = str(partition[0])
    end_range = str(partition[1])
    out = "{}_{}".format(start_range, end_range)
    log = logging + out + ".log"

    # The following command gets executed on remote machine: 
    # python3 process.py <version> <mat> <start> <end> <out> <bucket_id> <results> <reference> <date> <logging> <num_descendants>
    command = parse_command(mat, start_range, end_range, out, logging)

    info = gcloud_run(command, log)
    processes.append({'partition': partition, 'operation_id': info['operation_id']})
    completed.append(False)

# Start total runtime for all instances running in parallel 
start_time = timeit.default_timer()

while not all(completed):
   i = 0
   for process in processes:
     partition = process['partition']
     operation_id = process['operation_id']
     info = gcloud_describe(operation_id, key_file)
     done = info['done']
     if done:
       completed[i] = True
     print("partition: {}, operation_id: {}, done: {}".format(partition, operation_id, done))
     i+=1
   # Query active GCP jobs every 2 mins
   time.sleep(120)

print("All instance jobs have finished.  Aggregating results from remote machines.")

# All job have completed, log total runtime, copy to GCP logging directory
stop_time = timeit.default_timer()
runtime_log = open("aggregate_runtime.log", "w")
runtime_log.write("Timing for recombination detection tree date:{}{}{}".format('\t', date, '\n'))
runtime_log.write("Total runtime searching {} tree with {} long branches:{}{}  (Hours:Minutes:Seconds){}".format(date, long_branches, '\t', str(timedelta(seconds=stop_time - start_time)), '\n'))
runtime_log.close()
subprocess.run(["gsutil", "cp", "aggregate_runtime.log", "gs://{}/{}".format(bucket_id, logging)])

current = str(os.getcwd())

local_results = current + "/{}".format(config["results"])
# Create local output directory 
subprocess.run(["mkdir", "-p", local_results])

# Create local temporary directory to copy remote results and aggregate 
temp = current + "/merge_results/"
subprocess.run(["mkdir", "-p", temp])

# Make sure temp directory was created correctly
if (os.path.isdir(temp) == False):
    print("Local results directory not created. Check naming error")
    raise FileNotFoundError

# Copy over all subdirectories from GCP Bucket to local temp directory
remote_results = "gs://{}/{}/*".format(bucket_id, config["results"])
subprocess.run(["gsutil", "cp", "-r", remote_results, temp])

# File to aggregate all detected recombination events
recombinants = open(local_results + "/recombinants_{}.txt".format(date), "w")
unfiltered_recombinants = open(local_results + "/unfiltered_recombinants_{}.txt".format(date), "w")
descendants = open(local_results + "/descendants_{}.txt".format(date), "w")

# Aggregate the results from all remote machines and combine into one final file in 'results/' dir
for directory in os.listdir(temp):
    subdir = temp + directory 
    print("SUBDIR: ", subdir)
    # Guarantee to be only 3 files in each subdirectory
    files = os.listdir(subdir)
    for file in files:
      # Skip over recombination.tsv and descendents.tsv files for aggregating recombination events
      if "filtered_recombinants.txt" in file:
        f1 = open(subdir + "/" +  file, "r")
        for line in f1:
          # One detected recombinant per line, aggregate all lines in each file
          recombinants.write(line)
        f1.close()
      if "recombination.tsv" in file:
        f2 = open(subdir + "/" +  file, "r")
        for line in f2:
          # One detected recombinant per line, aggregate all lines in each file
          unfiltered_recombinants.write(line)
        f2.close()
      if "descendants.tsv" in file:
        f3 = open(subdir + "/" +  file, "r")
        for line in f3:
          # One detected recombinant per line, aggregate all lines in each file
          descendants.write(line)
        f3.close()
recombinants.close()
descendents.close()
unfiltered_recombinants.close()

# Remove temp directory 
#subprocess.run(["rm", "-r", temp])
print("Filtered recombination events results written to {}/recombinants_{}.txt".format(local_results,date))

# Check to make sure Chronumental job finished successfully
while(p1.is_alive()):
    print("Chronumental job not finished running yet.")
    time.sleep(20)

# Rank recombinants and generate final recombinant node information output file
filtration_results_file = "{}/recombinants_{}.txt".format(local_results,date)
# Assumes Chronumental inferred dates file output to current directory
chron_dates_file = "chronumental_dates_{}.tsv".format(metadata)
recomb_output_file = "{}/final_recombinants_{}.txt".format(local_results, date)

cmd = post_process(mat, filtration_results_file, chron_dates_file, str(config["date"]), recomb_output_file)
print(cmd)
subprocess.run(cmd)

# Copy over final results file to GCP storage bucket
subprocess.run(["gsutil", "cp", recomb_output_file, results])

print("Final recombination event results written to {}/recombinants_{}.txt".format(local_results,date))
