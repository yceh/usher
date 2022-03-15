#!/bin/python3

import subprocess
import sys
import time
import os
import pathlib
import re
import json
import time
import datetime
import yaml

def get_config():
    with open('ripples.yaml') as f:
      data = yaml.load(f, Loader=yaml.FullLoader)
      return data 

def setup():
    cmd1 = ["gcloud", "config", "set", "project", PROJECT_ID]
    cmd2 = ["gcloud", "auth", "activate-service-account", "--key-file", KEY_FILE]
    print("Setting up GCP access through gcloud.")
    subprocess.run(cmd1)
    subprocess.run(cmd2)

def get_partitions(nodes, instances):
    w = nodes//instances
    partitions = [ (i,min(i+w-1,nodes-1)) for i in range(0,nodes,w) ]
    return partitions

def convert(n):
    return str(datetime.timedelta(seconds = n))

def parse_command(MAT, start, end, out):
    #TODO: Verify correct path inside of docker container here
    command = "python3 recombination/filtering/process.py {} {} {} {} {} {}".format(VERSION,
            MAT, start, end, out, BUCKET_ID, RESULTS)
    return command

def gcloud_run(command):
    cmd = ["gcloud", "beta", "lifesciences", "pipelines", "run",
        "--location", "us-central1",
        "--regions", "us-central1",
        "--machine-type", MACHINE_TYPE, 
        "--boot-disk-size", BOOT_DISK_SIZE,
        "--logging", LOGGING,
        "--docker-image", docker_image,
        "--command-line", command,
        #"--outputs", RESULTS,
        "--format=json"]

    print(" ".join(cmd))
    out = subprocess.check_output(cmd)
    result = json.loads(out)
    name = result['name']
    id = re.search("operations\/(\d+)", name).groups()[0]
    return {'operation_id': id, 'result': result}

def gcloud_describe(operation_id):
     cmd = ["gcloud", "beta", "lifesciences",
            "operations", "describe", "--format=json", operation_id]
     out = subprocess.check_output(cmd)
     result = json.loads(out)
     done = 'done' in result and result['done'] == True
     return {'done': done, 'result': result}

# Configs and credentials from yaml
config = get_config()

# Set global config options
BUCKET_ID = config["BUCKET_ID"]
PROJECT_ID = config["PROJECT_ID"]  
KEY_FILE = config["KEY_FILE"]
VERSION = config["VERSION"]
BOOT_DISK_SIZE = config["BOOT_DISK_SIZE"]
MACHINE_TYPE = config["MACHINE_TYPE"]
LOGGING = "gs://{}/logging/{}".format(BUCKET_ID, config["LOGGING"])
# Remote location in GCP Storge Bucket to copy filtered recombinants
#TODO: -> Make sure you have them create a directory of this name on GCP
RESULTS = "gs://{}/{}".format(BUCKET_ID, config["RESULTS"])
docker_image = "mrkylesmith/ripples-pipeline:latest"

# Activate credentials for GCP Console and gcloud util
setup()

#Using setup script
#project_id = os.getenv('PROJECT_ID')     
#key_file = os.getenv('KEY_FILE')
#docker_image = os.getenv('DOCKER_IMAGE') 

#TODO: Need this information from ripples ahead of time
# Total num of long branches to search
nodes = 47000               
# Num of paritions to split of p instances 
INSTANCES = config["INSTANCES"]
# Num of long branches to search per instance
w = nodes//instances        

partitions = get_partitions(nodes, instances)
print("nodes: {}, instances: {}, w: {}".format(nodes, instances, w))
print("partitions: {}".format(partitions))

processes = []
for partition in partitions:

    start = str(partition[0])
    end = str(partition[1])
    out = "{}_{}".format(start, end)
    
    command = parse_command(MAT, start, end, out)
    # python3 /data/process.py optimized-large-radius-pruneCatExcludeB30.usher.no-long-branches.pb 0 3 0_3

    info = gcloud_run(command)
    processes.append({'partition': partition, 'operation_id': info['operation_id']})


completed = 0
while completed < instances:
   for process in processes:
     partition = process['partition']
     operation_id = process['operation_id']
     info = gcloud_describe(operation_id)
     done = info['done']
     if done:
       completed += 1
     print("partition: {}, operation_id: {}, done: {}".format(partition, operation_id, done))
   print("{} of {} completed".format(completed, instances))
   time.sleep(1)

#print("All instance jobs have finished.  Merging final ouput.")
#TODO: Deal with merging outputs


