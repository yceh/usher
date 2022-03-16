#!/bin/python3
# Simple script to set configuration arguments for GCP
import subprocess
import sys
import os
import pathlib
import json
import yaml


def get_config():
    with open('setup.yaml') as f:
      data = yaml.load(f, Loader=yaml.FullLoader)
      return data 


def cloud_setup():
    cmd1 = ["gcloud", "config", "set", "project", project_id]
    cmd2 = ["gcloud", "auth", "activate-service-account", "--key-file", key_file]
    print("Setting up GCP access through gcloud.")
    subprocess.run(cmd1)
    subprocess.run(cmd2)


config = get_config()
bucket_id = config["bucket_id"]
project_id = config["project_id"]
key_file = config["key_file"]
cloud_setup()


