![Ripples](images/ripples_logo.png)

# **Recombination Inference using Phylogenetic PLacEmentS** (RIPPLES)

RIPPLES is a program to detect recombination events in large mutation annotated tree (MAT) files.  This repo contains a workflow for running RIPPLES on Google Cloud Platform.

Please also so this tutorial for various applications: [RIPPLES Tutorial](https://usher-wiki.readthedocs.io/en/latest/ripples.html)

# RIPPLES on Google Cloud Platform

<br>

## Setup your Google Cloud Platform Account
___

1. **Setup Cloud Console:**
	- If you don't already have a Google Cloud Platform account, please follow these instructions to open a Cloud Console, create a project name (`bucket_id`) and project ID (`project_id`):
       	[Installation and Setup](https://cloud.google.com/deployment-manager/docs/step-by-step-guide/installation-and-setup)

<br>

2. **Add a service account**:
	- TODO: Add instructions for setting up service account.  ADD link also.

<br>

3. **Download Keys JSON file**
	- Once you have created a service account, you need to add keys to this serivce account.
	Click the Navagation Menu side bar on the Web Console and go to `IAM & Admin` -> `Service Accounts` and click on the active service account you just created from the previous step.

	- Click the tab the says `Keys` and click `ADD KEY` and `Create new key`.  Select `JSON` Key type.  A new `<key>.json` file will automatically be downloaded from your browser.

	- Move this downloaded `<key>.json` file to the following location:

		```
		~/.config/gcloud/<key>.json
		```

	-  Then run the following command in your terminal to set environment variable path to location where you just placed your downloaded `<keys>.json` file.

		```
		KEY=~/.config/gcloud/<keys>.json
		```

<br>


## Run RIPPLES Docker Container
___

Pull and run public RIPPLES Docker image with the following command:
```
docker run -it -e GOOGLE_APPLICATION_CREDENTIALS=/tmp/keys/<keys>.json -v ${KEY}:/tmp/keys/<keys>.json:ro mrkylesmith/ripples-pipeline:latest
```

This will drop you into an interactive (`-it` flag) Docker container shell where you will be able to launch RIPPLES jobs to GCP from. The docker image is configured for you with all the necessary packages to run RIPPLES, so no additional manual installs will be required by you.

<br>


## Setup Access to GCP Account and APIs
___

- Add your bucket ID, project ID, and the downloaded access keys file to the `ripples.yaml` configuration file like this: 
	- **Make sure the access credentials JSON file you downloaded is at the path below, or update the path in `ripples.yaml`**
```
bucket_id: <your_bucket_id> 
project_id: <your_project_id> 
key_file: /tmp/keys/<your_key_file.json>
```

3. Run the following setup program, which will activate your access credentials and enable all the necessary GCP APIs
```
python3 setup.py
```

<br>

## Configure RIPPLES job to execute and GCP Machine Type
___

Set configurations for the current RIPPLES job you want to run in `ripples.yaml` , shown below:
```
 # Ripples parameters config
 version: ripples
 mat: <mat.pb>
 results: <results/>
 
 # GCP machine and Storage Bucket config
 instances: 4
 boot_disk_size: 30
 machine_type: e2-standard-16
 logging: <example.log>

```
### RIPPLES Options:
- `version`: Set as `ripples` or `ripples-fast` to run fast verison of RIPPLES program (see trade-offs of `ripples-fast` below).

- `mat`: The Mutation Annotated Tree (MAT) protobuf that you want to search for recombination events. 

- `results`: The output directory remote path on GCP Storage Bucket where RIPPLES will output results. 

### GCP Machine Type Options:
- `instances`: Number of GCP machines that RIPPLES will be automatically parallelized across.  Results will be automatically aggregated into `results` directory on your GCP Storage Bucket when all RIPPLES jobs are complete.

- `boot_disk_size`: Instance startup disk size. **Leave as 30GB**.

- `machine_type`: Select the GCP machine type you would like to launch.  

- `logging`: Name of the logging file for this particular RIPPLES job (**change for different jobs**) that will be output into your GCP Storage bucket under `bucket_id/logging/<logging>`.

<br>

## Running Ripples on GCP
___

Execute the following command to launch your RIPPLES job according to the set configurations in `ripples.yaml`.
```
python3 run.py
```
Once all jobs are complete, the dectected recombinants will be placed in the specified `results/` directory in your GCP Storage Bucket.





