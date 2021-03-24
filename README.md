## Running the pipeline
The typical command for running the pipeline is as follows:

```bash
nextflow run dolphinnext/barcodseq:1.0 --profile docker --reads '*_R{1,2}.fastq.gz' --mate 'pair'
```

This will launch the pipeline with the `docker` configuration profile. See below for more information about profiles.

Note that the pipeline will create the following files in your working directory:

```bash
work            # Directory containing the nextflow working files
results         # Finished results 
.nextflow_log   # Log file from Nextflow
```

### Updating the pipeline
When you run the above command, Nextflow automatically pulls the pipeline code from GitHub and stores it as a cached version. When running the pipeline after this, it will always use the cached version. In order to download latest version of the pipeline you need to run following command:

```bash
nextflow pull dolpinnext/barcodeseq:1.0
```
