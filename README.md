# PlastoGram-analysis

This repository contains the data and code necessary to reproduce the results from the article  *Where do they come from, where do they go? Efficient prediction of protein subplastid localization and sequence origin with PlastoGram.* (2022) by Sidorczuk K., Gagat P., Kała J., Nielsen H., Pietluch F., Mackiewicz P., Burdukiewicz M.


## Getting started

This repository uses [renv](https://CRAN.R-project.org/package=renv) and [targets](https://CRAN.R-project.org/package=targets) packages to control the workflow and assure the reproducibility. 

Some of the data files are too large to store them on GitHub but they can be downloaded using the links below:

- [All_sequences.fasta](https://www.dropbox.com/s/y9bo9xl9dgqapb0/All_sequences.fasta?dl=0) - All downloaded sequences of proteins localizing to different plastid compartments. See the supplementary materials of the paper for exact queries used to obtain them.

- [Dataset_annotations_references.xlsx](https://www.dropbox.com/s/wuewgtk6c0mnwrd/Dataset_annotations_references.xlsx?dl=0) - File with UniProt annotations of downloaded proteins, along with our curated localization and references. See the Methods section of the article for detailed description of manual curation procedure.

Part of the analysis requires HMMER software. Please make sure that you have installed HMMER before reproducing the pipeline. Please see the [HMMER documentation](http://hmmer.org/documentation.html) for installation guidelines.

To reproduce the results clone the repo, set your path to the directories with data files and:

``` r
renv::restore()
targets::tar_make()
```


## Content

**\_targets.R** - reproducible pipeline for generation of all data sets and results processing, 

**data** - data files used during the study, e.g. for creation of the positive dataset,

**drafts** - draft codes used for initial exploratory analyses,

**functions** - all functions used for running the pipeline and obtaining results,

**renv** - renv package files,

**third-party** - third-party executables used in the pipeline.

## Pre-calculated results

If you are interested in seeing some intermediate results but do not want to run the whole pipeline, we provide links to the most important directories and files with results obtained during the study.

- [Replication 1](https://www.dropbox.com/s/799p0hehrzn8ys1/All_models_predictions_envelope_rep1.csv?dl=0), [Replication 2](https://www.dropbox.com/s/p2vm6wfp714j06t/All_models_predictions_envelope_rep2.csv?dl=0), [Replication 3](https://www.dropbox.com/s/s6f2w96a4ll1kc1/All_models_predictions_envelope_rep3.csv?dl=0), [Replication 4](https://www.dropbox.com/s/n91w3yzputp8awe/All_models_predictions_envelope_rep4.csv?dl=0), [Replication 5](https://www.dropbox.com/s/x8j135uvk38wja4/All_models_predictions_envelope_rep5.csv?dl=0) - prediction results of all lower-level models in 10-fold CV repeated 5 times.

- [Model_architectures_envelope](https://www.dropbox.com/sh/kjklnuq7tx995ll/AABvZBzURqgN4T8GFOAkE0Gpa?dl=0) - directory with files describing each of the considered ensembles.

- [Model_architectures_envelope_results](https://www.dropbox.com/sh/ocxb587buov1y07/AAB99qhDz2fTtcqttFYdiigra?dl=0) - directory containing lower-level models prediction results filtered according to each ensemble.

- [Architectures_envelope_performance.csv](https://www.dropbox.com/s/92n6q3ahundfmtm/Architectures_envelope_performance.csv?dl=0) - file with performance measures of all ensembles calculated for each replication and fold separately. 

- [Architectures_envelope_mean_performance.csv](https://www.dropbox.com/s/3iv1acs18ofe7lf/Architectures_envelope_mean_performance.csv?dl=0) - file with averaged performance measures of all ensembles.


