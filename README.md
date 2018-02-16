# sea_zrh_rs

Code for post-processing fmri data that has been pre-processed with
fmriprep.

Docker container is [here](https://hub.docker.com/r/fliem/sea_zrh_rs/builds/)

## Usage


### The basic command:

    docker run --rm -ti \
    fliem/sea_zrh_rs:{version} \
    fmriprep_dir output_dir analysis_level


### analysis_levels
The tool has several analysis levels

* **participant_sbc_motion_cc**:
produces maps with seedbased correlation to test confound regression
with FD + motion and CompCor regressors
* **group_collect_motion**:
aggregates mean and max FD and creates distribution plots


### A full docker command could look like:
Let's assume that the fmriprep data is in `/project/fmriprep` and
we want so save post-processing output to `/project/rs_postprocessing`.

    docker run --rm -ti \
    -v /project/fmriprep:/data/in \
    -v /project/rs_postprocessing:/data/out \
    fliem/sea_zrh_rs:dev3 \
    /data/in /data/out participant_sbc_motion_cc \
    --TR 2 --n_cpus 4

* For denoising TR in seconds is required
* By default all subjects found are processed. To restrict subjects,
specify `--participant_label s01 s02`
* Analsysis are run in parallel if `--n_cpus` is given.


## Full usage
    usage: run.py [-h]
                  [--participant_label PARTICIPANT_LABEL [PARTICIPANT_LABEL ...]]
                  [--n_cpus N_CPUS] [--TR TR]
                  fmriprep_dir output_dir
                  {participant_sbc_motion_cc,group_collect_motion}

    SEA ZRH RS analysis code

    positional arguments:
      fmriprep_dir          The directory with the fmriprep-preprocessed data.
      output_dir            The directory where the output files will be written
                            to. Can be the same base dir for all analysis levels,
                            as subdirectories are created
      {participant_sbc_motion_cc,group_collect_motion}
                            Level of the analysis that will be performed.
                            *"participant_sbc_motion_cc": produces maps with
                            seedbased correlation to test confound regression with
                            FD + motion and CompCor regressors *"group_collect_motion":
                            aggregates mean and max FD and creates distribution
                            plots

    optional arguments:
      -h, --help            show this help message and exit
      --participant_label PARTICIPANT_LABEL [PARTICIPANT_LABEL ...]
                            The label of the participant that should be analyzed.
                            The label corresponds to sub-<participant_label> from
                            the BIDS spec (so it does not include "sub-"). If this
                            parameter is not provided all subjects should be
                            analyzed. Multiple participants can be specified with
                            a space separated list. (default: None)
      --n_cpus N_CPUS       Number of CPUs/cores available to use. If not defined
                            use 1 core (default: 1)
      --TR TR               TR of your data in seconds. (default: None)
