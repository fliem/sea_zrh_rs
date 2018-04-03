# sea_zrh_rs

Code for post-processing fmri data that has been pre-processed with
fmriprep.

Docker container is [here](https://hub.docker.com/r/fliem/sea_zrh_rs/builds/)

## Usage


### The basic command:

    docker run --rm -ti \
    fliem/sea_zrh_rs:{version} \
    fmriprep_dir output_dir analysis_level


### A full docker command could look like:
Let's assume that the fmriprep data is in `/project/fmriprep` and
we want so save post-processing output to `/project/rs_postprocessing`.

`{version}` takes the version of the docker container you want to run.
E.g., `fliem/sea_zrh_rs:v1.0.0`


    docker run --rm -ti \
    -v /project/fmriprep:/data/in \
    -v /project/rs_postprocessing:/data/out \
    fliem/sea_zrh_rs:{version} \
    /data/in /data/out participant_1_sbc_pcc \
    --TR 2 --n_cpus 4

* For denoising TR in seconds is required
* By default all subjects found are processed. To restrict subjects,
specify `--participant_label s01 s02`
* Analysis are run in parallel if `--n_cpus` is given.


## Processing steps (analysis levels)
### participant_1_sbc_pcc
Extracts PCC SBC maps for each participant/session.

    docker run --rm -ti \
    -v /project/fmriprep:/data/in \
    -v /project/rs_postprocessing:/data/out \
    fliem/sea_zrh_rs:{version} \
    /data/in /data/out participant_1_sbc_pcc \
    --TR 2 --n_cpus 4

### group_1_sbc_pcc
Calculates mean SBC maps.

    docker run --rm -ti \
    -v /project/fmriprep:/data/in \
    -v /project/rs_postprocessing:/data/out \
    fliem/sea_zrh_rs:{version} \
    /data/in /data/out group_1_sbc_pcc \
    --n_cpus 4

## participant_2_conmats
Extracts conmats for each subject/session.

    docker run --rm -ti \
    -v /project/fmriprep:/data/in \
    -v /project/rs_postprocessing:/data/out \
    fliem/sea_zrh_rs:{version} \
    /data/in /data/out participant_2_conmats \
    --TR 2 --n_cpus 4

## group_2_collect_motion
Collects motion time series for all subjects in one file.

    docker run --rm -ti \
    -v /project/fmriprep:/data/in \
    -v /project/rs_postprocessing:/data/out \
    fliem/sea_zrh_rs:{version} \
    /data/in /data/out group_2_collect_motion

## Full usage
    usage: run.py [-h]
                  [--participant_label PARTICIPANT_LABEL [PARTICIPANT_LABEL ...]]
                  [--n_cpus N_CPUS] [--TR TR]
                  fmriprep_dir output_dir
                  {participant_1_sbc_pcc,group_1_sbc_pcc,participant_2_conmats,group_2_collect_motion}

    SEA ZRH RS analysis code

    positional arguments:
      fmriprep_dir          The directory with the fmriprep-preprocessed data.
      output_dir            The directory where the output files will be written
                            to. Can be the same base dir for all analysis levels,
                            as subdirectories are created
      {participant_1_sbc_pcc,group_1_sbc_pcc,participant_2_conmats,group_2_collect_motion}
                            Level of the analysis that will be performed.
                            *"participant_1_sbc_pcc": produces maps with seedbased
                            correlation *"group_1_sbc_pcc": mean sbc maps
                            *"participant_2_conmats": connectivity matrix
                            extraction on the subject level
                            *"group_2_collect_motion": motion ts. Info from all
                            participants are collected into one file

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


