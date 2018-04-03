#!/usr/bin/env python3
import argparse
import os
from joblib import Parallel, delayed
from utils import get_subject_sessions, get_motion_ts_one_subject
from conmats import conmat_one_session
from sbc import sbc_one_session, sbc_group
import itertools
import pandas as pd

if __name__ == "__main__":
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                     description='SEA ZRH RS analysis code')
    parser.add_argument('fmriprep_dir', help='The directory with the fmriprep-preprocessed data.')
    parser.add_argument('output_dir', help='The directory where the output files will be written to. '
                                           'Can be the same base dir for all analysis levels, as subdirectories are '
                                           'created')
    parser.add_argument('analysis_level', help='''Level of the analysis that will be performed. 
                                               *"participant_1_sbc_pcc": produces maps with seedbased correlation
                                               *"group_1_sbc_pcc": mean sbc maps
                                               *"participant_2_conmats": connectivity matrix extraction on the subject 
                                               level
                                               *"group_2_collect_motion": motion ts. 
                                               Info from all participants are collected into one file 
                                                '''
                        , choices=['participant_1_sbc_pcc',
                                   'group_1_sbc_pcc',
                                   'participant_2_conmats',
                                   'group_2_collect_motion'])

    parser.add_argument('--participant_label',
                        help='The label of the participant that should be analyzed. The label '
                             'corresponds to sub-<participant_label> from the BIDS spec '
                             '(so it does not include "sub-"). If this parameter is not '
                             'provided all subjects should be analyzed. Multiple '
                             'participants can be specified with a space separated list.', nargs="+")

    parser.add_argument('--n_cpus', help='Number of CPUs/cores available to use. If not defined use 1 core', default=1,
                        type=int)

    parser.add_argument('--TR', help='TR of your data in seconds.', type=float)

    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)
    os.chdir(args.output_dir)

    subjects, subjects_sessions = get_subject_sessions(args.fmriprep_dir, args.participant_label)
    print("Processing {} subjects and a total of {} sessions".format(len(subjects), len(subjects_sessions)))

    parc_list = ["msdl", "gordon", "basc197", "basc444"]

    if args.analysis_level == "participant_1_sbc_pcc":

        if not args.TR:
            raise Exception("TR required for this step. Stopping")

        output_dir = os.path.join(args.output_dir, "sbc", "participant")
        os.makedirs(output_dir, exist_ok=True)

        _ = Parallel(n_jobs=args.n_cpus)(
            delayed(sbc_one_session)(*suse, args.fmriprep_dir, output_dir, args.TR) for suse in
            subjects_sessions)

    elif args.analysis_level == "group_1_sbc_pcc":
        input_dir = os.path.join(args.output_dir, "sbc", "participant")
        output_dir = os.path.join(args.output_dir, "sbc", "group")
        sbc_group(input_dir, output_dir)

    elif args.analysis_level == "participant_2_conmats":
        if not args.TR:
            raise Exception("TR required for this step. Stopping")

        output_dir = os.path.join(args.output_dir, "conmats", "participant")
        os.makedirs(output_dir, exist_ok=True)

        subjects_sessions_parc = list(itertools.product(subjects_sessions, parc_list))
        _ = Parallel(n_jobs=args.n_cpus)(
            delayed(conmat_one_session)(suse[0][0], suse[0][1], args.fmriprep_dir, output_dir, args.TR, suse[1])
            for suse in subjects_sessions_parc)

    elif args.analysis_level == "group_2_collect_motion":
        output_dir = os.path.join(args.output_dir, "motion", "group")
        os.makedirs(output_dir, exist_ok=True)

        subjects, subjects_sessions = get_subject_sessions(args.fmriprep_dir, args.participant_label)

        # collect motion time series
        motion_dfs = []
        for subject, session in subjects_sessions:
            motion_dfs.append(get_motion_ts_one_subject(subject, session, args.fmriprep_dir))
        motion_df = pd.concat(motion_dfs)
        motion_df.reset_index(inplace=True, drop=True)
        out_file = os.path.join(output_dir, "group_motion_ts.tsv")
        motion_df.to_csv(out_file, sep="\t")

    else:
        raise NotImplementedError(args.analysis_level)
