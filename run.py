#!/usr/bin/env python3
import argparse
import os
from joblib import Parallel, delayed
from utils import get_subject_sessions, collect_motion
from connectivity import sbc_one_subject

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                 description='SEA ZRH RS analysis code')
parser.add_argument('fmriprep_dir', help='The directory with the fmriprep-preprocessed data.')
parser.add_argument('output_dir', help='The directory where the output files will be written to. '
                                       'Can be the same base dir for all analysis levels, as subdirectories are '
                                       'created')
parser.add_argument('analysis_level', help='''Level of the analysis that will be performed. 
                                           *"participant_sbc_motion_cc": produces maps with seedbased correlation to 
                                               test confound regression with FD and CompCor regressors
                                           *"group_collect_motion": aggregates mean and max FD and creates 
                                           distribution plots'''
                    , choices=['participant_sbc_motion_cc', 'group_collect_motion'])

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

os.chdir(args.output_dir)

subjects, subjects_sessions = get_subject_sessions(args.fmriprep_dir, args.participant_label)
print("Processing {} subjects and a total of {} sessions".format(len(subjects), len(subjects_sessions)))

if args.analysis_level == "participant_sbc_motion_cc":

    if not args.TR:
        raise Exception("TR required for this step. Stopping")

    full_out_dir = os.path.join(args.output_dir, "participant_1_sbc__motion_cc")
    os.makedirs(full_out_dir, exist_ok=True)

    regressors = ['FramewiseDisplacement', 'aCompCor01', 'aCompCor02', 'aCompCor03', 'aCompCor04', 'aCompCor05', 'X',
                  'Y', 'Z', 'RotX', 'RotY', 'RotZ']
    _ = Parallel(n_jobs=args.n_cpus)(
        delayed(sbc_one_subject)(*suse, args.fmriprep_dir, full_out_dir, args.TR, regressors) for suse in
        subjects_sessions)

elif args.analysis_level == "group_collect_motion":
    full_out_dir = os.path.join(args.output_dir, "group_1_motion")
    os.makedirs(full_out_dir, exist_ok=True)

    collect_motion(subjects_sessions, args.fmriprep_dir, full_out_dir, args.n_cpus)

else:
    raise NotImplementedError(args.analysis_level)