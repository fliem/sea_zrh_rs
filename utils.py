import os
from glob import glob
import itertools
import matplotlib

matplotlib.use('Agg')

import pandas as pd
import seaborn as sns
from joblib import Parallel, delayed
from matplotlib import pyplot as plt


def get_subject_sessions(fmriprep_dir, participant_label, raise_if_empty=True):
    """
    returns a list of subjects and a list of subject, session tuples
    """
    os.chdir(fmriprep_dir)

    # get subjects
    if participant_label:
        subjects = list(map(lambda s: glob("sub-{}".format(s)), participant_label))
        subjects = list(itertools.chain.from_iterable(subjects))
    else:
        subjects = glob("sub-*")

    subjects = list(filter(lambda s: not s.endswith(".html"), subjects))
    subjects = list(map(lambda s: s.split("-")[-1], subjects))
    subjects.sort()

    # get subject_sessions
    def get_subject_sessions_one_sub(subject):
        sessions = glob("sub-{}/ses-*".format(subject))
        sessions.sort()
        return list(map(lambda s: (subject, s.split("-")[-1]), sessions))

    # subjects_sessions list of tupels: [('s01', 'tp1'), ('s01', 'tp2'), ('s02', 'tp1')]
    subjects_sessions = list(map(get_subject_sessions_one_sub, subjects))
    subjects_sessions = list(itertools.chain.from_iterable(subjects_sessions))

    # filter out sessions that only have anat, but no func:
    subjects_sessions = list(filter(lambda s: os.path.isdir(os.path.join(fmriprep_dir, "sub-" + s[0], "ses-" + s[1],
                                                                         "func")), subjects_sessions))

    print(subjects, subjects_sessions)
    if raise_if_empty:
        if not subjects:
            raise Exception("No subjects found")
        if not subjects_sessions:
            raise Exception("No Sessions found")

    return subjects, subjects_sessions


def check_and_return_path(search_str):
    l = glob(search_str)
    if len(l) != 1:
        raise Exception("glob found more or less than one file: {}. {}".format(search_str, l))
    else:
        return l[0]


def get_files(fmriprep_dir, subject, session):
    subject_dir = os.path.join(fmriprep_dir, "sub-" + subject)
    session_dir = os.path.join(subject_dir, "ses-" + session)

    confounds_file = check_and_return_path(
        os.path.join(session_dir, "func", "sub-*_task-rest_run-1_bold_confounds.tsv"))
    brainmask_file = check_and_return_path(
        os.path.join(session_dir, "func", "sub-*_task-rest_run-1_bold_space-MNI152NLin2009cAsym_brainmask.nii.gz"))
    rs_file = check_and_return_path(
        os.path.join(session_dir, "func", "sub-*_task-rest_run-1_bold_space-MNI152NLin2009cAsym_preproc.nii.gz"))
    try:
        anat_file = check_and_return_path(
            os.path.join(subject_dir, "anat", "sub-*_T1w_space-MNI152NLin2009cAsym_preproc.nii.gz"))
    except:  # might be a case with only 1 session, there the anat files are under sub-xx/ses-xx/anat/...
        anat_file = check_and_return_path(
            os.path.join(session_dir, "anat", "sub-*_T1w_space-MNI152NLin2009cAsym_preproc.nii.gz"))
    return confounds_file, brainmask_file, rs_file, anat_file


def get_motion_one_subject(subject, session, fmriprep_dir):
    # returns data frame with FD_mean, FD_median, FD_max, subject, session
    confounds_file, brainmask_file, rs_file, anat_file = get_files(fmriprep_dir, subject, session)
    df = pd.read_csv(confounds_file, sep="\t")
    agg = df[['FramewiseDisplacement']].agg(["mean", "median", "max"])
    agg.rename({"mean": "FD_mean", "median": "FD_median", "max": "FD_max"}, inplace=True)
    agg = agg.T
    agg["subject"] = subject
    agg["session"] = session
    agg.reset_index(drop=True, inplace=True)
    return agg


def get_motion_ts_one_subject(subject, session, fmriprep_dir):
    # returns data frame with FD time series
    confounds_file, brainmask_file, rs_file, anat_file = get_files(fmriprep_dir, subject, session)
    df = pd.read_csv(confounds_file, sep="\t")

    frames = df[['FramewiseDisplacement']].copy()
    frames["subject"] = subject
    frames["session"] = session
    frames.reset_index(inplace=True)
    frames.rename(columns={"index": "tr"}, inplace=True)
    return frames


def collect_motion(subjects_sessions, fmriprep_dir, full_out_dir, n_cpus):
    # aggreagate FD
    dfs = Parallel(n_jobs=n_cpus)(
        delayed(get_motion_one_subject)(*suse, fmriprep_dir) for suse in subjects_sessions)

    df = pd.concat(dfs)
    df = df[['subject', 'session', 'FD_mean', 'FD_median', 'FD_max']]
    df.sort_values(by=["subject", "session"], inplace=True)
    out_file = os.path.join(full_out_dir, "group_motion.tsv")
    df.to_csv(out_file, sep="\t", index=False)
    print("Writing to {}".format(out_file))

    # plots
    for m in ["FD_mean", "FD_median", "FD_max"]:
        sns.distplot(df[m])
        plt.savefig(os.path.join(full_out_dir, m + "_hist.pdf"))
        plt.close()
        sns.boxplot(y=df[m])
        plt.savefig(os.path.join(full_out_dir, m + "_box.pdf"))
        plt.close()

    # export FD time series
    dfs = Parallel(n_jobs=n_cpus)(
        delayed(get_motion_ts_one_subject)(*suse, fmriprep_dir) for suse in subjects_sessions)
    df = pd.concat(dfs)
    df.sort_values(by=["subject", "session"], inplace=True)
    out_file = os.path.join(full_out_dir, "group_motion_time_series.tsv")
    df.to_csv(out_file, sep="\t", index=False)
    print("Writing time series to {}".format(out_file))
