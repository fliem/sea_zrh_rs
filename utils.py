import os
from glob import glob
import itertools
import matplotlib

matplotlib.use('Agg')

import pandas as pd


def get_subject_sessions(input_dir, participant_label, raise_if_empty=True, fmriprep_dir=True):
    """
    returns a list of subjects and a list of subject, session tuples
    * input_dir can be an fmriprep_dir; set fmriprep_dir=True to filter out sessions that only have anat, but no func
    * input_dir can be any hierarchy of type sub-01/ses-11; then set fmriprep_dir=False
    """
    os.chdir(input_dir)

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

    if fmriprep_dir:
        # filter out sessions that only have anat, but no func:
        subjects_sessions = list(filter(lambda s: os.path.isdir(os.path.join(input_dir, "sub-" + s[0], "ses-" + s[1],
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


def get_spikereg_confounds(motion_ts, threshold):
    """
    motion_ts = [0.1, 0.7, 0.2, 0.6, 0.3]
    threshold = 0.5
    get_spikereg_confounds(motion_ts, threshold)

    returns
    1.) a df with spikereg confound regressors (trs with motion > thershold)
       outlier_1  outlier_2
    0          0          0
    1          1          0
    2          0          0
    3          0          1
    4          0          0

    2.) a df with counts of outlier and non-outlier trs
       outlier  n_tr
    0    False     3
    1     True     2
    """
    df = pd.DataFrame({"motion": motion_ts})
    df.fillna(value=0, inplace=True)  # first value is nan
    df["outlier"] = df["motion"] > threshold
    outlier_stats = df.groupby("outlier").count().reset_index().rename(columns={"motion": "n_tr"})

    df["outliers_num"] = 0
    df.loc[df.outlier, "outliers_num"] = range(1, df.outlier.sum() + 1)
    outliers = pd.get_dummies(df.outliers_num, dtype=int, drop_first=True, prefix="outlier")

    return outliers, outlier_stats


def get_confounds(confounds_file, kind="36P", spikereg_threshold=None):
    """
    takes a fmriprep confounds file and creates data frame with regressors.
    kind == "36P" returns Satterthwaite's 36P confound regressors
    kind == "9P" returns CSF, WM, Global signal + 6 motion parameters (used in Ng et al., 2016)

    if spikereg_threshold=None, no spike regression is performed

    Satterthwaite, T. D., Elliott, M. A., Gerraty, R. T., Ruparel, K., Loughead, J., Calkins, M. E., et al. (2013).
    An improved framework for confound regression and filtering for control of motion artifact in the preprocessing
    of resting-state functional connectivity data. NeuroImage, 64, 240–256.
    http://doi.org/10.1016/j.neuroimage.2012.08.052

    Ng et al. (2016). http://doi.org/10.1016/j.neuroimage.2016.03.029
    """
    if kind not in ["36P", "9P", "6P"]:
        raise Exception("Confound type unknown {}".format(kind))

    df = pd.read_csv(confounds_file, sep="\t")

    p6 = df[['X', 'Y', 'Z', 'RotX', 'RotY', 'RotZ']]
    p9 = df[['CSF', 'WhiteMatter', 'GlobalSignal', 'X', 'Y', 'Z', 'RotX', 'RotY', 'RotZ']]
    p9_der = p9.diff().fillna(0)
    p9_der.columns = [c + "_der" for c in p9_der.columns]
    p18 = pd.concat((p9, p9_der), axis=1)
    p18_2 = p18 ** 2
    p18_2.columns = [c + "_2" for c in p18_2.columns]
    p36 = pd.concat((p18, p18_2), axis=1)

    if kind == "36P":
        confounds = p36
    elif kind == "9P":
        confounds = p9
    elif kind == "6P":
        confounds = p6

    if spikereg_threshold:
        threshold = spikereg_threshold
    else:
        # if no spike regression still call get_spikereg_confounds to get count of available trs
        threshold = 99999
    outliers, outlier_stats = get_spikereg_confounds(df["FramewiseDisplacement"].values, threshold)

    if spikereg_threshold:
        confounds = pd.concat([confounds, outliers], axis=1)

    return confounds, outlier_stats


def test_spikereg():
    confounds_file = os.path.join("test_data/sub-1_ses-1_task-rest_run-1_bold_confounds.tsv")
    df = pd.read_csv(confounds_file, sep="\t")
    spikereg_threshold = 0.4
    outliers, outlier_stats = get_spikereg_confounds(df["FramewiseDisplacement"].values, spikereg_threshold)
    assert outliers.columns.tolist() == ['outlier_1', 'outlier_2'], "spikereg outlier count not correct"


def test_spikereg_outlier_stats():
    confounds_file = os.path.join("test_data/sub-1_ses-1_task-rest_run-1_bold_confounds.tsv")
    df = pd.read_csv(confounds_file, sep="\t")
    spikereg_threshold = 0.4
    outliers, outlier_stats = get_spikereg_confounds(df["FramewiseDisplacement"].values, spikereg_threshold)
    assert outlier_stats.sum()["n_tr"] == df.shape[0], "outlier stats count not correct"


def test_36p():
    confounds_file = os.path.join("test_data/sub-1_ses-1_task-rest_run-1_bold_confounds.tsv")
    confounds, outlier_stats = get_confounds(confounds_file)
    assert confounds.shape == (225, 36), "Shape of 36P df is wrong {}".format(confounds.shape)


def test_36p_spikereg():
    confounds_file = os.path.join("test_data/sub-1_ses-1_task-rest_run-1_bold_confounds.tsv")
    spikereg_threshold = 0.4
    confounds, outlier_stats = get_confounds(confounds_file, spikereg_threshold=spikereg_threshold)
    assert confounds.shape == (225, 38), "Shape of 36P+spikereg df is wrong {}".format(confounds.shape)


def test_36p_spikereg2():
    confounds_file = os.path.join("test_data/sub-1_ses-1_task-rest_run-1_bold_confounds.tsv")
    spikereg_threshold = 0.5
    confounds, outlier_stats = get_confounds(confounds_file, spikereg_threshold=spikereg_threshold)
    assert confounds.shape == (225, 36), "Shape of 36P+spikereg df is wrong {}".format(confounds.shape)


def test_9p():
    confounds_file = os.path.join("test_data/sub-1_ses-1_task-rest_run-1_bold_confounds.tsv")
    confounds, outlier_stats = get_confounds(confounds_file, "9P")
    assert confounds.shape == (225, 9), "Shape of 9P df is wrong {}".format(confounds.shape)


def save_feather(df, filename):
    "saves df to feather file after resetting index"
    if not filename.endswith(".feather"):
        filename += ".feather"
    df.index.name = "index"
    df.reset_index(inplace=True)
    df.columns = df.columns.astype(str)
    df.to_feather(filename)
