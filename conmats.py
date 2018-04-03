import os
from nilearn import input_data, connectome, datasets, plotting
from utils import get_36P_confounds, get_files
import numpy as np
import pandas as pd
import matplotlib

matplotlib.use('Agg')
from matplotlib import pyplot as plt


def _get_roi_info(parc):
    if parc == "msdl":
        atlas = datasets.fetch_atlas_msdl()
        roi_file = atlas['maps']
        df_labels = pd.DataFrame({"roi_labels": atlas['labels']})
        if isinstance(df_labels["roi_labels"][0], bytes):
            df_labels["roi_labels"] = df_labels.roi_labels.apply(bytes.decode)
        roi_names = df_labels["roi_labels"].values
        roi_type = "maps"

    elif parc == "gordon":
        # fixme
        gordon_dir = "/parcs/Gordon/Parcels"
        roi_file = os.path.join(gordon_dir, "Parcels_MNI_111.nii")
        labs_df = pd.read_excel(os.path.join(gordon_dir, "Parcels.xlsx"))
        roi_names = labs_df.ParcelID.values
        roi_type = "labels"

    elif parc == "basc197":
        atlas = datasets.fetch_atlas_basc_multiscale_2015(version='sym')
        roi_file = atlas['scale197']
        roi_names = np.arange(1, 198).astype(int)
        roi_type = "labels"

    elif parc == "basc444":
        atlas = datasets.fetch_atlas_basc_multiscale_2015(version='sym')
        roi_file = atlas['scale444']
        roi_names = np.arange(1, 445).astype(int)
        roi_type = "labels"

    else:
        raise Exception("Parcellation not known {}".format(parc))
    return roi_file, roi_names, roi_type


raw_mat = np.array([
    (1, .3, .7, .2),
    (.3, 1, .6, .5),
    (.7, .6, 1, .9),
    (.2, .5, .9, 1),
])
roi_names = ["r1", "r2", "r3", "r4"]


def test_get_con_df():
    df = _get_con_df(raw_mat, roi_names)
    assert df.shape == (4, 4), "shape issue"
    assert (df.columns == roi_names).all(), "name issue"


def _get_con_df(raw_mat, roi_names):
    """
    takes a symmetrical connectivity matrix (e.g., numpy array) and a list of roi_names (strings)
    returns data frame with roi_names as index and column names
    e.g.
         r1   r2   r3   r4
    r1  0.0  0.3  0.7  0.2
    r2  0.3  0.0  0.6  0.5
    r3  0.7  0.6  0.0  0.9
    r4  0.2  0.5  0.9  0.0
    """
    # sanity check if matrix is symmetrical
    assert np.allclose(raw_mat, raw_mat.T), "matrix not symmetrical"

    np.fill_diagonal(raw_mat, 0)
    con_df = pd.DataFrame(raw_mat, index=roi_names, columns=roi_names)
    return con_df


def conmat_one_session(subject, session, fmriprep_dir, output_dir, tr, parc):
    full_out_dir = os.path.join(output_dir, "sub-{}".format(subject), "ses-{}".format(session))
    os.makedirs(full_out_dir, exist_ok=True)

    out_stub_short = "{}_{}".format(subject, session)
    out_stub = "sub-{}_ses-{}".format(subject, session)
    conmat_file = os.path.join(full_out_dir, "{}_parc-{}_conf-36P_conmat.tsv".format(out_stub, parc))
    conmat_plot_file = os.path.join(full_out_dir, "{}_parc-{}_conf-36P_conmat.png".format(out_stub, parc))
    report_file = os.path.join(full_out_dir, "{}_parc-{}_conf-36P_report.txt".format(out_stub, parc))

    if not (os.path.exists(conmat_file) and os.path.exists(conmat_plot_file) and os.path.exists(report_file)):
        print("*** Calc conmats for {} {} {} ***".format(subject, session, parc))

        confounds_file, brainmask_file, rs_file, anat_file = get_files(fmriprep_dir, subject, session)
        roi_file, roi_names, roi_type = _get_roi_info(parc)

        conmat, report_str = extract_mat(rs_file, brainmask_file, roi_file, confounds_file, roi_type, tr)
        conmat_df = _get_con_df(conmat, roi_names)

        conmat_df.to_csv(conmat_file, sep="\t")
        with open(report_file, "w") as fi:
            fi.write(report_str)

        plotting.plot_matrix(conmat, labels=roi_names, figure=(9, 7), vmax=1, vmin=-1, title=out_stub_short + " r")
        plt.savefig(conmat_plot_file, bbox_inches='tight')
        plt.close()
    else:
        print("*** Conmats for {} {} {} already computed. Do nothing. ***".format(subject, session, parc))


def extract_mat(rs_file, brainmask_file, roi_file, confounds_file, roi_type, tr):
    """
    36 P
    """

    # Masker
    masker_pars = {"mask_img": brainmask_file, "detrend": True, "standardize": True, "low_pass": 0.1, "high_pass": \
        0.01, "t_r": tr}

    if roi_type == "maps":
        # for msdl type probablistic rois
        masker = input_data.NiftiMapsMasker(roi_file, **masker_pars)
    elif roi_type == "labels":
        # for binary rois
        masker = input_data.NiftiLabelsMasker(roi_file, **masker_pars)
    else:
        raise Exception("roi type not known {}".format(roi_type))

    # Extract time series
    confounds = get_36P_confounds(confounds_file)
    time_series = masker.fit_transform(rs_file, confounds=confounds.values)

    con_measure = connectome.ConnectivityMeasure(kind='correlation')
    conmat = con_measure.fit_transform([time_series])[0]

    report_str = "rs_file\t{}\n".format(rs_file)
    report_str += "roi_file\t{}\n".format(roi_file)

    keys = list(masker_pars.keys())
    keys.sort()
    for k in keys:
        report_str += "{}\t{}\n".format(k, masker_pars[k])

    report_str += "\n\n"
    report_str += "confounds\t{}".format(", ".join(confounds.columns))
    report_str += "\n\n"
    report_str += confounds.to_string()

    return conmat, report_str


