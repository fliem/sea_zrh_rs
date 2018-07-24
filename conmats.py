import os
from nilearn import input_data, connectome, datasets, plotting
from utils import get_confounds, get_files, save_feather
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
        atlas_dir = "/parcs/Gordon/Parcels"
        roi_file = os.path.join(atlas_dir, "Parcels_MNI_111.nii")
        labs_df = pd.read_excel(os.path.join(atlas_dir, "Parcels.xlsx"))
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

    elif parc == "schaefer200":
        atlas_dir = "/parcs/Schaefer"
        schaefer_cols = "roi community c1 c2 c3 c4".split(" ")
        roi_file = os.path.join(atlas_dir, "Schaefer2018_200Parcels_17Networks_order_FSLMNI152_1mm.nii.gz")
        labs_df = pd.read_csv(os.path.join(atlas_dir, "Schaefer2018_200Parcels_17Networks_order.txt"), sep="\t",
                              names=schaefer_cols)
        roi_names = labs_df.roi
        roi_type = "labels"

    elif parc == "schaefer400":
        atlas_dir = "/parcs/Schaefer"
        schaefer_cols = "roi community c1 c2 c3 c4".split(" ")
        roi_file = os.path.join(atlas_dir, "Schaefer2018_400Parcels_17Networks_order_FSLMNI152_1mm.nii.gz")
        labs_df = pd.read_csv(os.path.join(atlas_dir, "Schaefer2018_400Parcels_17Networks_order.txt"), sep="\t",
                              names=schaefer_cols)
        roi_names = labs_df.roi
        roi_type = "labels"

    elif parc == "yeo17":
        atlas = datasets.fetch_atlas_yeo_2011()
        roi_file = atlas['thick_17']
        yeo_cols = "roi roi_labels c1 c2 c3 c4".split(" ")
        df_labels = pd.read_csv(atlas["colors_17"], sep=r"\s*", engine="python", names=yeo_cols, skiprows=1)
        roi_names = df_labels["roi_labels"].values
        roi_type = "labels"

    elif parc == "yeo17thin":
        atlas = datasets.fetch_atlas_yeo_2011()
        roi_file = atlas['thin_17']
        yeo_cols = "roi roi_labels c1 c2 c3 c4".split(" ")
        df_labels = pd.read_csv(atlas["colors_17"], sep=r"\s*", engine="python", names=yeo_cols, skiprows=1)
        roi_names = df_labels["roi_labels"].values
        roi_type = "labels"

    elif parc == "yeo17split":
        atlas_dir = "/parcs/Yeo_splithemi"
        roi_file = os.path.join(atlas_dir, "yeo_2011_thick17_splithemi.nii.gz")
        labs_df = pd.read_csv(os.path.join(atlas_dir, "yeo_2011_thick17_splithemi.tsv"), sep="\t")
        roi_names = labs_df.full_roi_name.values
        roi_type = "labels"

    elif parc == "yeo7":
        atlas = datasets.fetch_atlas_yeo_2011()
        roi_file = atlas['thick_7']
        yeo_cols = "roi roi_labels c1 c2 c3 c4".split(" ")
        df_labels = pd.read_csv(atlas["colors_7"], sep=r"\s*", engine="python", names=yeo_cols, skiprows=1)
        roi_names = df_labels["roi_labels"].values
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


def conmat_one_session(subject, session, fmriprep_dir, output_dir, tr, conf, parc, spikereg_threshold=None):
    full_out_dir = os.path.join(output_dir, "sub-{}".format(subject), "ses-{}".format(session))
    os.makedirs(full_out_dir, exist_ok=True)

    out_stub_short = "{}_{}".format(subject, session)
    out_stub = "sub-{}_ses-{}_parc-{}_conf-{}_spikereg-{}".format(subject, session, parc, conf, spikereg_threshold)

    conmat_file = os.path.join(full_out_dir, "{}_conmat.tsv".format(out_stub))
    conmat_file_feather = os.path.join(full_out_dir, "{}_conmat.feather".format(out_stub))
    conmat_plot_file = os.path.join(full_out_dir, "{}_conmat.png".format(out_stub))
    report_file = os.path.join(full_out_dir, "{}_report.txt".format(out_stub))
    outlier_stats_file = os.path.join(full_out_dir, "{}_outlier_stats.txt".format(out_stub))

    if not (os.path.exists(conmat_file) and os.path.exists(conmat_file_feather) and os.path.exists(conmat_plot_file)
            and os.path.exists(report_file) and os.path.exists(outlier_stats_file)):
        print("*** Calc conmats for {} {} {} {} {} ***".format(subject, session, parc, conf, spikereg_threshold))

        confounds_file, brainmask_file, rs_file, anat_file = get_files(fmriprep_dir, subject, session)
        roi_file, roi_names, roi_type = _get_roi_info(parc)

        conmat, report_str, outlier_stats = extract_mat(rs_file, brainmask_file, roi_file, confounds_file, conf,
                                                        roi_type, tr, spikereg_threshold)
        conmat_df = _get_con_df(conmat, roi_names)

        conmat_df.to_csv(conmat_file, sep="\t")
        with open(report_file, "w") as fi:
            fi.write(report_str)
        outlier_stats.to_csv(outlier_stats_file, sep="\t")

        plotting.plot_matrix(conmat, labels=roi_names, figure=(9, 7), vmax=1, vmin=-1, title=out_stub_short + " r")
        plt.savefig(conmat_plot_file, bbox_inches='tight')
        plt.close()

        save_feather(conmat_df, conmat_file_feather)

    else:
        print("*** Conmats for {} {} {} {} {} already computed. Do nothing. ***".format(subject, session, parc,
                                                                                        conf, spikereg_threshold))


def extract_mat(rs_file, brainmask_file, roi_file, confounds_file, conf, roi_type, tr, spikereg_threshold=None):
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
    confounds, outlier_stats = get_confounds(confounds_file, kind=conf, spikereg_threshold=spikereg_threshold)
    time_series = masker.fit_transform(rs_file, confounds=confounds.values)

    con_measure = connectome.ConnectivityMeasure(kind='correlation')
    conmat = con_measure.fit_transform([time_series])[0]

    report_str = "rs_file\t{}\n".format(rs_file)
    report_str += "roi_file\t{}\n".format(roi_file)

    keys = list(masker_pars.keys())
    keys.sort()
    for k in keys:
        report_str += "{}\t{}\n".format(k, masker_pars[k])

    report_str += "spike regression \t {}".format(spikereg_threshold)
    report_str += "\n\n"
    report_str += "confounds\t{}".format(", ".join(confounds.columns))
    report_str += "\n\n"
    report_str += confounds.to_string()

    return conmat, report_str, outlier_stats
