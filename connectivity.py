import os

import numpy as np
import pandas as pd
from nilearn import input_data, plotting

from utils import get_files


def sbc_one_subject(subject, session, fmriprep_dir, output_dir, tr, regressors):

    confounds_file, brainmask_file, rs_file, anat_file = get_files(fmriprep_dir, subject, session)
    out_stub = "{}_{}".format(subject, session)

    # look for output files
    out_file_thresh = os.path.join(output_dir, "sbc_1_thrsh0.3_" + out_stub + ".png")
    out_file_unthresh = os.path.join(output_dir, "sbc_2_unthrsh_" + out_stub + ".png")
    out_file_report = os.path.join(output_dir, "sbc_9_info_" + out_stub + ".txt")

    if not (os.path.exists(out_file_thresh) and os.path.exists(out_file_unthresh) and os.path.exists(out_file_report)):
        print("*** Running SBC for {} {} ***".format(subject, session))
        # see http://nilearn.github.io/auto_examples/03_connectivity/plot_seed_to_voxel_correlation.html
        pcc_coords = [(0, -52, 18)]
        lp_freq = 0.1
        hp_freq = 0.01

        # load confounds
        df = pd.read_csv(confounds_file, sep="\t")
        confounds = df[regressors].copy()
        #  replace the NaN at the beginning of FramewiseDisplacement with 0
        if 'FramewiseDisplacement' in regressors:
            confounds.loc[0, ['FramewiseDisplacement']] = 0

        # extract data from seed ROI
        seed_masker = input_data.NiftiSpheresMasker(
            pcc_coords, radius=8,
            mask_img=brainmask_file,
            detrend=True, standardize=True,
            low_pass=lp_freq, high_pass=hp_freq, t_r=tr)
        seed_time_series = seed_masker.fit_transform(rs_file, confounds=confounds.values)

        #  extract data for the entire brain
        brain_masker = input_data.NiftiMasker(
            smoothing_fwhm=6,
            mask_img=brainmask_file,
            detrend=True, standardize=True,
            low_pass=lp_freq, high_pass=hp_freq, t_r=tr)
        brain_time_series = brain_masker.fit_transform(rs_file, confounds=confounds.values)

        #  calculate correlation
        seed_based_correlations = np.dot(brain_time_series.T, seed_time_series) / \
                                  seed_time_series.shape[0]

        # plot
        seed_based_correlation_img = brain_masker.inverse_transform(seed_based_correlations.T)

        display = plotting.plot_stat_map(seed_based_correlation_img, threshold=0.3, bg_img=anat_file,
                                         cut_coords=pcc_coords[0], title=out_stub)
        display.add_markers(marker_coords=pcc_coords, marker_color='g', marker_size=300);
        display.savefig(out_file_thresh)
        display.close()

        display = plotting.plot_stat_map(seed_based_correlation_img, bg_img=anat_file,
                                         cut_coords=pcc_coords[0], title=out_stub)
        display.add_markers(marker_coords=pcc_coords, marker_color='g', marker_size=300);
        display.savefig(out_file_unthresh)
        display.close()

        # report
        report = "Confounds:\n"
        report += "\n".join(regressors) + "\n\n"
        report += "lp_freq: {}\n".format(lp_freq)
        report += "hp_freq: {}\n".format(hp_freq)

        report += "confounds_file: {}\n".format(confounds_file)
        report += "brainmask_file: {}\n".format(brainmask_file)
        report += "rs_file: {}\n".format(rs_file)
        report += "anat_file: {}\n".format(anat_file)

        with open(out_file_report, "w") as fi:
            fi.write(report)
    else:
        print("*** SBC for {} {} already computed. Do nothing. ***".format(subject, session))
