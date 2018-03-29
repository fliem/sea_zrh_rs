import os

import numpy as np
from nilearn import input_data, plotting
from glob import glob
from utils import get_files, get_36P_confounds


def sbc_one_session(subject, session, fmriprep_dir, output_dir, tr):

    confounds_file, brainmask_file, rs_file, anat_file = get_files(fmriprep_dir, subject, session)
    out_stub = "{}_{}".format(subject, session)

    # look for output files
    out_file_nii = os.path.join(output_dir, "sbc_pcc_1_" + out_stub + ".nii.gz")
    out_file_thresh = os.path.join(output_dir, "sbc_pcc_2_fisherz_thrsh0.5_" + out_stub + ".png")
    out_file_report = os.path.join(output_dir, "sbc_pcc_9_info_" + out_stub + ".txt")

    if not (os.path.exists(out_file_thresh) and os.path.exists(out_file_nii) and os.path.exists(out_file_report)):
        print("*** Running SBC for {} {} ***".format(subject, session))
        # see http://nilearn.github.io/auto_examples/03_connectivity/plot_seed_to_voxel_correlation.html
        pcc_coords = [(0, -52, 18)]
        lp_freq = 0.1
        hp_freq = 0.01

        # load confounds
        confounds = get_36P_confounds(confounds_file)

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
        seed_based_correlations_fisher_z = np.arctanh(seed_based_correlations)

        # plot
        seed_based_correlation_img = brain_masker.inverse_transform(seed_based_correlations_fisher_z.T)
        seed_based_correlation_img.to_filename(out_file_nii)

        display = plotting.plot_stat_map(seed_based_correlation_img, threshold=0.5, bg_img=anat_file,
                                         cut_coords=pcc_coords[0], title=out_stub + " (fisher z)")
        display.add_markers(marker_coords=pcc_coords, marker_color='g', marker_size=300);
        display.savefig(out_file_thresh)
        display.close()

        # report
        report = "Confounds:\n"
        report += "\n".join(confounds.columns) + "\n\n"
        report += "lp_freq: {}\n".format(lp_freq)
        report += "hp_freq: {}\n".format(hp_freq)

        report += "confounds_file: {}\n".format(confounds_file)
        report += "brainmask_file: {}\n".format(brainmask_file)
        report += "rs_file: {}\n".format(rs_file)
        report += "anat_file: {}\n\n".format(anat_file)

        report += confounds.to_string()

        with open(out_file_report, "w") as fi:
            fi.write(report)
    else:
        print("*** SBC for {} {} already computed. Do nothing. ***".format(subject, session))

def sbc_group(in_dir, out_dir):
    """
    for each session: load sbc (z transformed) of all subjects and calculate mean & sd
    """
    os.makedirs(out_dir, exist_ok=True)
    niis = glob(os.path.join(in_dir, "*.nii.gz"))
    sessions = list(set(map(lambda f: os.path.basename(f).split("_")[-1].split(".")[0], niis)))
    sessions.sort()
    print("calc mean sbc for {}".format(sessions))
    for ses in sessions:
        niis = glob(os.path.join(in_dir, "*{}.nii.gz".format(ses)))
        niis.sort()
        print(niis)
        from nilearn import image, plotting

        out_filename_mean_nii = os.path.join(out_dir, "sbc_1_mean_ses-{}.nii.gz".format(ses))
        out_filename_sd_nii = os.path.join(out_dir, "sbc_1_sd_ses-{}.nii.gz".format(ses))
        out_filename_mean_png = os.path.join(out_dir, "sbc_2_mean_ses-{}.png".format(ses))
        out_filename_sd_png = os.path.join(out_dir, "sbc_2_sd_ses-{}.png".format(ses))
        out_filename_list = os.path.join(out_dir, "sbc_ses-{}_scans.txt".format(ses))

        mean_sbc = image.mean_img(niis)
        mean_sbc.to_filename(out_filename_mean_nii)
        concat_img = image.concat_imgs(niis)
        sd_sbc = image.math_img("np.std(img, axis=-1)", img=concat_img)
        sd_sbc.to_filename(out_filename_sd_nii)

        with open(out_filename_list, "w") as fi:
            fi.write("\n".join(niis))

        display = plotting.plot_stat_map(mean_sbc, threshold=0.5,
                                         cut_coords=(0, -52, 18), title="mean sbc {} (fisher z)".format(ses))
        display.add_markers(marker_coords=[(0, -52, 18)], marker_color='g', marker_size=300);
        display.savefig(out_filename_mean_png)
        display.close()

        display = plotting.plot_stat_map(sd_sbc,
                                         cut_coords=(0, -52, 18), title="sd sbc {} (fisher z)".format(ses))
        display.add_markers(marker_coords=[(0, -52, 18)], marker_color='g', marker_size=300);
        display.savefig(out_filename_sd_png)
        display.close()