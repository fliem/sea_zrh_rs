from io import StringIO
import pandas as pd
from nilearn import image, datasets
import numpy as np
import os
import argparse


def split_yeo(out_path):
    # take yeo atlas and offset rois in rh +100
    atlas_yeo_2011 = datasets.fetch_atlas_yeo_2011()
    atlas_yeo = atlas_yeo_2011.thick_17
    atlas_yeo_img = image.reorder_img(image.image.check_niimg_3d(atlas_yeo))

    split_x = 128
    offset = np.zeros_like(atlas_yeo_img.dataobj)
    offset[split_x:, :, :] = 100
    offset[atlas_yeo_img.dataobj == 0] = 0
    offset_img = image.new_img_like(atlas_yeo_img, offset, copy_header=True, affine=atlas_yeo_img.affine)
    yeo_splithem_img = image.math_img("im1 + im2", im1=atlas_yeo_img, im2=offset_img)

    # create updated roi list
    yeo17_networks_str = StringIO("""roi,roi_name,dmn_subnetwork
    1,Visual A,
    2,Visual B,
    3,SomMot A,
    4,SomMot B,
    5,DorsAttn A,
    6,DorsAttn B,
    7,SalVentAttn A,
    8,SalVentAttn B,
    9,Limbic B,
    10,Limbic A,
    11,Control C,
    12,Control A,
    13,Control B,
    14,TempPar,
    15,Default C,mtl
    16,Default A,core
    17,Default B,dorsal medial""")
    df_labels = pd.read_csv(yeo17_networks_str, header=0)

    df_labels_lh = df_labels.copy()
    df_labels_rh = df_labels.copy()
    df_labels_lh["hemi"] = "lh"
    df_labels_rh["hemi"] = "rh"
    df_labels_rh["roi"] += 100

    df_combined = pd.concat((df_labels_lh, df_labels_rh))

    # output
    os.makedirs(out_path, exist_ok=True)
    out_file = os.path.join(out_path, "yeo_2011_thick17_splithemi")
    yeo_splithem_img.to_filename(out_file + ".nii.gz")

    df_combined.to_csv(out_file + ".tsv", sep="\t", index=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='splits rois of yeo atlas into one roi per hemi')
    parser.add_argument('out_path', help='path to save split nii files')
    args = parser.parse_args()

    split_yeo(args.out_path)
