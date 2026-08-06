import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from pathlib import Path
from scipy.ndimage import shift
from scipy.ndimage import rotate
from skimage.registration import phase_cross_correlation
import cv2
from skimage.transform import warp_polar, rotate, rescale
from skimage.util import img_as_float

#USER INPUTS _____________________________________________________________________________________________

level1_folder = Path(r"D:\2026\FinalMaterials\ReProcessedData\WAVES\October2nd\1sExposures\Level1Data")
level2_folder = Path(r"D:\2026\FinalMaterials\ReProcessedData\WAVES\October2nd\1sExposures\Level2Data")
date='October2nd'
exposure_time='1'
xshifts=[0, 0, -1, -1, -1, 0, 0, -1, 0, 0, -1, 0, -1, -1, -1, -1, -1, -1, 0, -1, -1, -1, -1, -1, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 0, -1, 0, 0, -1, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 0, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 0, -1, -1, -1, 0, -1, -1, -1, -1, -1, 0, 0, -1, 0, -1, -1, -1, -1, -1, -1, -1, -1, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1]
yshifts=[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 0, 0, 0, 1, 0, 0, 1, 0, 1, 1, 0, 0, 0, 1, 0, 1, 1, 1, 1, 1, 0, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 1, 1, 1]

# _______________________________________________________________________________________________________


level1_files=sorted(level1_folder.glob("*.fits*"))
level1_cube = np.memmap(rf'{level1_folder}\{date}_{exposure_time}_level1.dat',dtype=np.float32,mode='w+',shape=(len(level1_files),2048,2048))
level2_files=sorted(level2_folder.glob("*.fits*"))
level2_cube = np.memmap(rf'{level2_folder}\{date}_{exposure_time}_level2.dat',dtype=np.float32,mode='w+',shape=(len(level2_files),2048,2048))

badcolumns=[96,152,157,167,185,223,227,261,333,346,366,375,459,472,473,490,496,
           508,525,539,626,713,762,810,834,835,956,1024,1026,1038,1100,1102,1143,
           1155,1166,1209,1231,1257,1284,1301,1302,1416,1444,1490,1495,1504,1543,
           1602,1621,1626,1654,1818,1871,1885,1893,1897,1906,1955,1961]

for i, file in enumerate(level1_files):

    print(f"Shifting level 1 image {i}")

    image_data=fits.getdata(file)
    image_data[~np.isfinite(image_data)]=np.nan
    image_data[image_data<0]=np.nan

    for column in badcolumns:

        image_data[:,column]=np.nan

    image_data = shift(image_data,shift=(-yshifts[i],-xshifts[i]),order=0,mode='constant',cval=0)
    level1_cube[i]=image_data

level1_cube.flush()


for i, file in enumerate(level2_files):

    print(f"Shifting level 2 image {i}")

    image_data=fits.getdata(file)
    image_data[~np.isfinite(image_data)]=np.nan
    image_data[image_data<0]=np.nan

    for column in badcolumns:

        image_data[:,column]=np.nan

    image_data = shift(image_data,shift=(yshifts[i],xshifts[i]),order=0,mode='constant',cval=0)
    level2_cube[i]=image_data

level2_cube.flush()
