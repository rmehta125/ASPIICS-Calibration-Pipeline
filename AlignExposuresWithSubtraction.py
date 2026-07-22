import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from scipy.ndimage import shift
from pathlib import Path

shift_direction='yshift'

folder = Path(r"D:\2026\FinalMaterials\ReProcessedData\CMEWATCH\July18th\1sExposures\Level2Data")
files = sorted(folder.glob("*.fits*"))

image0=fits.getdata(files[0])
image0[~np.isfinite(image0)]=np.nan
image0[image0<0]=np.nan

shifts=[]
testshifts=[-3,-2,-1,0,1,2,3]

for i,file in enumerate(files):

    image=fits.getdata(file)
    image[~np.isfinite(image)]=np.nan
    image[image<0]=np.nan

    totaldifferences=[]

    for s in testshifts:

        if shift_direction == 'xshift':

            image_shifted=shift(image,shift=(0,s),order=1,mode='constant',cval=0)

        elif shift_direction == 'yshift':

            image_shifted=shift(image,shift=(s,0),order=1,mode='constant',cval=0)

        totaldifference=0

        region=image_shifted[520:555,870:1020]
        region0=image0[520:555,870:1020]
        difference=np.nansum(np.abs(region-region0))
        totaldifference+=difference

        region=image_shifted[868:1168,1417:1442]
        region0=image0[868:1168,1417:1442]
        difference=np.nansum(np.abs(region-region0))
        totaldifference+=difference

        region=image_shifted[925:1177,542:580]
        region0=image0[925:1177,542:580]
        difference=np.nansum(np.abs(region-region0))
        totaldifference+=difference

        region=image_shifted[1450:1525,850:1130]
        region0=image0[1450:1525,850:1130]
        difference=np.nansum(np.abs(region-region0))
        totaldifference+=difference

        totaldifferences.append(totaldifference)

    idx=totaldifferences.index(min(totaldifferences))
    ideal_shift=testshifts[idx]
    shifts.append(ideal_shift)
    print(f"Ideal shift for image {i}: {ideal_shift}")

if shift_direction=='xshift':

    print(f'xshifts={shifts}')

if shift_direction=='yshift':

    print(f'yshifts={shifts}')
