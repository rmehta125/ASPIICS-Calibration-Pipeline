import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from pathlib import Path
from scipy.ndimage import shift

level2_folder = Path(r"D:\2026\FinalMaterials\ReProcessedData\CMEWATCH\July18th\10sExposures\Level2Data")
date='July18th'
exposure_time='10'


level2_files=sorted(level2_folder.glob("*.fits*"))
cube = np.array(np.memmap(rf"{level2_folder}\{date}_{exposure_time}_level2.dat",dtype=np.float32,mode='r',shape=(len(level2_files),2048,2048)))
ShotNoise = np.array(np.memmap(rf"{level2_folder}\{date}_{exposure_time}_shotnoise.dat",dtype=np.float32,mode='r',shape=(2048,2048)))

mean_I = np.nanmean(cube,axis=0)
DeltaI = np.nanstd(cube, axis=0)
DeltaI_Calibrated = np.sqrt(DeltaI**2-ShotNoise**2)
DeltaIoverI = DeltaI_Calibrated / mean_I

shells=np.linspace(450,1000,50)
ny, nx = DeltaIoverI.shape
y, x = np.indices((ny, nx))
xc = 1023.5
yc = 1023.5
x = x - xc
y = y - yc
r = np.sqrt(x**2 + y**2)
theta = np.arctan2(y, x)
fluctuationvalues=[]

for j in range(len(shells[:-1])):

    inner = shells[j]
    outer = shells[j+1]
    mask = (r >= inner) & (r < outer)

    DeltaIoverI_in = DeltaIoverI[mask]
    mean_I_in = mean_I[mask]
    shotnoise_in = ShotNoise[mask]
    measuredfluctuations_in = DeltaI[mask]

    fluctuationvalues.append(np.nanmean(DeltaIoverI_in))

CDELT1=2.8172
RSUN_ARC=958.669244220219
shells=np.linspace(450,1000,49)
shells=(shells*CDELT1)/RSUN_ARC
plt.figure(figsize=(6,4))
plt.xlabel('Solar Radii')
plt.ylabel(r'$\Delta$ I / I')
plt.plot(shells,fluctuationvalues,color='blue')
plt.title('Intensity Fluctuations as a Function of Height')

plt.style.use('dark_background')
fig,ax=plt.subplots(figsize=(8,8))
im = plt.imshow(DeltaIoverI,cmap='viridis',origin='lower',vmin=0,vmax=0.08)
colorbar=plt.colorbar()
colorbar.set_label(r'$\Delta$ I/I')
plt.title('Intensity Fluctuation Heat Map')
plt.xlabel('X Pixel')
plt.ylabel('Y Pixel')

fluctuations=np.memmap(rf'{level2_folder}\{date}_{exposure_time}_fluctuations.dat',dtype=np.float32,mode='w+',shape=(2048,2048))
fluctuations[:]=DeltaIoverI
fluctuations.flush()

plt.show()