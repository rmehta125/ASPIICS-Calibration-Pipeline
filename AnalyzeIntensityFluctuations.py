import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from pathlib import Path
from scipy.ndimage import shift
from astropy.wcs import WCS
from scipy.ndimage import rotate

#USER INPUTS _________________________________________________________________________________________________

level2_folder = Path(r"D:\2026\FinalMaterials\ReProcessedData\WAVES\October2nd\1sExposures\Level2Data")
level3_folder = Path(r"D:\2026\FinalMaterials\ReProcessedData\WAVES\October2nd\MergedExposures\Level3Data")
date='October2nd'
exposure_time='1'

#_____________________________________________________________________________________________________________

level2_files=sorted(level2_folder.glob("*.fits*"))
cube = np.array(np.memmap(rf"{level2_folder}\{date}_{exposure_time}_level2.dat",dtype=np.float32,mode='r',shape=(len(level2_files),2048,2048)))
ShotNoise = np.array(np.memmap(rf"{level2_folder}\{date}_{exposure_time}_shotnoise.dat",dtype=np.float32,mode='r',shape=(2048,2048)))

level3_files=sorted(level3_folder.glob("*.fits*"))
cube_level3 = np.array(np.memmap(rf'{level3_folder}\{date}_{exposure_time}_level3.dat',dtype=np.float32,mode='r',shape=(len(level3_files),2048,2048)))


'''#NORMALIZE BY INTENSITY FROM LEVEL 2
mean_I = np.nanmean(cube,axis=0)
DeltaI = np.nanstd(cube, axis=0)
DeltaI_Calibrated = np.sqrt(DeltaI**2-ShotNoise**2)
DeltaIoverI = DeltaI_Calibrated / mean_I'''

#NORMALIZE BY INTENSITY FROM LEVEL 3
file0=level2_files[0]
header=fits.getheader(file0)
wcs = WCS(header)
sun_x, sun_y = wcs.world_to_pixel_values(0.0, 0.0)
xshift=sun_x-1023.5
yshift=sun_y-1023.5

mean_I=np.nanmean(cube_level3,axis=0)
rotation_angle=fits.getheader(file0)['CROTA']
mean_I=rotate(mean_I,angle=rotation_angle,reshape=False,cval=np.nan,order=0)
#mean_I=shift(mean_I,shift=(0,-20),order=0,mode='constant',cval=np.nan)
DeltaI = np.nanstd(cube, axis=0)
DeltaI_Calibrated = np.sqrt(DeltaI**2-ShotNoise**2)
DeltaI_Calibrated = shift(DeltaI_Calibrated,shift=(-yshift-10,-xshift-10),order=0,mode='constant',cval=np.nan)
#DeltaI_Calibrated = rotate(DeltaI_Calibrated,angle=-rotation_angle,reshape=False,cval=np.nan,order=0)
DeltaIoverI = DeltaI_Calibrated / mean_I

plt.figure(figsize=(7,7))
plt.imshow(mean_I,origin='lower',cmap='magma')
plt.imshow(DeltaI_Calibrated,origin='lower',cmap='viridis')
plt.show()


#EDIT THIS SECTION ONCE YOU FIGURE OUT THE SHIFTING

file0=level2_files[0]
image_data=fits.getdata(file0)
header=fits.getheader(file0)
observation_date=header['DATE-OBS']
wcs = WCS(header)
sun_x, sun_y = wcs.world_to_pixel_values(0.0, 0.0)
xshift=sun_x-1023.5
yshift=sun_y-1023.5
CDELT1=header['CDELT1']
solarradius=header['RSUN_ARC']
#DeltaIoverI=shift(DeltaIoverI,shift=(-yshift,-xshift),order=0,mode='constant',cval=np.nan)

smallestshell=(1.25*solarradius)/CDELT1
shells=np.linspace(smallestshell,1000,61)
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

    #fluctuationvalues.append(np.nanmean(DeltaIoverI_in))
    fluctuationvalues.append(np.nanmedian(DeltaIoverI_in))

shells=np.linspace(int(smallestshell),1000,60)
shells=(shells*CDELT1)/solarradius
plt.figure(figsize=(6,4))
plt.xlabel('Solar Radii')
plt.ylabel(r'$\Delta$ I / I')
plt.plot(shells,fluctuationvalues,color='blue')
plt.title('Intensity Fluctuations as a Function of Height')

plt.style.use('dark_background')
fig,ax=plt.subplots(figsize=(8,8))
ax.set_facecolor(color="#010129")
im = plt.imshow(DeltaIoverI,cmap='afmhot',origin='lower',vmin=0,vmax=0.5)
colorbar=plt.colorbar()
colorbar.set_label(r'$\Delta$ I/I')
plt.title(f'Intensity Fluctuation Heat Map, {str(observation_date)[:10]}')
plt.xlabel('X Pixel')
plt.ylabel('Y Pixel')

fluctuations=np.memmap(rf'{level2_folder}\{date}_{exposure_time}_fluctuations.dat',dtype=np.float32,mode='w+',shape=(2048,2048))
fluctuations[:]=DeltaIoverI
fluctuations.flush()

#print(np.array(fluctuationvalues).tolist())
print(np.array(shells).tolist())


plt.show()

