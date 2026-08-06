import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path
from astropy.io import fits

level1_folder = Path(r"D:\2026\FinalMaterials\ReProcessedData\WAVES\October2nd\1sExposures\Level1Data")
level2_folder = Path(r"D:\2026\FinalMaterials\ReProcessedData\WAVES\October2nd\1sExposures\Level2Data")
calibration_folder = Path(r"D:\2026\CalibrationData\release\calibration_data_20241205000000")
date='October2nd'
exposure_time='1'

level1_files=sorted(level1_folder.glob("*.fits*"))
level1_cube = np.memmap(rf'{level1_folder}\{date}_{exposure_time}_level1.dat',dtype=np.float32,mode='r',shape=(len(level1_files),2048,2048))
level2_files=sorted(level2_folder.glob("*.fits*"))
level2_cube = np.memmap(rf'{level2_folder}\{date}_{exposure_time}_level2.dat',dtype=np.float32,mode='r',shape=(len(level2_files),2048,2048))
shotnoise=np.memmap(rf'{level2_folder}\{date}_{exposure_time}_shotnoise.dat',dtype=np.float32,mode='w+',shape=(2048,2048))

Gain = 0.119    #DN/electron
readout = 7     #DN
N=1
Exptime=float(exposure_time)
Conversion=23600000000.0
biasA=fits.getdata(rf'{calibration_folder}\bias_A.fits')
biasB=fits.getdata(rf'{calibration_folder}\bias_B.fits')
temporary_file=level2_files[0]
temperature=fits.getheader(temporary_file)['APS_TEMP']
bias=biasA+biasB*temperature


badcolumns=[96,152,157,167,185,223,227,261,333,346,366,375,459,472,473,490,496,
           508,525,539,626,713,762,810,834,835,956,1024,1026,1038,1100,1102,1143,
           1155,1166,1209,1231,1257,1284,1301,1302,1416,1444,1490,1495,1504,1543,
           1602,1621,1626,1654,1818,1871,1885,1893,1897,1906,1955,1961]

for column in badcolumns:

    bias[:,column]=np.nan

contributing_data=np.nanmean(np.array(level1_cube),axis=0)-bias
ShotNoise_photon_dark = np.sqrt(contributing_data / Gain)
ShotNoise_photon_dark_DN = ShotNoise_photon_dark * Gain
ShotNoise_DN = np.sqrt(ShotNoise_photon_dark_DN**2+readout**2)
ShotNoise_MSB = (ShotNoise_DN / Exptime) * Conversion**-1

#testing old method
'''contributing_data=np.nanmean(np.array(level2_cube),axis=0)
ShotNoise_MSB = np.sqrt((contributing_data * Gain)/(Conversion*Exptime))
readout=(readout/Exptime)*Conversion**-1
ShotNoise_MSB = np.sqrt(ShotNoise_MSB**2+readout**2)'''

shotnoise[:]=ShotNoise_MSB
shotnoise.flush()

plt.style.use('dark_background')
plt.imshow(ShotNoise_MSB,cmap='viridis',origin='lower',vmin=np.nanpercentile(ShotNoise_MSB,1),vmax=np.nanpercentile(ShotNoise_MSB,99))
plt.title('Heat Map of Shot Noise')
colorbar=plt.colorbar()
colorbar.set_label('MSB')
plt.xlabel('X Pixel')
plt.ylabel('Y Pixel')
plt.show()


'''shells=np.linspace(450,1000,50)
ny, nx = DeltaIoverI.shape
y, x = np.indices((ny, nx))
xc = 1023.5
yc = 1023.5
x = x - xc
y = y - yc
r = np.sqrt(x**2 + y**2)
theta = np.arctan2(y, x)
fluctuationvalues=[]
intensityvalues=[]
shotnoisevalues=[]
measuredfluctuationvalues=[]

for j in range(len(shells[:-1])):

    inner = shells[j]
    outer = shells[j+1]
    mask = (r >= inner) & (r < outer)

    DeltaIoverI_in = DeltaIoverI[mask]
    mean_I_in = mean_I[mask]
    shotnoise_in = ShotNoise[mask]
    measuredfluctuations_in = DeltaI[mask]

    fluctuationvalues.append(np.nanmean(DeltaIoverI_in))
    intensityvalues.append(np.nanmean(mean_I_in))
    shotnoisevalues.append(np.nanmean(shotnoise_in))
    measuredfluctuationvalues.append(np.nanmean(measuredfluctuations_in))

gradient=np.gradient(intensityvalues,shells[:-1])
shiftingnoise=np.abs(gradient*(1/np.sqrt(2)))

#Noise Sources

shells=np.linspace(450,1000,49)
plt.style.use('default')
plt.figure(figsize=(6,4))
plt.xlabel('Pixel')
plt.ylabel(r'$\Delta I$')

plt.plot(shells,measuredfluctuationvalues,color='blue',label='Measured fluctuations')
plt.plot(shells,shiftingnoise,color='red',label='Misalignment noise',ls='--')
plt.plot(shells,shotnoisevalues,color='green',label='Shot noise')
plt.plot(shells,readoutnoise_MSB_array,color='orange',label='Read-out Noise')
plt.plot(shells,darkcurrent_noise_MSB_array,color='purple',label='Dark current noise')
plt.title('Measured Fluctuations and Noise Sources as a Function of Height')
plt.legend()

plt.show()'''