import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from pathlib import Path
plt.style.use('dark_background')

level2_folder = Path(r"D:\2026\FinalMaterials\ReProcessedData\WAVES\October8th\1sExposures\Level2Data")
date='October8th'
cadence=5
exposure_time=1
pixel = [1724,1531]

level2_files=sorted(level2_folder.glob("*.fits*"))
cube=np.array(np.memmap(rf'{level2_folder}\{date}_{exposure_time}_level2.dat',dtype=np.float32,mode='r',shape=(len(level2_files),2048,2048)))
shotnoise=np.array(np.memmap(rf'{level2_folder}\{date}_{exposure_time}_shotnoise.dat',dtype=np.float32,mode='r',shape=(2048,2048)))
fluctuations=np.array(np.memmap(rf'{level2_folder}\{date}_{exposure_time}_fluctuations.dat',dtype=np.float32,mode='r',shape=(2048,2048)))

fig2,ax2=plt.subplots(figsize=(5,5))
ax2.imshow(fluctuations,origin='lower',cmap='viridis',vmin=np.nanpercentile(fluctuations,1),vmax=0.08)
fig,ax=plt.subplots(figsize=(13,3))
plt.subplots_adjust(bottom=0.2)

#pixels = [[1066,1879],[123,851],[1571,1669],[1328,274]]
pixels = [pixel]       
masterpowers=np.zeros(len(level2_files))

for k,pixel in enumerate(pixels):

    pixel1series=[]
    pixel2series=[]
    pixel3series=[]
    pixel4series=[]
    pixel5series=[]
    pixel6series=[]
    pixel7series=[]
    pixel8series=[]
    pixel9series=[]

    image=np.nanmean(cube,axis=0)

    for i in range(len(cube[:,0,0])):

        data=cube[i,:,:]

        pixel1series.append(data[pixel[1],pixel[0]])
        pixel2series.append(data[pixel[1],pixel[0]+1])
        pixel3series.append(data[pixel[1]+1,pixel[0]+1])
        pixel4series.append(data[pixel[1]+1,pixel[0]])
        pixel5series.append(data[pixel[1]+1,pixel[0]-1])
        pixel6series.append(data[pixel[1],pixel[0]-1])
        pixel7series.append(data[pixel[1]-1,pixel[0]-1])
        pixel8series.append(data[pixel[1]-1,pixel[0]])
        pixel9series.append(data[pixel[1]-1,pixel[0]+1])

    pixelseries=np.array([pixel1series,pixel2series,pixel3series,pixel4series,pixel5series,
                pixel6series,pixel7series,pixel8series,pixel9series])

    powers=np.zeros(len(pixelseries[0]))

    for series in pixelseries:

        series=np.array(series)
        series=series-np.mean(series)
        number = np.sum(np.isnan(series))
        power=(np.abs(np.fft.fft(series)))**2
        powers+=power
        freq=np.fft.fftfreq(len(series),d=cadence)
        #ax.plot(freq[freq>0],power[freq>0],alpha=0.1)

    powers=powers/len(pixelseries)
    ax.plot(freq[freq>0],powers[freq>0],label='5s Cadence',alpha=1)
    ax.set_title(f'Power Spectrum of Pixel Intensity, {pixel}')
    ax.set_xlabel('Frequency (Hz)')
    ax.set_ylabel('Power')
    masterpowers+=powers

    ax2.plot(pixel[0],pixel[1],color='crimson',marker='o',alpha=1)

masterpowers=masterpowers/len(pixels)
#ax.plot(freq[freq>0],masterpowers[freq>0],label='5s Cadence',alpha=1,color='greenyellow',linewidth=2)

plt.show()