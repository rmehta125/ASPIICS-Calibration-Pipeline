import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path
from astropy.io import fits
from scipy.ndimage import rotate
from astropy.wcs import WCS
from matplotlib.colors import ListedColormap

# USER INPUTS ____________________________________________________________________________________________________

level2_folder = Path(r"D:\2026\FinalMaterials\ReProcessedData\WAVES\October2nd\1sExposures\Level2Data")
date='October2nd'
exposure_time='1'

# ________________________________________________________________________________________________________________

level2_files=sorted(level2_folder.glob("*.fits*"))
fluctuations = np.array(np.memmap(rf"{level2_folder}\{date}_{exposure_time}_fluctuations.dat",dtype=np.float32,mode='r',shape=(2048,2048)))

def rotate_point(origin, point, angle):

    origin_x, origin_y = origin
    point_x, point_y = point
    angle = np.radians(angle)

    new_point_x = origin_x + np.cos(angle) * (point_x - origin_x) - np.sin(angle) * (point_y - origin_y)
    new_point_y = origin_y + np.sin(angle) * (point_x - origin_x) + np.cos(angle) * (point_y - origin_y)

    return new_point_x, new_point_y

sun_radii=[]
rotation_angles=[]
occulter_x_positions=[]
occulter_y_positions=[]

for file in level2_files:

    header=fits.getheader(file)
    sun_radii.append(header['RSUN_ARC'])
    rotation_angles.append(header['CROTA'])
    occulter_x_positions.append(header['X_IO'])
    occulter_y_positions.append(header['Y_IO'])
    CDELT1=header['CDELT1']

sun_radii=np.array(sun_radii)
rotation_angles=np.array(rotation_angles)
occulter_x_positions=np.array(occulter_x_positions)
occulter_y_positions=np.array(occulter_y_positions)
RSUN_ARC=np.mean(sun_radii)
rotation_angle=np.mean(rotation_angles)

fluctuations=rotate(fluctuations,angle=-rotation_angle,reshape=False,cval=np.nan,order=0)

file0=level2_files[0]
image_data=fits.getdata(file0)
header=fits.getheader(file0)
observation_date=header['DATE-OBS']
wcs = WCS(header)
sun_x, sun_y = wcs.world_to_pixel_values(0.0, 0.0)
xshift=sun_x-1023.5
yshift=sun_y-1023.5
occulter_x=np.mean(occulter_x_positions)-xshift
occulter_y=np.mean(occulter_y_positions)-yshift
occulter_x,occulter_y=rotate_point((1023.5,1023.5),(occulter_x,occulter_y),rotation_angle)

plt.style.use('dark_background')
fig,ax=plt.subplots(figsize=(8,7))
im=ax.imshow(fluctuations,cmap='magma',origin='lower',vmin=0,vmax=0.3)






colorbar=fig.colorbar(im,label=r'$\Delta I / I$',shrink=0.9)
occulter_radius=390
occulter = plt.Circle((occulter_x,occulter_y), radius=occulter_radius, fill=True, color='black')
occulter_outline = plt.Circle((occulter_x,occulter_y), radius=occulter_radius, fill=False, color='lightsteelblue', linewidth=0.5)
sun_circle = plt.Circle((1023.5, 1023.5), radius=(RSUN_ARC)/CDELT1, color='lightsteelblue',fill=False,linewidth=0.5,ls=(0,(5,10)))
ax.add_patch(sun_circle)

theta_ticks=[0,30,60,120,150]
ticks=np.linspace(0,2047,2048)

for tick in theta_ticks:

    ax.plot(ticks,1023.5+(ticks-1023.5)*np.tan(np.radians(tick)),color='lightsteelblue',linewidth=0.5)

ax.set_xlim(0,2048)
ax.set_ylim(0,2048)
ax.vlines(1023.5,0,2048,color='lightsteelblue',linewidth=0.5)
ax.set_title(rf'Intensity Fluctuation Heat Map, {str(observation_date)[:10]}')

plt.show()

angle=float(int(input('Please enter the angle for lineout analysis, measured counterclockwise from up.  Enter as an integer: ')))
angle+=90

fluctuations=rotate(fluctuations,angle=angle,reshape=False,cval=np.nan,order=0)
occulter_x,occulter_y=rotate_point((1023.5,1023.5),(occulter_x,occulter_y),-angle)



ny, nx = fluctuations.shape
y, x = np.indices((ny, nx))
xc = 1023.5
yc = 1023.5
x = x - xc
y = y - yc
r = np.sqrt(x**2 + y**2)
theta = np.arctan2(y, x)
theta[theta<0]+=2*np.pi
theta=np.degrees(theta)

mask = (x > 0) & ((y > -2) & (y < 2))

totalvalues=fluctuations[mask]
cutoff=np.nanpercentile(totalvalues,99.5)
print(cutoff)

radialbins=np.arange(350,1030,5)
radialmidpoints=[]
meanfluctuations=[]
errors=[]

for i,radius in enumerate(radialbins[:-1]):


    endradius=radialbins[i+1]
    midpoint=(radius+endradius)/2
    values=fluctuations[(radius < r) & (r <= endradius) & mask & (fluctuations < cutoff)]

    value=np.nanmean(values)
    error=np.nanstd(values)

    errors.append(error)
    radialmidpoints.append(midpoint)
    meanfluctuations.append(value)

radialmidpoints=(np.array(radialmidpoints)*CDELT1)/RSUN_ARC
meanfluctuations=np.array(meanfluctuations)
errors=np.array(errors)

plt.style.use('dark_background')
fig1,ax1=plt.subplots(figsize=(9,4))
plt.subplots_adjust(bottom=0.15)
fig2,ax2=plt.subplots(figsize=(8,7))

ax1.plot(radialmidpoints,meanfluctuations,color='#ff9408')
ax1.fill_between(radialmidpoints,meanfluctuations-0.5*errors,meanfluctuations+0.5*errors,alpha=0.3,color="#ff9408")
ax1.set_xlabel(r'Solar Radii')
ax1.set_ylabel(r'$\Delta$/I')
ax1.set_title(rf'$\Delta I/I$ as a function of $r$ at $\phi={int(angle-90)}^\circ$, {str(observation_date)[:10]}')
ax1.set_facecolor(color="#010129")
ax1.set_xlim(1.25,3)
ax1.set_ylim(0,np.nanmax(meanfluctuations)+0.15)

im=ax2.imshow(fluctuations,origin='lower',vmin=0,vmax=0.3,cmap=new_cmap)
colorbar=fig2.colorbar(im,label=r'$\Delta I/I$',shrink=0.9)
occulter_radius=390
occulter = plt.Circle((occulter_x,occulter_y), radius=occulter_radius, fill=True, color='black')
occulter_outline = plt.Circle((occulter_x,occulter_y), radius=occulter_radius, fill=False, color='lightsteelblue', linewidth=0.5)
sun_circle = plt.Circle((1023.5, 1023.5), radius=(RSUN_ARC)/CDELT1, color='lightsteelblue',fill=False,linewidth=0.5,ls=(0,(5,10)))
ax2.add_patch(sun_circle)
ax2.set_xlim(0,2048)
ax2.set_ylim(0,2048)
ax2.set_xticks([])
ax2.set_yticks([])
ax2.set_facecolor(color="#010129")
ax2.set_title(rf'Intensity Fluctuation Heat Map, {str(observation_date)[:10]}')
ax2.hlines(1023.5,1023.5,2048,color='white',linewidth=0.5)

print(f'radialmidpoints={radialmidpoints[~np.isnan(meanfluctuations)].tolist()}')
print(f'meanfluctuations={meanfluctuations[~np.isnan(meanfluctuations)].tolist()}')
print(f'errors={errors[~np.isnan(meanfluctuations)].tolist()}')

plt.show()


