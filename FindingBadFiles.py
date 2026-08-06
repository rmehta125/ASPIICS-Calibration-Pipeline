import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from pathlib import Path

# USER INPUTS ___________________________________________

folder = Path(r"D:\2026\FinalMaterials\ReProcessedData\WAVES\October2nd\MergedExposures\Level3Data")

# _______________________________________________________

files = sorted(folder.glob("*.fits*"))

for i,file in enumerate(files):

    image = fits.getdata(file)
    number = np.sum(np.isnan(image))
    header = fits.getheader(file)
    print(header['DATE-OBS'],number)
