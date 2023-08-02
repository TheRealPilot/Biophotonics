"""
Decomposition of glucose solution spectral data into constituent spectra.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from math import pi

#Use these ratios if plotting 0.0555 M spectra
GLUCOSE_RATIO_1 = 0.00625
WATER_RATIO_1 = 0.99375

#Use these ratios if plotting 2 M spectra
GLUCOSE_RATIO_2 = 0.26488
WATER_RATIO_2 = 0.73512

#Use these ratios if plotting saturated glucose
GLUCOSE_SAT = 0.75
WATER_SAT = 0.25

#Load the spectral datasets

#Load Saturated spectral datasets
data2 = pd.read_csv('A:\PythonStuff\ActuallyPython\BIophotonics\Data\Saturated.csv')
df2 = data2.drop(['Unnamed: 6'], axis=1)
wavelengths = df2['empty'][1:][::-1].values.astype(float)
water2 = df2['Unnamed: 3'][1:][::-1].values.astype(float)
glucoseSat = df2['Unnamed: 5'][1:][::-1].values.astype(float)

dataf = pd.read_csv("A:\PythonStuff\ActuallyPython\BIophotonics\Data\pathlength-test.csv")
wavelengthsf = dataf['Wavelength (nm)'][1:][::-1].values.astype(float)
waterf = 

# Load ham and optica data using pandas and sort by ascending wavelengths
data3 = pd.read_csv('A:\PythonStuff\ActuallyPython\BIophotonics\Data\HamGreat.csv')
data3_sorted = data3.sort_values(by='x')
hamWavelengths = data3_sorted['x'].values.astype(float)
hamAbsobance = data3_sorted['y'].values.astype(float)

data4 = pd.read_csv('A:\PythonStuff\ActuallyPython\BIophotonics\Data\OpticaGreat.csv')
data4_sorted = data4.sort_values(by='x')
opticaWavelengths = data4_sorted['x'].values.astype(float)
opticaAbsorbance = data4_sorted['y'].values.astype(float)

def isolate(water_r, glucose_r, absorbance_g, absorbance_w, path_length, molarity, transmittance=False):
    """
    Function demixes the spectral components of a mixture.  

    Input variables:

        water_r: mass ratio of water in glucose solution
        glucose_r: mass ratio of glucose in glucose solution
        epsilon_w: set of extinction coefficients for water
        epsilon_g: set of extinction coeffciencts for glucose
        path_length: the pathlength of container (cm)
        molarity: the molarity of glucose in the solution

    Returns: 

        Array of separated glucose spectra. 
    """
    epsilon_w, epsilon_g = extinction(WATER_RATIO_2, GLUCOSE_RATIO_2, absorbance_w, absorbance_g, 0.02, 2, transmittance=transmittance)
    ag = (water_r*epsilon_w + glucose_r*epsilon_g)*path_length*molarity 
    water_component = (water_r)*epsilon_w*path_length*molarity
    return np.abs(ag - water_component)

def division(water, glucose):

    return norm(glucose/water)

def extinction(ratio_w, ratio_g, absorbance_w, absorbance_g, path_length, concentration, transmittance=False):
    """
    Creates a set of molar absorptivity coefficients for water and glucose.

    Input variables:
        ratio_w: mass ratio of water in glucose solution
        ratio_g: mass reatio of glucose in glucose solution
        absorbance_w: absorbance spectra of pure DI
        absorbance_g: absorbance spectra of glucose solution
        path_length: the pathlength of the container (cm)
        concentration: the molarity of glucose in the solution (mol/liter)

    Returns:
        Molar absorbtivity coefficients for water and glucose (L/mol*cm) from 
        350 nm to 2500 nm.
    """
    absorbance_w, absorbance_g = norm(absorbance_w, absorbance_g)
    if transmittance:
        absorbance_w = transform(absorbance_w)
        absorbance_g = transform(absorbance_g)
    epsilon_w = absorbance_w / (path_length)
    epsilon_g = (absorbance_g/(path_length*concentration) - ratio_w*epsilon_w) / ratio_g
    return epsilon_w, epsilon_g

def CDtoORD(wavelengths, spectra):
    #NOTE: Not Finished. Likely obsolete. 
    """
    Function converts CD spectra from either delta extinction or degree of ellipticity to ORD. 
    This implementation follows the numerical KK-transform described by Polavarapu, which was 
    rewritten from Ohta and Ishida's method. 
    doi: 10.1021/jp0524328

    Innput variables:
        wavelengths: set of wavelengths of spectra
        spectra: spectral data

    Returns:
        Molar rotation data. 
    """
    try:
        units = input("Enter a number to select the type of CD Spectra:\n 1. Extinction Difference\n 2. Degree of Ellipticity")
        if not 0 < units < 3:
            raise ValueError("Invalid option. Please enter a valid number.")
    except ValueError as error:
        print(error)

    if units == 1:
        pass

def simga(wavelengths, intensities, h=0.01):
    """
    Implementation of the weird sigma operator described in Polavarapu's paper. 

    If the wavelength λ, where molar rotation is to be calculated, corresponds 
    to an odd data number, then the summation is carried over even data numbers.

    If the wavelength where molar rotation is to be calculated corresponds to an 
    even data number, then the summation is carried over odd data numbers.

    μj is wavelength, and [θ(μj)] is spectral intensity.
    In other words, μj is x and [θ(μj)] is y. 

    You need to calculate the entire sigma for each wavelength. 
    """
    spectra_even = []
    spectra_odd = []
    full_spectra = []

    for i in range(len(wavelengths)):
        if i % 2 != 0:
            sum = 0
            for j in wavelengths[1::2]:
                sum += (intensities[1::2][j]/(wavelengths[i] - wavelengths[1::2][j])) - (intensities[1::2][j]/(wavelengths[i] + wavelengths[1::2][j]))
            spectra_even += (2/pi)*(2*h)*(1/2)*sum
        if j % 2 == 0:
            sum = 0
            for j in wavelengths[::2]:
                sum += (intensities[::2][j]/(wavelengths[i] - wavelengths[::2][j])) - (intensities[::2][j]/(wavelengths[i] + wavelengths[::2][j]))
            spectra_odd += (2/pi)*(2*h)*(1/2)*sum
                
    for num1, num2 in zip(spectra_even, spectra_odd):
        full_spectra.append(num1)
        full_spectra.append(num2)

    if len(spectra_even) > len(spectra_odd):
        full_spectra.append(spectra_even[-1])
    elif len(spectra_odd) > len(spectra_even):
        full_spectra.append(spectra_odd[-1])

    return full_spectra

def plot(wavelengths, separated, water, glucose, pulled_wave, pulled, optica_wave, optica, nothing=None):
    """
    Plots spectra.
    """
    fig, ax = plt.subplots()
    optica = transform(norm(optica))
    separated = norm(separated)
    glucose = norm(glucose)
    water = norm(water)

    # Cut out wavelengths larger than 2400
    mask = wavelengths <= 2400
    wavelengths_cut = wavelengths[mask]
    separated_cut = separated[mask]
    water_cut = water[mask]
    glucose_cut = glucose[mask]

    if nothing is not None:
        #nothing, separated, glucose, water = norm(nothing, separated, glucose, water)
        ax.plot(wavelengths, nothing, label='No Sample')
        ax.plot(wavelengths_cut, separated_cut, label='Glucose Spectra')
        ax.plot(wavelengths_cut, glucose_cut, label='Solution Spectra')
        ax.plot(wavelengths_cut, water_cut, label='DI Spectra')
        ax.plot(list(pulled_wave), pulled, label='Extracted Data')  # Convert pulled_wave to a list
    else:
        #separated, glucose, water = norm(separated, glucose, water)
        ax.plot(wavelengths_cut, separated_cut, label='Glucose Spectra')
        #ax.plot(wavelengths_cut, glucose_cut, label='Solution Spectra')
        #ax.plot(wavelengths_cut, water_cut, label='DI Spectra')
        ax.plot(list(pulled_wave), pulled, label='Ham Data') # Convert pulled_wave to a list
        ax.plot(list(optica_wave), optica, label='Optica Data')      

    ax.set_xlabel('Wavelength (nm)')
    ax.set_ylabel('Normalized Absorbance')
    ax.set_title('Saturated Glucose Spectra (0.133 s integration time)')
    ax.legend()
    #ax.set_yscale('log')
    plt.show()


def norm(*args):
    """
    Normalizes a numpy array vector or a tuple of numpy array vectors. 
    Returns a single ndarray if only one vector is provided, and returns a tuple of ndarrays otherwise.
    """
    if len(args) == 1:
        vector = args[0]
        max_val = np.max(vector)
        min_val = np.min(vector)
        normalized_vector = (vector - min_val) / (max_val - min_val)
        return normalized_vector
    else:
        vectors = args
        concatenated = np.concatenate(vectors)
        max_val = np.max(concatenated)
        min_val = np.min(concatenated)
        normalized_vectors = [(vector - min_val) / (max_val - min_val) for vector in vectors]
        return tuple(normalized_vectors)


def transform(spectra):
    """Transfroms % transmittance to absorbance. Absorbance = 2 – log(%T)"""
    spectra = np.clip(spectra, a_min=1e-10, a_max=None)
    #return 1 - np.log10(spectra)
    return 1 - spectra

if __name__=="__main__":
    #print(data3_sorted)
    separated_spectra_10M = isolate(WATER_SAT, GLUCOSE_SAT, glucoseSat, water2, 0.02, 4.163)
    plot(wavelengths, transform(separated_spectra_10M), transform(water2), transform(glucoseSat), hamWavelengths, transform(hamAbsobance), opticaWavelengths, opticaAbsorbance)
