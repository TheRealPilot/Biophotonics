"""
Decomposition of glucose solution spectral data into constituent spectra.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from math import pi

#Use these ratios if plotting 179 mg/dL (0.00993561 M) spectra
GLUCOSE_RATIO_179 = 0.000446
WATER_RATIO_179 = 0.99955

#Use these ratios if plotting 0.0555 M spectra
GLUCOSE_RATIO_2500 = 0.00625
WATER_RATIO_2500 = 0.99375

#Use these ratios if plotting 2 M spectra
GLUCOSE_RATIO_2 = 0.26488
WATER_RATIO_2 = 0.73512

#Use these ratios if plotting saturated glucose
GLUCOSE_SAT = 0.75
WATER_SAT = 0.25

#Load the spectral datasets

data = pd.read_csv('/Users/davidl./Biophotonics/Data/200path.csv')
df = data.drop(['Unnamed: 4'], axis=1)
wavelengths2500 = df['sample1'][1:][::-1].values.astype(int)
water2500 = df['Unnamed: 1'][1:][::-1].values.astype(float)
glucose2500 = df['Unnamed: 3'][1:][::-1].values.astype(float)

dataf = pd.read_csv('/Users/davidl./Biophotonics/Data/2M Glucose.csv')
dff = dataf.drop(['Unnamed: 4'], axis=1)
wavelengths = dff['sample1'][1:][::-1].values.astype(int)
water2M = dff['Unnamed: 1'][1:][::-1].values.astype(float)
#solution179 = dff['Unnamed: 3'][1:][::-1].values.astype(float)
glucose2M = dff['Unnamed: 3'][1:][::-1].values.astype(float)

#4.16M not baselined RAW
data2 = pd.read_csv('//Users/davidl./Biophotonics/Data/Saturated.csv')
df2 = data2.drop(['empty'], axis=1)
wavelengthsSat = df2['water'][1:][::-1].values.astype(float)
water2 = df2['Unnamed: 3'][1:][::-1].values.astype(float)
glucoseSat = df2['Unnamed: 5'][1:][::-1].values.astype(float)

#4.16M fresh cuvette w/ baseline
datafresh = pd.read_csv('/Users/davidl./Biophotonics/Data/SaturatedNewCuvette.csv')
dffresh = datafresh.drop(['Unnamed: 4'], axis=1)
wavelengthsFresh = dffresh['water'][1:][::-1].values.astype(float)
waterfresh = dffresh['Unnamed: 1'][1:][::-1].values.astype(float)
glucose416 = dffresh['Unnamed: 3'][1:][::-1].values.astype(float)

# Read the data from the first CSV file
data3 = pd.read_csv('/Users/davidl./Biophotonics/Data/hamExtracted.csv')
# Sort the data by the 'x' (wavelength) column
data3.sort_values(by='x', inplace=True)
# Update the sorted wavelengths and absorbance in the same variables
hamWavelengths = data3['x'].values.astype(int)
hamAbsorbance = data3['y'].values.astype(float)

# Read the data from the second CSV file
data4 = pd.read_csv('/Users/davidl./Biophotonics/Data/opticaExtracted.csv')
data83 = pd.read_csv('/Users/davidl./Biophotonics/Data/Fuglerud83mm.csv')
data500 = pd.read_csv('/Users/davidl./Biophotonics/Data/Fuglerud500mm.csv')
# Sort the data by the 'x' (wavelength) column
data4.sort_values(by='x', inplace=True)
data83.sort_values(by='x', inplace=True)
data500.sort_values(by='x', inplace=True)

# Update the sorted wavelengths and absorbance in the same variables
opticaWavelengths = data4['x'].values.astype(int)
opticaAbsorbance = data4['y'].values.astype(float)
optica83wavelengths = data83['x'].values.astype(int)
optica83 = data83['y'].values.astype(float)
optica500wavelengths = data500['x'].values.astype(int)
optica500 = data500['y'].values.astype(float)

# Load the +45 polarized input data
datapol = pd.read_csv('/Users/davidl./Biophotonics/Data/VH_polarized.csv')
datapol.drop(['Unnamed: 10'], axis=1)
wavelengthspol = datapol['water'][1:][::-1].values.astype(float)
water_V = datapol['Unnamed: 1'][1:][::-1].values.astype(float)
glucose_V = datapol['Unnamed: 5'][1:][::-1].values.astype(float)
water_H = datapol['Unnamed: 7'][1:][::-1].values.astype(float)
glucose_H = datapol['Unnamed: 9'][1:][::-1].values.astype(float)

# Load the unpolarized input data
datanopol = pd.read_csv('/Users/davidl./Biophotonics/Data/VH_nopolarization.csv')
datanopol.drop(['Unnamed: 8'], axis=1)
water_V_nopol = datanopol['Unnamed: 1'][1:][::-1].values.astype(float)
glucose_V_nopol = datanopol['Unnamed: 3'][1:][::-1].values.astype(float)
water_H_nopol = datanopol['Unnamed: 5'][1:][::-1].values.astype(float)
glucose_H_nopol = datanopol['Unnamed: 7'][1:][::-1].values.astype(float)

#HSA confouding factor
dataHSA = pd.read_csv('/Users/davidl./Biophotonics/Data/aqueousHSA.csv')
dataHSA.sort_values(by='x', inplace=True)
HSAwavelengths = dataHSA['x'].values.astype(int)
HSA = dataHSA['y'].values.astype(float) 

#Globulin confounding factor
dataGlob = pd.read_csv('/Users/davidl./Biophotonics/Data/globulin.csv')
dataGlob.sort_values(by='x', inplace=True)
globwavelengths = dataGlob['x'].values.astype(int)
glob = dataGlob['y'].values.astype(float)

#Ascorbate confouding factor
dataC = pd.read_csv('/Users/davidl./Biophotonics/Data/VitaminC_FTNIR.csv')
dataC.sort_values(by='x', inplace=True)
Cwavelengths = dataC['x'].values.astype(int)
ascorbate = dataC['y'].values.astype(float)

#Cholesterol Confounding factor
dataChol = pd.read_csv('/Users/davidl./Biophotonics/Data/CrystalCholesterol.csv')
dataChol.sort_values(by='x', inplace=True)
Cholwavelengths = dataChol['x'].values.astype(int)
cholesterol = dataChol['y'].values.astype(float)


def isolate(water_r, glucose_r, absorbance_g, absorbance_w, path_length, molarity):
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
    epsilon_w, epsilon_g = extinction(WATER_RATIO_2, GLUCOSE_RATIO_2, absorbance_w, absorbance_g, 0.02, 2)
    ag = (water_r*epsilon_w + glucose_r*epsilon_g)*path_length*molarity 
    water_component = (water_r)*epsilon_w*path_length*molarity
    return np.abs(ag - water_component)


def extinction(ratio_w, ratio_g, absorbance_w, absorbance_g, path_length, concentration):
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
    epsilon_w = absorbance_w / (path_length)
    epsilon_g = (absorbance_g/(path_length*concentration) - ratio_w*epsilon_w) / ratio_g
    return epsilon_w, epsilon_g

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

def plot(wavelengths, separated, water, glucose):
    """
    Plots spectra.
    """
    fig, ax = plt.subplots()
    separated = norm(separated)
    glucose = norm(glucose)
    water = norm(water)

    # Cut out wavelengths larger than 2400
    mask = (wavelengths >= 1100) & (wavelengths <= 2400)
    wavelengths_cut = wavelengths[mask]
    separated_cut = separated[mask]
    water_cut = water[mask]
    glucose_cut = glucose[mask]
    
    #separated, glucose, water = norm(separated, glucose, water)
    ax.plot(wavelengths_cut, separated_cut, label='Glucose Spectra')
    ax.plot(wavelengths_cut, glucose_cut, label='Solution Spectra')
    ax.plot(wavelengths_cut, water_cut, label='DI Spectra')

    ax.set_xlabel('Wavelength (nm)')
    ax.set_ylabel('Absorbance')
    ax.set_title('4.16M at 200 um path length')
    ax.legend()
    #ax.set_yscale('log')
    plt.show()

def plot_nosep(wavelengths, water, glucose):
    """
    Plots spectra.
    """
    fig, ax = plt.subplots()
    glucose = norm(glucose)
    water = norm(water)

    # Cut out wavelengths larger than 2400
    mask = (wavelengths >= 1100) & (wavelengths <= 2400)
    wavelengths_cut = wavelengths[mask]
    water_cut = water[mask]
    glucose_cut = glucose[mask]

    #separated, glucose, water = norm(separated, glucose, water)
    ax.plot(wavelengths_cut, glucose_cut, label='Solution Spectra')
    ax.plot(wavelengths_cut, water_cut, label='DI Spectra')

    ax.set_xlabel('Wavelength (nm)')
    ax.set_ylabel('Normalized Absorbance')
    ax.set_title('179 mg/dL at 2mm path length')
    ax.legend()
    #ax.set_yscale('log')
    plt.show()


def plot_extra(wavelengths, separated, water, glucose, pulled_wave, pulled, optica_wave, optica, optica83_wave, optica83, optica500_wave, optica500):
    """
    Plots spectra.
    """
    fig, ax = plt.subplots()
    pulled = norm(pulled)
    separated = norm(separated)
    water = norm(water)
    optica = norm(optica)
    optica83 = norm(optica83)
    optica500 = norm(optica500)

    # Cut out wavelengths larger than 2400
    mask = (wavelengths >= 1100) & (wavelengths <= 2400)
    wavelengths_cut = wavelengths[mask]
    separated_cut = separated[mask]
    water_cut = water[mask]
    glucose_cut = glucose[mask]

    ax.plot(wavelengths_cut, separated_cut, label='Glucose Spectra')
    #ax.plot(wavelengths_cut, glucose_cut, label='Solution Spectra')
    ax.plot(wavelengths_cut, water_cut, label='DI Spectra')
    ax.plot(list(pulled_wave), pulled, label='Ham et al') # Convert pulled_wave to a list
    #ax.plot(list(optica_wave), optica, label='Fuglerud et al (Commerical)')
    #ax.plot(list(optica83_wave), optica83, label='Fuglerud et al (83 mM)')
    ax.plot(list(optica500_wave), optica500, label='Fuglerud et al (500 mM)')      

    ax.set_xlabel('Wavelength (nm)')
    ax.set_ylabel('Normalized Absorbance')
    ax.set_title('Spectra Comparison')
    ax.legend()
    #ax.set_yscale('log')
    plt.show()

def plot_confounding(wavelengths, wavelengths_vitC, wavelengths_HSA, wavelengths_glob, wavelengths_chol, water, glucose, ascorbate, albumin, globulin, cholesterol):
    """Plots glucose along with confouding factors."""

    mask = (wavelengths >=1100) & (wavelengths <= 2500)
    wavelengths_cut = wavelengths[mask]
    water_cut = water[mask]
    glucose_cut = glucose[mask]
    #ascorbate_cut = ascorbate[mask]
    #albumin_cut = albumin[mask]
    #globulin_cut = globulin[mask]

    fig, ax = plt.subplots()
    water = norm(water_cut)
    glucose = norm(glucose_cut)
    ascorbate = norm(ascorbate)
    albumin = norm(albumin)
    globulin = norm(globulin)
    cholesterol = norm(cholesterol)
    #water, glucose, ascorbate, albumin, globulin = norm(water, glucose, ascorbate, albumin, globulin)

    ax.plot(wavelengths_cut, water, label='Water')
    ax.plot(wavelengths_cut, glucose, label="Glucose")
    ax.plot(wavelengths_vitC, ascorbate, label='Ascorbate')
    ax.plot(wavelengths_HSA, albumin, label='Albumin')
    ax.plot(wavelengths_glob, globulin, label='Globulin')
    ax.plot(wavelengths_chol, cholesterol, label='Cholesterol')

    ax.set_xlabel('Wavelength (nm)')
    ax.set_ylabel('Normalized Absorbance')
    ax.set_title('Glucose & Confounding Factors')
    ax.legend()
    plt.show()

def plot_both(wavelengths1, wavelengths2, glucose1, glucose2, water):
    
    mask = (wavelengths1 >=1100) & (wavelengths1 <= 2500)
    wavelengths1 = wavelengths1[mask]
    glucose1 = glucose1[mask]

    fig, ax = plt.subplots()
    glucose1 = norm(glucose1)
    glucose2 = norm(glucose2)
    water = norm(water)

    ax.plot(wavelengths1, glucose1, label='Old Cuvette')
    ax.plot(wavelengths2, glucose2, label='New Cuvette')
    ax.plot(wavelengths2, water, label='Water')

    ax.set_xlabel('Wavelength (nm)')
    ax.set_ylabel('Normalized Absorbance')
    ax.set_title('4.16M Glucose Solutions')
    ax.legend()
    plt.show()

def plot_polarization(wavelengths1, wavelengths2, glucose1, glucose2, water1, water2):

    fig, ax = plt.subplots()
    glucose1 = norm(glucose1)
    glucose2 = norm(glucose2)
    water1 = norm(water1)
    water2 = norm(water2)

    ax.plot(wavelengths1, glucose1, label='Vertically Polarized Glucose')
    #ax.plot(wavelengths2, glucose2, label='Horizontally Polarized Glucose')
    ax.plot(wavelengths1, water1, label='Vertically Polarized Water')
    #ax.plot(wavelengths1, water2, label='Horizontally Polarized Water')

    ax.set_xlabel('Wavelength (nm)')
    ax.set_ylabel('Normalized Absorbance')
    ax.set_title('4.16M Glucose Solutions')
    ax.legend()
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
    """Transforms % transmittance to absorbance."""
    spectra_t = np.clip(spectra, a_min=1e-10, a_max=None)
    #return 1 - spectra
    return -np.log(spectra_t)

def transform_opt(spectra):
    return 1 - spectra

def toORD(V1, V2):
    """Returns observed rotation (2 alpha) in milidegrees."""
    return np.abs((V1-V2)/(V1+V2))/2

def plot_ORD(wavelengths, milidegrees):
    """Plots ORD (milidegrees) vs wavelength."""
    fig, ax = plt.subplots()

    ax.plot(wavelengths, milidegrees, label='ORD')

    ax.set_xlabel('Wavelength (nm)')
    ax.set_ylabel('Observed Rotation')
    ax.set_title('4.16M Glucose Solutions')
    ax.legend()
    plt.show() 

if __name__=="__main__":
    print(datanopol)
    #separated_spectra_2M = isolate(WATER_RATIO_2, GLUCOSE_RATIO_2, glucose2M, water2M, 0.02, 2)
    separated_V = isolate(WATER_SAT, GLUCOSE_SAT, glucose_V, water_V, 0.02, 4.16)
    separated_H = isolate(WATER_SAT, GLUCOSE_SAT, glucose_H, water_H, 0.02, 4.16)
    ORD = toORD(separated_V, separated_H)
    ORD_water = toORD(water_V, water_H)
    separated_V_nopol = isolate(WATER_SAT, GLUCOSE_SAT, glucose_V_nopol, water_V_nopol, 0.02, 4.16)
    separated_H_nopol = isolate(WATER_SAT, GLUCOSE_SAT, glucose_H_nopol, water_H_nopol, 0.02, 4.16)
    #seperated_spectra_179 = isolate(WATER_RATIO_179, GLUCOSE_RATIO_179, solution179, waterf, 0.02, 0.00993561)
    #separated_spectra_2500 = isolate(WATER_RATIO_2500, GLUCOSE_RATIO_2500, glucose2500, water2500, 0.02, 0.0555)
    #separated_spectra_4M = isolate(WATER_SAT, GLUCOSE_SAT, glucoseSat, water2, 0.02, 4.16)
    #separated_spectra_4M_2 = isolate(WATER_SAT, GLUCOSE_SAT, glucose416, waterfresh, 0.02, 4.16)
    #plot(wavelengthsSat, transform(separated_spectra_4M), transform(water2), transform(glucoseSat))
    #plot(wavelengths, separated_spectra_2M, water2M, glucose2M)
    #plot_extra(wavelengthsFresh, transform(separated_spectra_4M_2), transform(waterfresh), glucose416, hamWavelengths, transform_opt(hamAbsorbance), opticaWavelengths, transform_opt(opticaAbsorbance), optica83wavelengths, optica83, optica500wavelengths, transform_opt(optica500))
    #plot_both(wavelengthsSat, wavelengthsFresh, transform(separated_spectra_4M), transform(separated_spectra_4M_2), transform(waterfresh))
    #plot_confounding(wavelengthsFresh, Cwavelengths, HSAwavelengths, globwavelengths, Cholwavelengths, transform(waterfresh), transform(separated_spectra_4M_2), ascorbate, HSA, glob, cholesterol)
    #plot_polarization(wavelengthspol, wavelengthspol, transform(separated_V), transform(separated_H), transform(water_V), transform(water_H))
    #lot_polarization(wavelengthspol, wavelengthspol, transform(separated_V_nopol), transform(separated_H_nopol), transform(water_V_nopol), transform(water_H_nopol))
    plot_ORD(wavelengthspol, ORD)
    plot_ORD(wavelengthspol, ORD_water)
