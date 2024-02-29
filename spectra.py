"""Class implementation of spectral processing functions"""

"""

*** Design principles ***

An instance of the spectra class should contain all data collected during a single
data collection event, independent of the experimental procedure. The load_spectrum 
method should only called once after class object instatiantion. 

A core priniciple is the elimination of unnecessary data movement. Any method implemented
within the class should not require the user to pass copies of the data as parameters.

"""


"""NOTE: I MADE A SIGNIFICANT CHANGE TO THE CODE. NO NORMALIZATION FOR THE TIME BEING"""




import numpy as np
import pandas as pd
from typing import List, Tuple, Union
import matplotlib.pyplot as plt
class Spectra:
    """
    Contains the spectra for processing and display. An instance of 
    this class is designed to contain spectra collected in one experiment.

    Examples:
    # Create a Spectra instance with solution parameters
    Saturated_spectra = Spectra(water_ratio=0.9875, 
                                solution_ratio=0.0125, 
                                path_length=0.02, 
                                concentration=0.056)

    # Load preprocessed spectra and wavelengths into the class instance
    wavelengths = wavelengthsSat  
    Saturated_spectra.load_spectrum('water', wavelengths, water)
    Saturated_spectra.load_spectrum('solution solution', wavelengts, solution_solution)

    # Get the names of loaded spectra
    spectra_names = Saturated_spectra.get_spectra_names()
    print("Loaded spectra names:", spectra_names)

    # Access a specific spectrum
    spectrum_name = 'water'
    spectrum_info = Saturated_spectra.get_spectrum(spectrum_name)
    print(f"Spectrum '{spectrum_name}':\n", spectrum_info)

    # Calculate and set molar absorptivity coefficients
    Saturated_spectra.absorptivity(water2, solutionSat)

    # Demix the spectral components of solution solution
    separated_solute = Saturated_spectra.demix(water2, solutionSat)

    # Plot the demixed solution spectrum
    Spectra.plot_spectra(wavelengths, [separated_solute], ['Demixed solution'], 'Demixed solution Spectrum')
    """

    def __init__(self, water_ratio: float, solution_ratio: float, path_length: float, concentration: float):
        """
        Init variables:
            water_ratio: mass ratio of water in solution solution
            solution_ratio: mass ratio of solute in the solution.
            path_length: the pathlength of the container (cm)
            concentration: the molarity of solution in the solution (mol/liter)
        """

        """Spectral data needs to be populated with """

        self.spectra_data = {}
        self.water_ratio = water_ratio
        self.solution_ratio = solution_ratio
        self.path_length = path_length
        self.concentration = concentration

    def load_spectrum(self, name: str, wavelengths: np.ndarray, *spectra_data: np.ndarray):
        spectra_info = [{'wavelength': w, 'spectrum_data': s}
                        for w, s in zip(wavelengths, spectra_data)]
        self.spectra_data[name] = spectra_info

    def get_spectra_names(self):
        return list(self.spectra_data.keys())

    def get_spectrum(self, name: str) -> Union[None, List[dict]]:
        return self.spectra_data.get(name, None)

    def absorptivity(self, water: np.ndarray, solution: np.ndarray):
        """
        Creates a set of molar absorptivity coefficients for water and the solute.

        Input variables:
            water: Spectral data of water
            solution: Spectral data of solution.
            normalize: True if you want to normalize the coefficients. Useful for plotting.

        Returns:
            Molar absorptivity coefficients for water and solution.
        """
        if self.path_length == 0:
            raise ValueError(
                "Path length cannot be zero for absorptivity calculation.")
        if self.concentration == 0:
            raise ValueError(
                "Concentration cannot be zero for absorptivity calculation.")
        if np.isnan(self.path_length) or np.isnan(self.concentration) or np.isnan(self.water_ratio) or np.isnan(self.solution_ratio):
            raise ValueError(
                "NaN values detected in absorptivity calculation parameters.")

        epsilon_w = water / self.path_length
        epsilon_solution = (solution / (self.path_length * self.concentration) -
                            self.water_ratio * epsilon_w) / self.solution_ratio

        if np.isnan(epsilon_w).any() or np.isnan(epsilon_solution).any():
            raise ValueError(
                "NaN values detected in absorptivity calculation results.")

        self.epsilon_water = epsilon_w
        self.epsilon_solution = epsilon_solution

    def demix(self, water: np.ndarray, solution: np.ndarray) -> np.ndarray:
        """
        Function demixes the spectral components of solution solution. 

        Input variables:
            water: Spectral data of water.
            solution: Spectral data of solution solution.

        Returns:
            Array of separated solution spectra. Note that this is not automatically added to self.spectra_data. 
        """
        if not hasattr(self, 'epsilon_water') or not hasattr(self, 'epsilon_solution'):
            raise ValueError(
                "Epsilon values for water and solution are not set. Please execute absorptivity() first.")

        if water.shape != solution.shape:
            raise ValueError(
                "Water and solution spectra must have the same shape.")

        ag = (self.water_ratio * self.epsilon_water + self.solution_ratio *
              self.epsilon_solution) * self.path_length * self.concentration
        water_component = self.water_ratio * self.epsilon_water * \
            self.path_length * self.concentration

        if ag.shape != water_component.shape:
            raise ValueError(
                "Ag and water component arrays must have the same shape.")

        try:
            if np.abs(ag - water_component) is not None:
                self.demixed_spectra = np.abs(ag - water_component)
        except Exception as e:
            print(e)

    def _normalize(*args: np.ndarray) -> Union[np.ndarray, Tuple[np.ndarray, ...]]:
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
            try:
                vectors = args
                concatenated = np.concatenate(vectors)
                max_val = np.max(concatenated)
                min_val = np.min(concatenated)
                normalized_vectors = [
                    (vector - min_val) / (max_val - min_val) for vector in vectors]
            except Exception as e:
                print(e)
            else:
                print("Concatenation successful")
                return tuple(normalized_vectors)

    def filter_data(self):
        """
        Function removes NaN columns from spectra_data.
        """
        for name, spectra_info in self.spectra_data.items():
            for spectrum in spectra_info:
                spectrum_data = spectrum['spectrum_data']
                nan_indices = np.isnan(spectrum_data)
                non_nan_indices = ~nan_indices.any(axis=0)
                spectrum['spectrum_data'] = spectrum_data[:, non_nan_indices]

    def plot_spectra(self, spectrum_names: List[str], wavelengths: List[float], labels: List[str], title: str, transform: bool):
        """
        Plots spectra.

        Input variables:
            spectrum_names: List of names of the spectra to plot.
            wavelengths: List of wavelengths (x-axis). 
            labels: List of labels for each spectrum.
            title: Title for the plot.
        """
        fig, ax = plt.subplots()
        spectrum_data = None

        for spectrum_name, label in zip(spectrum_names, labels):
            spectrum_info = self.get_spectrum(spectrum_name)
            if self.demixed_spectra is not None:
                spectrum_data = self.demixed_spectra
            """RIGHT HERE IS WHERE NORMALIZATION WOULD GO"""
            if spectrum_data is not None:
                # spectrum_data = np.array(spectrum_data)
                if transform:
                    ax.plot(wavelengths, self.transform(
                        spectrum_data), label=label)
                else:
                    ax.plot(wavelengths, spectrum_data, label=label)
        ax.set_xlabel('Wavelength (nm)')
        ax.set_ylabel('Absorbance')
        ax.set_title(title)
        ax.legend()
        plt.show()

    def _plot_spectra(self):
        """
        A improved attempt at plotting spectra. No data should ever be supplied to this
        function. All data is contained within an instance of the class. 
        """
        fig, ax = plt.subplots()

        """
        TODO: You first need to fix the demix function. Do not create a new dataType there
        """

    def transform(self, spectra: np.ndarray) -> np.ndarray:
        """Transfroms betweeen transmittance and absorbance."""
        spectra = np.clip(spectra, a_min=1e-10, a_max=None)
        return 1 - spectra
