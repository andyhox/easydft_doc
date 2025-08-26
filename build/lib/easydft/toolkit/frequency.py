import re
import numpy as np
from ase.thermochemistry import HarmonicThermo, IdealGasThermo
from ase.build.molecule import molecule
from pymatgen.io.vasp.outputs import Outcar
from pprint import pprint


class AnalyticFreqProcessor:
    """
    VASP OUTCAR vibrational frequency post-processing tool.
    Supports automatic parsing of vibrational frequencies, imaginary frequencies, electronic energy, and calculation of thermochemical corrections for adsorbed and gas-phase molecules.
    """
    """
    Post-process vibrational frequency data from VASP OUTCAR files to compute thermochemical corrections,
    including symmetry number detection via PySCF.

    Attributes:
        outcar_path (str): Path to the OUTCAR file.
        temp (float): Temperature in Kelvin.
        pressure (float): Pressure in Pa.
        _epot (float): Electronic energy without entropy contributions.
        _vib_energies (np.ndarray): Array of real vibrational energies (eV).
        _imag_freqs (np.ndarray): Array of imaginary frequencies (cm^-1).
    """

    def __init__(self, outcar_path: str, temp: float = 298.15, pressure: float = 101325.0):
        """
        Initialize the frequency analysis processor.
        
        Args:
            outcar_path (str): Path to the OUTCAR file.
            temp (float): Temperature (K), default 298.15.
            pressure (float): Pressure (Pa), default 101325.0.
        """
        self.outcar_path = outcar_path
        self.temp = temp
        self.pressure = pressure
        self._epot = None
        self._vib_energies = None
        self._imag_freqs = None

    @property
    def epot(self) -> float:
        """
        Read the electronic energy (without entropy) from the OUTCAR file.
        
        Returns:
            float: Electronic energy (eV).
        """
        if self._epot is None:
            outcar = Outcar(self.outcar_path)
            self._epot = outcar.final_energy_wo_entrp
        return self._epot

    def _parse_frequencies(self) -> None:
        """
        Parse real and imaginary frequencies from the OUTCAR file.
        
        Sets:
            self._vib_energies (np.ndarray): Real vibrational energies (eV)
            self._imag_freqs (np.ndarray): Imaginary frequencies (cm^-1)
        """
        real_energies = []
        imag_freqs = []
        pattern_imag = re.compile(r"([-+]?[0-9]*\.?[0-9]+)\s+cm-1")

        with open(self.outcar_path, 'r') as f:
            for line in f:
                if 'f/i=' in line and 'THz' in line:
                    match = pattern_imag.search(line)
                    if match:
                        imag_freqs.append(float(match.group(1)))
                elif 'f  =' in line and 'cm-1' in line:
                    parts = line.split()
                    try:
                        energy_meV = float(parts[-2])
                        real_energies.append(energy_meV / 1000.0)
                    except (ValueError, IndexError):
                        continue

        self._vib_energies = np.array(real_energies)
        self._imag_freqs = np.array(imag_freqs)

    @property
    def vib_energies(self) -> np.ndarray:
        """
        Get real vibrational energies (eV), automatically parse if not already done.
        
        Returns:
            np.ndarray: Array of real vibrational energies (eV).
        """
        if self._vib_energies is None:
            self._parse_frequencies()
        return self._vib_energies

    @property
    def imag_freqs(self) -> np.ndarray:
        """
        Get imaginary frequencies (cm^-1), automatically parse if not already done.
        
        Returns:
            np.ndarray: Array of imaginary frequencies (cm^-1).
        """
        if self._imag_freqs is None:
            self._parse_frequencies()
        return self._imag_freqs

    def compute_adsorbate_thermo(self) -> dict:
        """
        Calculate the zero-point energy and Helmholtz free energy correction for adsorbates.
        
        Returns:
            dict: {'zpe': float, 'helmholtz_energy': float}
        """
        thermo = HarmonicThermo(
            vib_energies=self.vib_energies,
            potentialenergy=self.epot,
            ignore_imag_modes=True
        )
        return {
            'zpe': thermo.get_ZPE_correction(),
            'helmholtz_energy': thermo.get_helmholtz_energy(self.temp)
        }

    def compute_gas_thermo(
        self,
        spin: float,
        atoms,
        symmetry_number: int,
        geometry: str = 'linear'
    ) -> dict:
        """
        Calculate thermodynamic functions for ideal gas molecules.
        
        Args:
            spin (float): Number of unpaired electrons / 2.
            atoms: ASE Atoms object or list of atomic symbols.
            symmetry_number (int): Molecular symmetry number (can be obtained by parse_symmetrynumber).
            geometry (str): 'linear'/'nonlinear'/'monoatomic'.
        Returns:
            dict: {'zpe', 'gibbs_energy', 'enthalpy', 'entropy'}
        """
        thermo = IdealGasThermo(
            vib_energies=self.vib_energies,
            potentialenergy=self.epot,
            ignore_imag_modes=True,
            spin=spin,
            symmetrynumber=symmetry_number,
            geometry=geometry,
            atoms=atoms
        )
        return {
            'zpe': thermo.get_ZPE_correction(),
            'gibbs_energy': thermo.get_gibbs_energy(self.temp, self.pressure),
            'enthalpy': thermo.get_enthalpy(self.temp),
            'entropy': thermo.get_entropy(self.temp, self.pressure)
        }


# if __name__ == '__main__':
#     # 示例运行
#     h2_atoms = molecule('H2')
#     processor = AnalyticFreqProcessor('/public/home/39041fcf59/ehukai/AA/COFs/H2/zpe/OUTCAR')
#     gas_results = processor.compute_gas_thermo(spin=0, atoms=h2_atoms, symmetry_number=2, geometry='linear')
#     pprint(gas_results)
