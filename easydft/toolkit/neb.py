from pymatgen.core import Structure
from pymatgen.io.vasp.outputs import Oszicar
from pymatgen.io.vasp.inputs import Potcar
from pymatgen.analysis.transition_state import NEBAnalysis
from pymatgen.analysis.diffusion.neb.pathfinder import IDPPSolver
from pymatgen.analysis.diffusion.analyzer import get_conversion_factor
from pathlib import Path
from typing import List
import pandas as pd
import numpy as np
import math
import os
import re

class IDPP_path:
    """
    NEB path generation tool based on the IDPP algorithm.
    """
        
    @classmethod
    def gen_path(
        cls, 
        output_path: Path, 
        initial_structure: Structure, 
        final_structure: Structure, 
        nimages: int = 5,
        idpp_params: dict = None,
        save_output: bool = True) -> List[Structure]:
        """
        Generate the IDPP path and save as POSCAR files.
        
        Args:
            output_path (Path): Directory to save NEB path
            initial_structure (Structure): Initial structure
            final_structure (Structure): Final structure
            nimages (int): Number of intermediate images (total images = nimages+2)
            idpp_params (dict): IDPP algorithm parameters
        Returns:
            List[Structure]: List of optimized path structures
        """
        # Generate the IDPP path
        solver = IDPPSolver.from_endpoints(
            endpoints=[initial_structure, final_structure],
            nimages=nimages,
            sort_tol=1.0
        )

        # 设置默认IDPP参数
        default_params = {
            "maxiter": 1000,
            "tol": 1e-5,
            "gtol": 0.001,
            "step_size": 0.05,
            "max_disp": 0.05,
            "spring_const": 5.0,
        }
        if idpp_params:
            default_params.update(idpp_params)

        # Get optimized path
        optimized_path = solver.run(**default_params)
        
        if save_output:
            # Save each image
            for i, structure in enumerate(optimized_path):
                image_dir = os.path.join(output_path, f"{i:02d}")
                os.makedirs(image_dir, exist_ok=True)
                
                poscar_path = os.path.join(image_dir, "POSCAR")
                structure.to(filename=str(poscar_path), fmt="poscar")
        
        
        return optimized_path

class NEB_analysis:
    """
    Class for analyzing NEB (Nudged Elastic Band) paths, including energy barrier calculation and diffusion property analysis.

    This class provides methods to:
    - Load and analyze NEB calculation results from a directory.
    - Extract and plot the energy profile along the NEB path.
    - Calculate the energy barrier.
    - Optionally, compute diffusion-related properties such as hop displacement, pre-exponential factor (D0), diffusion coefficient (D), conversion factor, and ionic conductivity.

    Attributes:
        path (str or Path): The directory containing NEB calculation results.
        results (NEBAnalysis): The NEBAnalysis object containing parsed NEB results.

    Methods:
        __init__(path): Initialize the NEB_analysis object and load NEB results.
        _neb_results(path): Load NEB results from the specified directory.
        plot: Property to get the matplotlib figure of the NEB energy profile.
        get_pomass_from_potcar(potcar_path, element): Class method to extract atomic mass from a POTCAR file.
        sort_results(...): Analyze the NEB path to obtain the energy barrier and, optionally, diffusion properties.
    """

    def __init__(
        self,
        path,
        ):
        """
        Initialize the NEB_analysis object.

        Args:
            path (str or Path): Path to the directory containing NEB calculation results.
        """
        self.path = path
        self._neb_results()

    def _neb_results(self):
        """
        Load NEB results from the specified directory.

        Args:
            path (str or Path): Path to the NEB calculation directory.

        Returns:
            NEBAnalysis: The NEBAnalysis object containing parsed NEB results.
        """
        results = NEBAnalysis.from_dir(self.path)
        self.results = results
        return self.results

    def plot(self, title: str = None):
        """
        Plot the NEB (Nudged Elastic Band) energy profile.

        This method generates and returns a matplotlib Axes object showing the energy profile along the NEB path.
        Optionally, a title can be set for the plot.

        Args:
            title (str, optional): The title to display on the plot. Defaults to None.

        Returns:
            matplotlib.axes.Axes: The matplotlib Axes object containing the NEB energy profile plot.
        """
        ax = self.results.get_plot()
        ax.set_title(title, fontsize=26)
        return ax

    @classmethod
    def get_atomic_mass_from_potcar(cls, potcar_path: str, element: str):
        """
        Get the atomic mass of a specified element from a POTCAR file.

        This method parses the POTCAR file to extract the atomic mass (POMASS) for the given element symbol.
        The mass is returned in kilograms.

        Args:
            potcar_path (str): The file path to the POTCAR file.
            element (str): The element symbol (e.g., "Li", "Na", "O").

        Returns:
            float: The atomic mass of the specified element in kilograms.

        Raises:
            KeyError: If the specified element is not found in the POTCAR file.
            FileNotFoundError: If the POTCAR file does not exist.
        """
        with open(potcar_path, 'r') as f:
            content = f.read()

        pattern = re.compile(
            r'VRHFIN\s*=\s*(\w+):.*?POMASS\s*=\s*([\d.]+)',
            re.DOTALL
        )
        mass_dict = {m.group(1): float(m.group(2))
                for m in pattern.finditer(content)}
        return mass_dict[element]*1.660539e-27

    def sort_results(
        self,
        calc_D_factor: bool = False,
        nimages: int = None,
        temp: float = None,
        element: str = None,
        element_mass: float = None,
        ) -> pd.DataFrame:
        """
        Analyze the NEB path to obtain the energy barrier and, optionally, diffusion properties.

        Args:
            calc_D_factor (bool, optional): Whether to calculate diffusion properties (hop displacement, D0, D, factor, and conductivity).
            nimages (int, optional): Number of intermediate images (not including initial and final).
            temp (float, optional): Temperature in Kelvin, required if calc_D_factor is True.
            element (str, optional): Element symbol for diffusion coefficient calculation, required if calc_D_factor is True.
            element_mass (float, optional): Atomic mass of the element (in kg), required if calc_D_factor is True.

        Returns:
            pd.DataFrame: DataFrame containing the energy barrier and, if requested, diffusion properties.

        Raises:
            FileNotFoundError: If the specified path or required files do not exist.
            ValueError: If required parameters for diffusion calculation are missing.

        Notes:
            - The energy barrier is calculated as the difference between the maximum and minimum energies along the NEB path.
            - If calc_D_factor is True, the function also computes the hop displacement, pre-exponential factor (D0), diffusion coefficient (D), conversion factor, and ionic conductivity (sigma).
        """

        if calc_D_factor and (temp is None or element is None or element_mass is None or nimages is None):
            raise ValueError("nimages, temp, element and mass required for calc_D_factor")

        # Boltzmann constant
        kb = 1.380649e-23

        # dict for storing all results
        sort_data = {}

        energy_list = self.results.energies

        # diffusion energy barrier
        Ea = max(energy_list) - min(energy_list)
        sort_data["Energy_barrier(eV)"] = Ea

        if calc_D_factor:
            # get hop displacements
            initial_structure = Structure.from_file(os.path.join(self.path,"00","POSCAR"))
            final_structure = Structure.from_file(os.path.join(self.path,f"{nimages+1:02d}","POSCAR"))

            disp = np.array([s2.distance(s1)
                            for s1, s2 in zip(initial_structure, final_structure)])
            # unit: Angstrom
            mst_disp = np.sqrt(np.sum(disp**2))
            # unit: cm
            l = mst_disp*1e-8

            # eV to J
            Ea_J = Ea * 1.60218e-19

            # vibration frequency.
            v = ((2*Ea_J)/(element_mass*l**2))**(1/2)
            exp = math.exp(-Ea_J/(kb*temp))
            # diffusion constant. Unit: (cm^2/s)
            D0 = v*l**2
            # diffusion coefficient. Unit: (cm^2/s)
            D = D0*exp
            factor = get_conversion_factor(initial_structure, element, temp)
            # Conductivity. Unit: (mS/cm)
            sigam = D * factor
            sort_data.update({
                "hop_displacement(angstrom)": mst_disp,
                "D0(cm^2/s)": D0,
                "D(cm^2/s)": D,
                "factor": factor,
                "sigma(mS/cm)": sigam
                })

            df = pd.DataFrame([sort_data])
            return df
        else:
            df = pd.DataFrame([sort_data])
            return df