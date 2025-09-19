from pymatgen.analysis.surface_analysis import WorkFunctionAnalyzer
from pymatgen.io.vasp.outputs import Outcar
from pymatgen.io.vasp.inputs import Poscar
from typing import Optional
import pandas as pd
import matplotlib.pyplot as plt
import os


class WorkFuncAnalyzer:
    """
    A class for analyzing and visualizing the work function from VASP calculation outputs.

    This class extracts the planar-averaged electrostatic potential along the z-direction
    from VASP output files (POSCAR, LOCPOT, OUTCAR), computes the work function, and
    provides options to save the data and plot the results.

    """

    def __init__(
        self,
        calc_dir: Optional[str] = None,
        poscar_file: Optional[str] = None,
        locpot_file: Optional[str] = None,
        outcar_file: Optional[str] = None,
        output_path: Optional[str] = None,
        save_data: bool = True,
        plot_fig: bool = False,
    ):
        """
        Initialize the WorkFunction class.

        Args:
            calc_dir (str or None): The directory containing the calculation files (POSCAR, LOCPOT, OUTCAR).
            poscar_file (str or None): Path to the POSCAR file. Used if calc_dir is not provided.
            locpot_file (str or None): Path to the LOCPOT file. Used if calc_dir is not provided.
            outcar_file (str or None): Path to the OUTCAR file. Used if calc_dir is not provided.
            output_path (str): Path to save the output data or figures.
            save_data (bool, optional): Whether to save the extracted data. Default is True.
            plot_fig (bool, optional): Whether to plot the work function figure. Default is False.

        Raises:
            FileExistsError: If the required files are not found in the specified directory or at the given paths.
        """
        plt.switch_backend('agg')

        self.calc_dir = calc_dir
        self.poscar_file = poscar_file
        self.locpot_file = locpot_file
        self.outcar_file = outcar_file
        self.output_path = output_path
        self.save_data = save_data
        self.plot_fig = plot_fig
        self.df = None
        self._validate_input()
        self._parse_input()
        self._extract_data()

        if self.output_path:
            os.makedirs(self.output_path, exist_ok=True)
            if self.save_data and self.df is not None:
                self.save_to_csv()
            if self.plot_fig and self.df is not None:
                self.plot_wf()

    def _validate_input(self):
        """
        Validate the input parameters to ensure that either a calculation directory or all required file paths are provided.

        Raises:
            ValueError: If neither calc_dir nor all of poscar_file, locpot_file, and outcar_file are provided.
        """
        if not self.calc_dir and not (self.poscar_file and self.locpot_file and self.outcar_file):
            raise ValueError("Either 'calc_dir' or all of 'poscar_file', 'locpot_file', 'outcar_file' must be provided.")

    def _parse_input(self):
        """
        Parse and validate input file paths.

        If calc_dir is provided, construct file paths for POSCAR, LOCPOT, and OUTCAR from this directory.
        Otherwise, use the provided file paths directly.

        Raises:
            FileExistsError: If any required file does not exist.
        """
        if self.calc_dir:
            self.poscar_file = os.path.join(self.calc_dir, "POSCAR")
            self.locpot_file = os.path.join(self.calc_dir, "LOCPOT")
            self.outcar_file = os.path.join(self.calc_dir, "OUTCAR")

        required_files = [self.poscar_file, self.locpot_file, self.outcar_file]
        for file in required_files:
            if not file or not os.path.exists(file):
                raise FileExistsError(f"Required file not found: {file}")

    def _extract_data(self):
        """
        Extract the planar-averaged electrostatic potential along the z-direction.

        Returns:
            tuple: A tuple containing:
                - pd.DataFrame: DataFrame with columns 'Z_direction' and 'Electrostatic_potential'.
                - WorkFunctionAnalyzer: The analyzer object containing work function analysis results.

        Raises:
            Exception: If extraction fails.
        """
        try:
            wf_analyzer = WorkFunctionAnalyzer.from_files(self.poscar_file, self.locpot_file, self.outcar_file)
            distance = wf_analyzer.along_c
            es_potential = wf_analyzer.locpot_along_c
            lattice_c = Poscar.from_file(self.poscar_file).structure.lattice.c
            df = pd.DataFrame({
                'Z_direction': [d * lattice_c for d in distance], 
                'Electrostatic_potential': es_potential
                })
        except Exception as e:
            raise e
        self.df = df
        self.wf_analyzer = wf_analyzer
        return self.df, self.wf_analyzer

    @property
    def work_function(self) -> float:
        """
        Get the calculated work function value.

        Returns:
            float: The work function in eV.
        """
        return self.wf_analyzer.work_function
    
    @property
    def vac_level(self) -> float:
        """
        Get the vacuum electrostatic potential level.

        Returns:
            float: The vacuum level in eV.
        """
        return self.wf_analyzer.vacuum_locpot
    
    @property
    def efermi(self) -> float:
        """
        Get the Fermi energy level.

        Returns:
            float: The Fermi level in eV.
        """
        return self.wf_analyzer.efermi

    def save_to_csv(self):
        """
        Save the extracted electrostatic potential data to a CSV file.

        The file is saved as 'workfunction.csv' in the output_path directory.

        Raises:
            Exception: If saving fails.
        """
        if not self.output_path or self.df is None:
            return
        csv_path = os.path.join(self.output_path, 'workfunction.csv')
        try:
            self.df.to_csv(csv_path, index=False)
        except:
            raise Exception(f'fail to save file to {self.output_path}')

    def plot_wf(self):
        """
        Plot the planar-averaged electrostatic potential and work function.

        The plot includes the Fermi level, vacuum level, and the calculated work function.
        The figure is saved as 'planar_average_plot.png' in the output_path directory.

        Raises:
            Exception: If plotting fails.
        """
        if not self.plot_fig or self.df is None:
            return
        try:
            work_func = self.work_function
            vac_level = self.vac_level
            efermi = self.efermi
            plt.figure(figsize=(9, 6))
            plt.plot(self.df['Z_direction'], self.df['Electrostatic_potential'], 'b-', linewidth=2)
            plt.axhline(y=efermi, color='r', linestyle='--', linewidth=1.5,
                        label=f'E-fermi = {efermi:.4f} eV')
            plt.axhline(y=vac_level, color='g', linestyle='--', linewidth=1.5,
                        label=f'Vacuum level = {vac_level:.4f} eV')
            plt.annotate(f'WorkFunction = {work_func:.4f} eV',
                        xy=(self.df['Z_direction'].max() * 0.5, (vac_level + efermi) * 0.5),
                        xytext=(10, 10),
                        textcoords='offset points',
                        fontsize=12,
                        ha='center')
            plt.xlabel("Z direction (Å)", fontsize=20)
            plt.xticks(fontsize=18)
            plt.ylabel("Electrostatic potential (eV)", fontsize=20)
            plt.yticks(fontsize=18)
            # plt.grid(True, linestyle='--', alpha=0.7)
            plt.legend(fontsize=12)
            plt.tight_layout()

            output_file = os.path.join(self.output_path, 'planar_average_plot.png')
            plt.savefig(output_file, dpi=300, bbox_inches='tight')
            plt.close()
        except:
            raise Exception(f'fail to plot figure')