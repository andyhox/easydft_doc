from pymatgen.core import Structure
from pymatgen.io.vasp.outputs import Oszicar
from pymatgen.io.vasp.inputs import Potcar
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
        idpp_params: dict = None) -> List[Structure]:
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
        solver = IDPPSolver.from_endpoints(
            endpoints=[initial_structure, final_structure],
            nimages=nimages,
            sort_tol=1.0
        )

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

        optimized_path = solver.run(**default_params)
        
        for i, structure in enumerate(optimized_path):
            image_dir = os.path.join(output_path, f"{i:02d}")
            os.makedirs(image_dir, exist_ok=True)
            
            poscar_path = os.path.join(image_dir, "POSCAR")
            structure.to(filename=str(poscar_path), fmt="poscar")
        
        
        return optimized_path

class NEB_analysis:
    """
    NEB path energy barrier and diffusion property analysis tool.
    """
    
    @classmethod
    def get_pomass_from_potcar(cls, potcar_path: str, element: str):
        """
        Get the atomic mass of the specified element from the POTCAR file.
        
        Args:
            potcar_path (str): Path to the POTCAR file
            element (str): Element symbol
        Returns:
            float: Atomic mass (amu)
        """
        with open(potcar_path, 'r') as f:
            content = f.read()
        
        pattern = re.compile(
            r'VRHFIN\s*=\s*(\w+):.*?POMASS\s*=\s*([\d.]+)',
            re.DOTALL
        )
        mass_dict = {m.group(1): float(m.group(2)) 
                for m in pattern.finditer(content)}
        return mass_dict[element]
    @classmethod
    def analysis(
        cls,
        path: Path,
        nimages: int,
        temp: float = None,
        element: str = None,
        calc_D_factor: bool = False,
        ) -> pd.DataFrame:
        """
        NEB path energy barrier and diffusion property analysis.
        
        Args:
            path (Path): Directory containing NEB images
            nimages (int): Number of intermediate images
            temp (float): Temperature (K)
            element (str): Element for diffusion coefficient calculation
            calc_D_factor (bool): Whether to calculate the diffusion factor
        Returns:
            pd.DataFrame: Table of energy barrier and diffusion properties
        Raises:
            FileNotFoundError: If the path or OSZICAR file does not exist
            ValueError: If required parameters are missing
        """
        path = Path(path) if isinstance(path, str) else path
        if not path.exists():
            raise FileNotFoundError(f"Path {path} does not exist")
            
        if calc_D_factor and (temp is None or element is None):
            raise ValueError("temp, element and mass required for calc_D_factor")
        
        kb = 1.380649e-23
        
        sort_data = {}
        energy_list = []
        
        images_dirs = [path/f"{i:02d}" for i in range(nimages+2)]
        for i,folder in enumerate(images_dirs):
            oszicar_path = folder/"OSZICAR"
            if not oszicar_path.exists():
                raise FileNotFoundError(f"OSZICAR not found in {folder}")
            oszicar = Oszicar(oszicar_path)
            final_energy = oszicar.final_energy
            energy_list.append(final_energy)
            sort_data[f"E_{i}"] = final_energy
            
        Ea = max(energy_list) - min(energy_list)
        sort_data["Energy_barrier"] = Ea
        
        if calc_D_factor:
            initial_structure = Structure.from_file(path/"00"/"CONTCAR")
            final_structure = Structure.from_file(path/f"{nimages+1:02d}"/"CONTCAR")
            
            disp = np.array([s2.distance(s1) 
                            for s1, s2 in zip(initial_structure, final_structure)])
            mst_disp = np.sqrt(np.sum(disp**2))
            l = mst_disp*1e-8
            
            Ea_J = Ea * 1.60218e-19
            
            potcar = os.path.join(path, "00", "POTCAR")
            mass = cls.get_pomass_from_potcar(potcar, element)*1.660539e-27
                    
            v = ((2*Ea_J)/(mass*l**2))**(1/2)
            exp = math.exp(-Ea_J/(kb*temp))
            D0 = v*l**2 
            D = D0*exp 
            factor = get_conversion_factor(initial_structure, element, temp)
            sigam = D * factor
            sort_data.update({
                "hop_displacement": mst_disp,
                "D0": D0,
                "D": D,
                "factor": factor,
                "sigma": sigam
                })

            df = pd.DataFrame([sort_data])
            return df
        else:
            df = pd.DataFrame([sort_data])
            return df