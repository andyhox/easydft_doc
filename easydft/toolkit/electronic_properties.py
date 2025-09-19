from pymatgen.io.vasp.outputs import Vasprun, Chgcar
from pymatgen.electronic_structure.core import OrbitalType, Spin, Orbital
from pymatgen.electronic_structure.bandstructure import BandStructureSymmLine, BandStructure
from pymatgen.core.periodic_table import Element
from easydft.core.parser import VasprunParser
from pathlib import Path
import pandas as pd
import numpy as np
from typing import Optional, Any
import os
        

class Analyzer():
    def __init__(
        self, 
        work_dir: None,
        file_path: None, 
        save_data: bool = True, 
        output_path: Optional[str] = None,
        ):
        self.work_dir = work_dir
        self.file_path = file_path
        self.output_path = output_path
        self.save_data = save_data
        self._validate_input()
    def _validate_input(self):
        if not self.work_dir and not self.file_path:
            raise ValueError("Either 'work_dir' or 'input_file' must be provided.")
        if self.save_data:
            if self.output_path is None:
                raise ValueError("output_path MUST be provided")
            if not os.path.isdir(self.output_path):
                raise NotADirectoryError(f"{self.output_path} is not a directory")

    def save_to_csv(
        self,
        df: pd.DataFrame, 
        filename: str
        ) -> None:
        
        if self.save_data:
            if self.output_path is None:
                raise ValueError("output_path MUST be provided")
            save_path = os.path.join(str(self.output_path), filename)
            df.to_csv(save_path, index=False)
        
class DosAnalyzer(Analyzer):
    """
    VASP density of states (DOS) parser.
    Supports extraction of element/orbital projected DOS and output as pandas DataFrame.
    """
    
    ORBITAL_MAP = {
        's': OrbitalType.s,
        'p': OrbitalType.p,
        'd': OrbitalType.d,
        'f': OrbitalType.f,
    }

    SUB_ORBITAL_MAP = {
        "s": Orbital.s,
        "p_y": Orbital.py,
        "p_z": Orbital.pz,
        "p_x": Orbital.px,
        "d_{xy}": Orbital.dxy,
        "d_{yz}": Orbital.dyz,
        "d_{z^2}": Orbital.dz2,
        "d_{xz}": Orbital.dxz,
        "d_{x^2-y^2}": Orbital.dx2,
        "f_{y(3x^2-y^2)}": Orbital.f_3,
        "f_{xyz}": Orbital.f_2,
        "f_{yz^2}": Orbital.f_1,
        "f_{z^3}": Orbital.f0,
        "f_{xz^2}": Orbital.f1,
        "f_{z(x^2-y^2)}": Orbital.f2,
        "f_{x(x^2-3y^2)}": Orbital.f3,
    }

    def __init__(
        self, 
        work_dir: Optional[str] = None,
        vasprun_file: Optional[str] = None, 
        save_data: bool = True, 
        output_path: Optional[str] = None
        ):
        
        """
        Initialize the DOS parser.
        
        Args:
            vasprun_file (str): Path to the VASP vasprun.xml file.
        """
        super().__init__(work_dir=work_dir, file_path=vasprun_file, save_data=save_data, output_path=output_path)
        self._parse_vasprun_file()
        self.vasprun_parser = VasprunParser(self.vasprun_file)
        self.vasprun = self.vasprun_parser.vasprun
        self.dos = self.vasprun.complete_dos
        self.symbols = list(set(self.vasprun.atomic_symbols))
        self.is_spin = self.vasprun.is_spin
        self.lorbit = self.vasprun.parameters['LORBIT']

    def _parse_vasprun_file(self):
        if self.work_dir:
            vasprun_file = os.path.join(self.work_dir, 'vasprun.xml')
            if os.path.exists(vasprun_file):
                self.vasprun_file = vasprun_file
            else:
                raise FileNotFoundError(f"vasprun.xml not found in work_dir: {self.work_dir}")
        elif self.file_path:
            if os.path.exists(self.file_path):
                self.vasprun_file = self.file_path
            else:
                raise FileNotFoundError(f"vasprun_file not found: {self.file_path}")
        else:
            raise ValueError("Either work_dir or vasprun_file must be provided")
        
        return self.vasprun_file
        
    def _add_dos_data(self, data, key, densities):
        """
        Add DOS data to the data dictionary (supports spin polarization).
        
        Args:
            data (dict): Data dictionary
            key (str): Data key name
            densities (dict): Spin-resolved densities
        Returns:
            dict: Updated data dictionary
        """
        if self.is_spin:
            data[f"{key}_up"] = densities[Spin.up]
            data[f"{key}_down"] = densities[Spin.down]
        else:
            data[key] = densities[Spin.up]
        return data

    def _is_valid_dos(self, densities_dict):
        """
        Determine whether the DOS data is valid (not all zeros).
        
        Args:
            densities_dict (dict): Spin-resolved densities
        Returns:
            bool: Whether valid
        """
        if densities_dict is None:
            return False
        for spin in [Spin.up, Spin.down] if self.is_spin else [Spin.up]:
            if (spin in densities_dict and 
                densities_dict[spin] is not None and 
                not np.allclose(densities_dict[spin], 0)):
                return True
        return False

    def parse_element_spd_dos(self) -> pd.DataFrame:
        """
        Parse element-projected (s/p/d) DOS.
        
        Returns:
            pd.DataFrame: DataFrame containing energy, total DOS, and element-projected DOS.
        """
        data = {'Energy': self.dos.energies - self.dos.efermi}
        data = self._add_dos_data(data, 'tdos', self.dos.densities)

        for element in self.symbols:
            element_spd_dos = self.dos.get_element_spd_dos(Element(element))
            
            for orb_name, orb_type in self.ORBITAL_MAP.items():
                if orb_type in element_spd_dos:
                    densities = element_spd_dos[orb_type].densities
                    if self._is_valid_dos(densities):
                        data = self._add_dos_data(
                            data, 
                            f"{element}_{orb_name}", 
                            densities
                        )
        df = pd.DataFrame(data)
        self.save_to_csv(df, "element_spd_dos.csv")
        return df
    
    def parse_orbital_spd_dos(self) -> pd.DataFrame:
        """
        Parse orbital-projected DOS by summing over all atomic sites of each element.

        Returns:
            pd.DataFrame: A DataFrame with energy and summed orbital DOS for each element.
        """
        
        if self.lorbit < 11:
            raise ValueError(f"Orbital-projected DOS requires LORBIT >= 11.")

        data = {'Energy': self.dos.energies - self.dos.efermi}
        data = self._add_dos_data(data, 'tdos', self.dos.densities)
        
        # Initialize an empty dict to accumulate summed densities
        orbital_dos_sum = {
            element: {orb_label: {Spin.up: 0.0, Spin.down: 0.0} for orb_label in self.SUB_ORBITAL_MAP}
            for element in self.symbols
        }

        # Loop over all sites
        for site in self.vasprun.final_structure.sites:
            element = site.specie.symbol

            for orb_label, orb in self.SUB_ORBITAL_MAP.items():
                # Prefer newer API that requires explicit orbital; fallback to older mapping-based API
                try:
                    dos_obj = self.dos.get_site_orbital_dos(site=site, orbital=orb)
                except TypeError:
                    get_site_orbital_dos_fn = getattr(self.dos, "get_site_orbital_dos")
                    site_dos_map = get_site_orbital_dos_fn(site)
                    dos_obj = site_dos_map.get(orb, None)

                if dos_obj is None:
                    continue

                densities_dict = getattr(dos_obj, "densities", None)
                if densities_dict is None:
                    continue

                for spin in [Spin.up, Spin.down] if self.is_spin else [Spin.up]:
                    densities = densities_dict.get(spin)
                    if densities is None:
                        continue
                    if isinstance(orbital_dos_sum[element][orb_label][spin], float):
                        orbital_dos_sum[element][orb_label][spin] = np.zeros_like(densities)
                    orbital_dos_sum[element][orb_label][spin] += densities

        # Add each element-orbital's DOS into data
        for element, orb_dos in orbital_dos_sum.items():
            for orb_label, spin_dict in orb_dos.items():
                if self._is_valid_dos(spin_dict):
                    data = self._add_dos_data(
                        data,
                        f"{element}_{orb_label}",
                        spin_dict
                    )

        df = pd.DataFrame(data)
        self.save_to_csv(df, "orbital_spd_dos.csv")
        return df

    def parse_element_dos(self) -> pd.DataFrame:
        """
        Parse element total DOS.
        
        Returns:
            pd.DataFrame: DataFrame containing energy, total DOS, and element total DOS.
        """
        data = {'Energy': self.dos.energies - self.dos.efermi}
        
        data = self._add_dos_data(data, 'tdos', self.dos.densities)

        element_dos_dict = self.dos.get_element_dos()
        for element in self.symbols:
            el = Element(element)
            if el in element_dos_dict:
                densities = element_dos_dict[el].densities
                if self._is_valid_dos(densities):
                    data = self._add_dos_data(data, str(el), densities)
                    
        df = pd.DataFrame(data)
        self.save_to_csv(df, "element_dos.csv")
        return df

class BandStructureAnalyzer(Analyzer):
    def __init__(
        self,
        work_dir: Optional[str] = None, 
        vasprun_file: Optional[str] = None, 
        save_data: bool=True,
        output_path: Optional[str] = None
        ):
        
        super().__init__(work_dir=work_dir, file_path=vasprun_file, save_data=save_data, output_path=output_path)

        self.vasprun_parser = VasprunParser(vasprun_file)
        self.vasprun = self.vasprun_parser.vasprun
        self._bandstructure = self.vasprun.get_band_structure(line_mode=True)
        self._kpoint_distances = self._calculate_kpoint_distances()

    @property
    def bandstructure(self) -> BandStructureSymmLine | BandStructure:
        return self._bandstructure
    @property
    def is_spin(self):
        return self.vasprun.is_spin
    
    @property
    def kpoints(self):
        return self.bandstructure.kpoints
    
    @property
    def efermi(self):
        return self.bandstructure.efermi
    
    @property
    def bands(self):
        return self.bandstructure.bands
    
    @property
    def band_gap(self):
        return self.bandstructure.get_band_gap()['energy']
    
    @property
    def is_gap_direct(self):
        return self.bandstructure.get_band_gap()['direct']
    
    @property
    def vbm(self):
        return self.bandstructure.get_vbm()['energy']
    
    @property
    def cbm(self):
        return self.bandstructure.get_cbm()['energy']
    @property
    def nkpoints(self) -> int:
        return len(self.kpoints)
    def _calculate_kpoint_distances(self):
        distances = [0.0]
        for i in range(1, self.nkpoints):
            dist = float(np.linalg.norm(
                self.kpoints[i].frac_coords - self.kpoints[i-1].frac_coords
                )
            )
            distances.append(distances[-1] + dist)
        return distances
    
    def parse_normal_bandstructure(self) -> pd.DataFrame:
        """
        Parse normal bandstructure.
        """
        base_data = []
        for band_idx in range(len(self.bands)):
            for k_idx in range(self.nkpoints):
                kpoint = self.kpoints[k_idx]
                kpoint_info = kpoint.label if kpoint.label else ",".join(f"{c:.4f}" for c in kpoint.frac_coords)
                
                base_data.append({
                    "band_index": band_idx,
                    "k_idx": k_idx,
                    "kpoint": kpoint_info,
                    "kpoint_distance": self._kpoint_distances[k_idx]
                })
        
        base_df = pd.DataFrame(base_data)
        
        if not self.is_spin:
            energy_data = []
            for band_idx in range(len(self.bands)):
                for k_idx in range(self.nkpoints):
                    energy = self.bands[Spin.up][band_idx, k_idx] - self.efermi
                    energy_data.append({
                        "band_index": band_idx,
                        "k_idx": k_idx,
                        "energy": energy
                    })
            energy_df = pd.DataFrame(energy_data)
            combined_df = pd.merge(base_df, energy_df, on=["band_index", "k_idx"])
            
        else:
            up_data = []
            for band_idx in range(len(self.bands)):
                for k_idx in range(self.nkpoints):
                    energy_up = self.bands[Spin.up][band_idx, k_idx] - self.efermi
                    up_data.append({
                        "band_index": band_idx,
                        "k_idx": k_idx,
                        "energy_up": energy_up
                    })
            up_df = pd.DataFrame(up_data)
            
            down_data = []
            for band_idx in range(len(self.bands)):
                for k_idx in range(self.nkpoints):
                    energy_down = self.bands[Spin.down][band_idx, k_idx] - self.efermi
                    down_data.append({
                        "band_index": band_idx,
                        "k_idx": k_idx,
                        "energy_down": energy_down
                    })
            down_df = pd.DataFrame(down_data)
            
            combined_df = pd.merge(base_df, up_df, on=["band_index", "k_idx"])
            combined_df = pd.merge(combined_df, down_df, on=["band_index", "k_idx"])
        
        self.save_to_csv(combined_df, "band.csv")
        return combined_df
    
# TODO: 
#   - 实现投影能带（projected bandstructure）功能
#   - 支持元素/轨道分波能带的提取与输出

class ChargeAnalyzer:
    
    @staticmethod
    def _parse_dim(path: Path):
        vasprun_file = path / "vasprun.xml"
        vasprun = Vasprun(vasprun_file)
        return (
            vasprun.parameters["NGX"],
            vasprun.parameters["NGY"],
            vasprun.parameters["NGZ"],
        )
    
    @staticmethod
    def _check_shape(shapes: list[tuple[int, int, int]]) -> bool:
        if not shapes:
            return False
        reference = shapes[0]
        return all(shape == reference for shape in shapes)
    
    @classmethod
    def DifferenceChargeDensity(
        cls,
        root_path: Optional[str] = None,
        save_dir: Optional[str] = None,
    ):
        root = Path(root_path).resolve()
        save = Path(save_dir or root).resolve()
        save.mkdir(parents=True, exist_ok=True)

        dirs = {"AB": root / "AB", "A": root / "A", "B": root / "B"}

        shapes = [cls._parse_dim(d) for d in dirs.values()]
        if not cls._check_shape(shapes):
            raise ValueError(f"Inconsistent FFT grid dimensions: {shapes}")

        chg = {k: Chgcar.from_file(d / "CHGCAR") for k, d in dirs.items()}

        diff = (
            chg["AB"]
            .linear_add(chg["A"], scale_factor=-1)
            .linear_add(chg["B"], scale_factor=-1)
        )

        out_file = save / "chgdiff.cube"
        diff.to_cube(str(out_file))
        return out_file
    
    @classmethod
    def DeformationChargeDensity(
        cls,
        root_path: Optional[str] = None,
        save_dir: Optional[str] = None,
    ):
        root = Path(root_path).resolve()
        save = Path(save_dir or root).resolve()
        save.mkdir(parents=True, exist_ok=True)

        dirs = {"scf": root / "scf", "atom": root / "atom"}

        shapes = [cls._parse_dim(d) for d in dirs.values()]
        if not cls._check_shape(shapes):
            raise ValueError(f"Inconsistent FFT grid dimensions: {shapes}")

        chg = {k: Chgcar.from_file(d / "CHGCAR") for k, d in dirs.items()}

        diff = (
            chg["scf"]
            .linear_add(chg["atom"], scale_factor=-1)
        )

        out_file = save / "chgdiff.cube"
        diff.to_cube(str(out_file))
        return out_file