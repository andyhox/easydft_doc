from pymatgen.core import Element
from pymatgen.core.structure import Structure, Molecule
from pymatgen.core.lattice import Lattice
from pymatgen.core.surface import SlabGenerator, Slab
from pymatgen.analysis.interfaces.coherent_interfaces import CoherentInterfaceBuilder
from pymatgen.analysis.interfaces.zsl import ZSLGenerator
from pymatgen.transformations.advanced_transformations import SlabTransformation
from pymatgen.io.ase import AseAtomsAdaptor
from pymatgen.core.operations import SymmOp
from pymatgen.analysis.adsorption import AdsorbateSiteFinder, plot_slab
from ase.build.molecule import molecule
from typing import Union, Any, Optional, List
from itertools import combinations
import numpy as np
import random
import math
import os

class BulkModify:
    """
    Crystal structure strain processing tool.
    Supports applying various strains to the structure, including lattice, volume, and angle strains.
    """
    
    @classmethod
    def add_lattice_strain(cls, structure: Slab|Structure, strain_list: list, mode: str) -> dict:
        """
        Apply strain in the specified direction to the initial structure.
        
        Args:
            structure (Structure): Initial structure
            strain_list (list): List of strains, e.g. [-0.01, 0, 0.01]
            mode (str): Strain mode, such as 'a', 'b', 'c', or 'all'
        Returns:
            dict: Key is the strain value, value is the strained structure
        """
        lattice = structure.lattice
        strained_structures = {}

        if mode == 'a':
            strain_indices = [0]
        elif mode == 'b':
            strain_indices = [1]
        elif mode == 'c':
            strain_indices = [2]
        elif mode == 'all':
            strain_indices = [0, 1, 2]
        else:
            raise ValueError("Invalid mode. Use 'a', 'b', 'c', or 'all'.")

        for strain in strain_list:
            # 初始晶格参数
            new_lattice_params = [lattice.a, lattice.b, lattice.c]
            
            # 根据mode选择施加应变的轴
            for i in strain_indices:
                new_lattice_params[i] = new_lattice_params[i] * (1 + strain)

            new_lattice = Lattice.from_parameters(*new_lattice_params, lattice.alpha, lattice.beta, lattice.gamma)
            new_structure = Structure(new_lattice, structure.species, structure.frac_coords)
            strained_structures[strain] = new_structure

        return strained_structures
    
    @classmethod
    def add_volume_strain(cls, structure: Structure, strain_list: list) -> dict:
        """
        Apply isotropic volume strain to the initial structure.
        
        Args:
            structure (Structure): Initial structure
            strain_list (list): List of strains, e.g. [-0.01, 0, 0.01]
        Returns:
            dict: Key is the strain value, value is the strained structure
        """
        lattice = structure.lattice
        strained_structures = {}

        for strain in strain_list:
            # 计算体积应变因子，并取三次方根
            strain_factor = (1 + strain) ** (1/3)
            
            # 新的晶格参数，等比例缩放
            new_lattice_params = [lattice.a * strain_factor, 
                                  lattice.b * strain_factor, 
                                  lattice.c * strain_factor]

            new_lattice = Lattice.from_parameters(*new_lattice_params, lattice.alpha, lattice.beta, lattice.gamma)
            new_structure = Structure(new_lattice, structure.species, structure.frac_coords)
            strained_structures[strain] = new_structure

        return strained_structures
    
    @classmethod
    def add_angle_strain(cls, structure: Structure, strain_list: list, mode: str) -> dict:
        """
        Apply angle strain to the initial structure.
        
        Args:
            structure (Structure): Initial structure
            strain_list (List[float]): List of strains, e.g. [-0.01, 0, 0.01]
            mode (str): Strain mode, can be 'alpha', 'beta', or 'gamma'
        Returns:
            dict: Key is the strain value, value is the strained structure
        """
        lattice = structure.lattice
        strained_structures = {}

        # 选择施加应变的角度
        if mode == 'alpha':
            angle_index = 0
        elif mode == 'beta':
            angle_index = 1
        elif mode == 'gamma':
            angle_index = 2
        else:
            raise ValueError("Invalid mode. Use 'alpha', 'beta', 'gamma'.")

        # 获取初始的角度参数
        initial_angles = [lattice.alpha, lattice.beta, lattice.gamma]

        for strain in strain_list:
            # 复制初始角度
            new_angles = initial_angles[:]
            
            # 对指定角度施加应变
            new_angles[angle_index] *= (1 + strain)

            # 创建新的晶格
            new_lattice = Lattice.from_parameters(lattice.a, lattice.b, lattice.c, *new_angles)
            new_structure = Structure(new_lattice, structure.species, structure.frac_coords)
            strained_structures[strain] = new_structure

        return strained_structures

class SlabModify:
    """
    Slab structure generation and modification tool.
    Supports slab generation, splitting, fixing, strain, and other operations.
    """

    @staticmethod
    def fix_orthogonal(slab: Slab):
        """
        Orthogonalize the slab lattice.
        
        Args:
            slab (Slab): Slab structure
        Returns:
            Slab: Orthogonalized slab
        """
        return slab.get_orthogonal_c_slab()

    @classmethod
    def gen_all_slabs(cls,
        structure: Structure,
        miller_index: tuple,
        min_slab_size: Union[int, float],
        min_vacuum_size: Union[int, float],
    ):
        """
        Generate all symmetrically equivalent slab structures.
        
        Args:
            structure (Structure): Initial structure
            miller_index (tuple): Miller index
            min_slab_size (float): Slab thickness
            min_vacuum_size (float): Vacuum thickness
        Returns:
            list: List of slab structures
        """
        slabgen = SlabGenerator(initial_structure=structure, 
                                miller_index=miller_index, 
                                min_slab_size=min_slab_size, 
                                min_vacuum_size=min_vacuum_size,
                                center_slab=True,
                                in_unit_planes=True,)
        slabs = slabgen.get_slabs(symmetrize=True)
        return slabs
    
    @classmethod
    def amorphous_gen_slab(cls,
                           structure, 
                           miller_index: tuple,
                           min_slab_size: Union[int, float],
                           min_vacuum_size: Union[int, float], 
                           shift: float=0):
        """
        Generate an amorphous slab structure.
        
        Args:
            structure (Structure): Initial structure
            miller_index (tuple): Miller index
            min_slab_size (float): Slab thickness
            min_vacuum_size (float): Vacuum thickness
            shift (float): Slab shift
        Returns:
            Slab: Slab structure
        """
        slabgen = SlabTransformation(miller_index=miller_index,
                                     min_slab_size=min_slab_size,
                                     min_vacuum_size=min_vacuum_size,
                                     center_slab=True,
                                     in_unit_planes=True,
                                     shift=shift,
                                     )
        slab = slabgen.apply_transformation(structure=structure)
        if slab.lattice.alpha != 90 and slab.lattice.beta != 90:
            Warning("Slab is not orthogonal, fix it.")
            slab = SlabModify.fix_orthogonal(slab)
        return slab

    @classmethod
    def gen_slab(cls,structure: Structure, 
                 miller_index: tuple, 
                 min_slab_size: Union[int, float], 
                 min_vacuum_size: Union[int, float],
                 shift: float=0.0,
                 in_unit_planes: bool = True,
                 ):
        """
        Generate a slab structure.
        
        Args:
            structure (Structure): Initial structure
            miller_index (tuple): Miller index
            min_slab_size (float): Slab thickness
            min_vacuum_size (float): Vacuum thickness
            shift (float): Slab shift
            in_unit_planes (bool): Whether to generate in unit planes
        Returns:
            Slab: Slab structure
        """
        slabgen = SlabGenerator(initial_structure=structure,
                                miller_index=miller_index,
                                min_slab_size=min_slab_size,
                                min_vacuum_size=min_vacuum_size,
                                in_unit_planes=in_unit_planes,)

        slab = slabgen.get_slab(shift=shift)
        if slab.lattice.alpha != 90 and slab.lattice.beta != 90:
            Warning("Slab is not orthogonal, fix it.")
            slab = SlabModify.fix_orthogonal(slab)
        return slab

    @classmethod
    def fix_slab(cls, 
                 slab: Structure, 
                 min_fix_ratio: Union[float, 0.0] = 0.0,
                 max_fix_ratio: Union[float, 1.0] = 1.0,
                 ):
        """
        Fix part of the atoms in the slab structure (selective dynamics).
        
        Args:
            slab (Structure): Slab structure
            min_fix_ratio (float): Fractional coordinate threshold for fixing lower atoms
            max_fix_ratio (float): Fractional coordinate threshold for fixing upper atoms (optional)
        Returns:
            Structure: Structure with selective_dynamics property added
        """

        fixed_atoms = []
        if max_fix_ratio is None:
            # 标记固定的原子
            for site in slab.sites:
                # 如果原子的分数坐标在要固定的层数范围内，则固定原子
                if site.frac_coords[2] <= min_fix_ratio:
                    fixed_atoms.append([False, False, False])  # 固定
                else:
                    fixed_atoms.append([True, True, True])  # 不固定

            # 为 slab 添加选择性动态属性
            slab.add_site_property("selective_dynamics", fixed_atoms)
            return slab
        else:
            for site in slab.sites:
                # 如果原子的分数坐标在要固定的层数范围内，则固定原子
                if site.frac_coords[2] <= min_fix_ratio or site.frac_coords[2] >= max_fix_ratio:
                    fixed_atoms.append([False, False, False])  # 固定
                else:
                    fixed_atoms.append([True, True, True])  # 不固定

            # 为 slab 添加选择性动态属性
            slab.add_site_property("selective_dynamics", fixed_atoms)
            return slab

    @classmethod
    def add_strain(cls, structure: Slab|Structure, strain_list: list, mode: str):
        """
        Apply strain to the slab structure.
        
        Args:
            structure (Slab|Structure): Slab or structure
            strain_list (list): List of strains
            mode (str): Strain mode
        Returns:
            dict: Dictionary of strained structures
        """
        return BulkModify.add_lattice_strain(structure=structure, strain_list=strain_list, mode=mode)
    
    @classmethod
    def split_mol(cls,
                   slab: Slab|Structure,
                   boundary_frac: float,
                   cubic_size: int
                   ):
        """
        Split the slab by the interface to generate molecule and base structures.
        
        Args:
            slab (Slab|Structure): Slab structure
            boundary_frac (float): Interface fraction
            cubic_size (int): Molecule box size
        Returns:
            (mol_structure, base_structure): Molecule structure and base structure
        """
        # 分解结构
        mol_sites = [site for site in slab.sites if site.frac_coords[2] > boundary_frac]
        base_sites = [site for site in slab.sites if site.frac_coords[2] <= boundary_frac]
        
        # 生成mol盒子，base结构
        mol_structure = Structure.from_sites(mol_sites)
        new_lattice = Lattice.cubic(cubic_size)
        # 计算分子在笛卡尔坐标下的中心
        cart_coords = mol_structure.cart_coords
        center = np.mean(cart_coords, axis=0)
        # 计算平移向量，将分子中心移到盒子中心
        shift_vector = np.array([cubic_size/2, cubic_size/2, cubic_size/2]) - center
        # 对所有原子进行平移
        new_cart_coords = cart_coords + shift_vector
        # 将新坐标转换为在新晶格下的分数坐标
        new_frac_coords = new_lattice.get_fractional_coords(new_cart_coords)
        # 用新的晶格、原有的原子种类和新分数坐标构造新的结构
        mol_structure = Structure(new_lattice, mol_structure.species, new_frac_coords)
        
        base_structure = Structure.from_sites(base_sites)
        
        return mol_structure, base_structure
    
    @classmethod
    def split_heter(
        cls,
        slab: Slab|Structure,
        boundary_frac: float,
        ):
        """
        Split the heterojunction slab by the interface to generate film and substrate structures.
        
        Args:
            slab (Slab|Structure): Slab structure
            boundary_frac (float): Interface fraction
        Returns:
            (film_structure, substrate_structure): Film structure and substrate structure
        """
        # 分解结构
        film_sites = [site for site in slab.sites if site.frac_coords[2] > boundary_frac]
        substrate_sites = [site for site in slab.sites if site.frac_coords[2] <= boundary_frac]
        
        # 生成mol盒子，base结构
        film_structure = Structure.from_sites(film_sites)
        substrate_structure = Structure.from_sites(substrate_sites)
        
        return film_structure, substrate_structure
    
class HeterojunctionModify:
    """
    Heterojunction interface generation and modification tool.
    """
    def __init__(
        self, 
        film_structure: Structure,
        substrate_structure: Structure,
        film_miller_index: tuple,
        substrate_miller_index: tuple,
        zslgen = ZSLGenerator()):
        """
        Initialize the heterojunction modifier.
        
        Args:
            film_structure (Structure): Film structure
            substrate_structure (Structure): Substrate structure
            film_miller_index (tuple): Film Miller index
            substrate_miller_index (tuple): Substrate Miller index
            zslgen: ZSLGenerator instance
        """
        self.film_structure = film_structure
        self.substrate_structure = substrate_structure
        self.film_miller_index = film_miller_index
        self.substrate_miller_index = substrate_miller_index
        self.zslgen = zslgen
        self.cib = self.gen_cib()
        self.terminations = self.cib.terminations

    def gen_cib(self):
        """
        Generate a CoherentInterfaceBuilder object.
        
        Returns:
            CoherentInterfaceBuilder: Interface builder
        """
        return CoherentInterfaceBuilder(substrate_structure=self.substrate_structure, 
                                        film_structure=self.film_structure, 
                                        film_miller=self.film_miller_index, 
                                        substrate_miller=self.substrate_miller_index,
                                        zslgen=self.zslgen)

    def gen_heter_interface(self,
                            termination: tuple,
                            gap: float = 2.0,
                            vacuum_over_film: float = 20.0,
                            film_thickness: float = 1,
                            substrate_thickness: float = 1,
                            in_layers: bool = True):
        """
        Generate heterojunction interface structures.
        
        Args:
            termination (tuple): Termination
            gap (float): Gap
            vacuum_over_film (float): Vacuum thickness above the film
            film_thickness (float): Film thickness
            substrate_thickness (float): Substrate thickness
            in_layers (bool): Whether in units of layers
        Returns:
            list: List of interface structures
        Raises:
            ValueError: Raised if no interface is generated
        """
        
        termination = tuple(termination)

        interfaces = self.cib.get_interfaces(
            termination=termination, gap=gap, vacuum_over_film=vacuum_over_film, 
            film_thickness=film_thickness, substrate_thickness=substrate_thickness, 
            in_layers=in_layers)
        
        if not interfaces:
            raise ValueError("No interfaces generated by gen_heter_interface")

        return interfaces
    
class AbsorptionModify:
    """
    Molecular adsorption related structure processing tool.
    """
    @classmethod
    def ase2pmg(cls, absorption_molecule: str) -> Molecule:
        '''
        Convert an ASE molecule object to a pymatgen Molecule.
        
        Args:
            absorption_molecule (str): Molecule name
        Returns:
            Molecule: pymatgen Molecule object
        '''
        ase_molecule = molecule(absorption_molecule)
        pmg_molecule = AseAtomsAdaptor.get_molecule(ase_molecule)
        return pmg_molecule

    @classmethod
    def find_adsorption_sites(cls, slab: Slab|Structure, adsorp_gap: float = 2) -> dict:
        '''
        Find adsorption sites on the slab surface.
        
        Args:
            slab (Slab|Structure): Slab structure
            adsorp_gap (float): Adsorption height
        Returns:
            dict: Adsorption site information
        '''
        asf = AdsorbateSiteFinder(slab)
        adsorption_sites = asf.find_adsorption_sites(distance=adsorp_gap)
        return adsorption_sites

def gen_substitute_structures(
    structure: Union[Slab, Structure],
    substitute_element: Union[Element, str],
    dopant: Union[Element, str, None],
    dopant_num: int,
    max_structures: int = 20, 
    save: bool = True,
    save_path: str = None,
) -> List[Structure]:
    """
    Generate structures by substituting or removing specific elements in a given structure.

    Args:
        structure (Union[Slab, Structure]): The input structure (pymatgen Slab or Structure object).
        substitute_element (Union[Element, str]): The element to be substituted (can be Element or string).
        dopant (Union[Element, str, None]): The dopant element to substitute in, or None/"vac"/"vacancy" for vacancy.
        dopant_num (int): Number of sites to substitute or remove.
        max_structures (int, optional): Maximum number of generated structures to return. Default is 20.
        save (bool, optional): Whether to save the generated structures as CIF files. Default is True.
        save_path (str, optional): Directory to save the CIF files. Required if save is True.

    Returns:
        List[Structure]: A list of generated structures with the specified substitutions or vacancies.

    Notes:
        - Symmetry operations are ignored.
        - If dopant is None, "vac", "vacancy", or "none", the selected sites will be removed (vacancy).
        - If dopant_num is 0, returns a copy of the original structure.
        - If the number of possible combinations exceeds max_structures, a random subset is returned.
        - If save is True and save_path is provided, the generated structures are saved as CIF files in the specified directory.
    """
    if isinstance(substitute_element, str):
        substitute_element = Element(substitute_element)

    if isinstance(dopant, str):
        if dopant.lower() in ["vac", "vacancy", "none"]:
            dopant = None
        else:
            dopant = Element(dopant)

    substitute_indices = [
        i for i, site in enumerate(structure)
        if site.specie == substitute_element
    ]
    total_substitute = len(substitute_indices)

    if total_substitute < dopant_num:
        raise ValueError(f"number of {substitute_element} sites are less than {dopant_num}")

    if dopant_num == 0:
        return [structure.copy()]  

    structures = []
    seen_combos = set()

    if save and save_path:
        os.makedirs(save_path, exist_ok=True)

    max_possible = math.comb(total_substitute, dopant_num)
    n_to_generate = min(max_structures, max_possible)

    while len(structures) < n_to_generate:
        combo = tuple(sorted(random.sample(substitute_indices, dopant_num)))
        if combo in seen_combos:
            continue
        seen_combos.add(combo)

        new_struct = structure.copy()
        if dopant is None:
            for idx in sorted(combo, reverse=True):
                del new_struct[idx]
        else:
            for idx in combo:
                new_struct.replace(idx, dopant)

        if save and save_path:
            filename = f"struct{len(structures)+1}.cif"
            new_struct.to(filename=os.path.join(save_path, filename))

        structures.append(new_struct)

    return structures