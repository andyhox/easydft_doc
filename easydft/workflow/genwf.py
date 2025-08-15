"""
This module is used to generate VASP workflows, supporting various tasks such as structure optimization, static calculation, non-self-consistent calculation, adsorption energy, and adhesion energy.
"""
from pymatgen.core import Structure, Lattice
from pymatgen.core.surface import Slab
from atomate2.vasp.sets.core import RelaxSetGenerator, StaticSetGenerator, NonSCFSetGenerator
from atomate2.vasp.jobs.core import RelaxMaker, StaticMaker, NonSCFMaker
from atomate2.vasp.flows.core import DoubleRelaxMaker
from atomate2.vasp.run import JobType
from atomate2.vasp.powerups import add_metadata_to_flow
from jobflow import Flow, Maker, job
from jobflow.managers.fireworks import flow_to_workflow
from jobflow_remote import submit_flow
from fireworks import LaunchPad
from datetime import datetime
from easydft.toolkit.structure_modify import SlabModify
from typing import Optional
import numpy as np
import warnings

class GenWF:
    """
    VASP workflow generator

    Used to generate different types of VASP calculation workflows based on the input structure, including double relaxation, relaxation+static, SCF+NonSCF, adsorption energy, adhesion energy, etc.
    """
    
    def __init__(
        self, 
        structure: Structure | Slab,
        use_custodian: bool = True
        ):
        """
        Initialize the VASP workflow generator.

        Args:
            structure (Structure): Input crystal structure
            use_custodian (bool): Whether to use custodian for error handling
        """
        self.structure = structure
        self.use_custodian = use_custodian
        self.jobtype = JobType.NORMAL if use_custodian else JobType.DIRECT
    
    def _default_name(self) -> str:
        """
        Generate the default workflow name in the format "formula_date".

        Returns:
            str: Default name
        """
        formula = self.structure.composition.reduced_formula
        # 获取当前日期
        date_str = datetime.now().strftime("%Y%m%d")
        return f"{formula}_{date_str}"

    def _get_flow_name(self, suffix: str, flow_name: Optional[str] = None) -> str:
        """
        Get the workflow name. If not specified, generate automatically and issue a warning.

        Args:
            suffix (str): Name suffix
            flow_name (str): User-specified name
        Returns:
            str: Workflow name
        """
        if flow_name is None:
            flow_name = f"{self._default_name()}_{suffix}"
            warnings.warn(f"flow_name is not defined, use default name: {flow_name}")
        return flow_name
    
    def double_relax_flow(
        self, 
        flow_name: Optional[str] = None,
        set1: Optional[RelaxSetGenerator] = None,
        set2: Optional[RelaxSetGenerator] = None,
        ) -> Flow:
        """
        Double relaxation workflow.

        Args:
            flow_name (str): Workflow name
            job_name (str): Job name
            set1 (VaspInputSet): VaspInputSet for the first relaxation
            set2 (VaspInputSet): VaspInputSet for the second relaxation
        Returns:
            Flow: Double relaxation workflow
        """
        flow_name = self._get_flow_name("double_relax", flow_name)
        
        if set1 is None:
            set1 = RelaxSetGenerator()
            warnings.warn(f'set1 is not defined, use default MPRelaxSet()')
        if set2 is None:
            set2 = RelaxSetGenerator()
            warnings.warn(f'set2 is not defined, use default MPRelaxSet()')
        # 初始化Maker
        relax1_maker = RelaxMaker(input_set_generator=set1,run_vasp_kwargs={"job_type":self.jobtype})
        relax2_maker = RelaxMaker(input_set_generator=set2,run_vasp_kwargs={"job_type":self.jobtype})
        doublerelax_maker = DoubleRelaxMaker(relax_maker1=relax1_maker, relax_maker2=relax2_maker)
        
        # 创建任务链
        doublerelax_job = doublerelax_maker.make(self.structure)
        flow = Flow(doublerelax_job, name=flow_name)
        flow = add_metadata_to_flow(flow, {"flow_name": flow_name})
        return flow
    
    def relax2scf_flow(
        self,
        flow_name: Optional[str] = None,
        relaxset: Optional[RelaxSetGenerator] = None,
        scfset: Optional[StaticSetGenerator] = None,
        ) -> Flow:
        """
        Relaxation + static calculation workflow.

        Args:
            flow_name (str): Workflow name
            relaxset (VaspInputSet): Input set for relaxation
            scfset (VaspInputSet): Input set for static calculation
        Returns:
            Flow: Relaxation + static calculation workflow
        """
        flow_name = self._get_flow_name("relax2scf", flow_name)
        
        if relaxset is None:
            relaxset = RelaxSetGenerator()
            warnings.warn(f'relaxset is not defined, use default MPRelaxSet()')
        if scfset is None:
            scfset = StaticSetGenerator()
            warnings.warn(f'staticset is not defined, use default MPStaticSet()')
        # 初始化Maker
        relax_maker = RelaxMaker(input_set_generator=relaxset,run_vasp_kwargs={"job_type":self.jobtype})
        static_maker = StaticMaker(input_set_generator=scfset,run_vasp_kwargs={"job_type":self.jobtype})
        
        # 创建任务链
        relax_job = relax_maker.make(self.structure)
        static_job = static_maker.make(structure=relax_job.output.structure, prev_dir=relax_job.output.dir_name)
        jobs = [relax_job, static_job]
        flow = Flow(jobs, output=static_job.output, name=flow_name, )
        flow = add_metadata_to_flow(flow, {"flow_name": flow_name})
        
        return flow
    
    def scf2nonscf_flow(
        self, 
        flow_name: Optional[str] = None,
        scfset: Optional[StaticSetGenerator] = None,
        nonscfset: Optional[NonSCFSetGenerator] = None,
        ) -> Flow:
        """
        SCF + NonSCF calculation workflow.

        Args:
            flow_name (str): Workflow name
            scfset (VaspInputSet): Input set for SCF calculation
            nonscfset (VaspInputSet): Input set for NonSCF calculation
        Returns:
            Flow: SCF + NonSCF workflow
        """
        flow_name = self._get_flow_name("scf2nonscf", flow_name)
        
        if scfset is None:
            scfset = StaticSetGenerator()
            warnings.warn(f'relaxset is not defined, use default MPRelaxSet()')
        if nonscfset is None:
            nonscfset = NonSCFSetGenerator()
            warnings.warn(f'staticset is not defined, use default MPStaticSet()')
        # 初始化Maker
        static_maker = StaticMaker(input_set_generator=scfset,run_vasp_kwargs={"job_type":self.jobtype})
        nonscf_maker = NonSCFMaker(input_set_generator=nonscfset,run_vasp_kwargs={"job_type":self.jobtype})
        
        # 创建任务链
        static_job = static_maker.make(self.structure)
        nonscf_job = nonscf_maker.make(structure=static_job.output.structure, prev_dir=static_job.output.dir_name)
        flow = Flow([static_job, nonscf_job], name=flow_name, output=nonscf_job.output)
        flow = add_metadata_to_flow(flow, {"flow_name": flow_name})
        
        return flow
    
    def relax2scf2nonscf_flow(
        self,
        flow_name: Optional[str] = None,
        relaxset: Optional[RelaxSetGenerator] = None,
        scfset: Optional[StaticSetGenerator] = None,
        nonscfset: Optional[NonSCFSetGenerator] = None,
        ) -> Flow:
        """
        Relaxation + SCF + NonSCF calculation workflow.

        Args:
            flow_name (str): Workflow name
            relaxset (VaspInputSet): Input set for relaxation
            scfset (VaspInputSet): Input set for SCF calculation
            nonscfset (VaspInputSet): Input set for NonSCF calculation
        Returns:
            Flow: Relaxation + SCF + NonSCF workflow
        """
        
        flow_name = self._get_flow_name("relax2scf2nonscf", flow_name)
        
        if relaxset is None:
            relaxset = RelaxSetGenerator()
            warnings.warn(f'relaxset is not defined, use default MPRelaxSet()')
        if scfset is None:
            scfset = StaticSetGenerator()
            warnings.warn(f'relaxset is not defined, use default MPStaticSet()')
        if nonscfset is None:
            nonscfset = NonSCFSetGenerator()
            warnings.warn(f'staticset is not defined, use default MPNonSCFSet()')
        # 初始化Maker
        relax_maker = RelaxMaker(input_set_generator=relaxset,run_vasp_kwargs={"job_type":self.jobtype})
        static_maker = StaticMaker(input_set_generator=scfset,run_vasp_kwargs={"job_type":self.jobtype})
        nonscf_maker = NonSCFMaker(input_set_generator=nonscfset,run_vasp_kwargs={"job_type":self.jobtype})
        
        # 创建任务链
        relax_job = relax_maker.make(self.structure)
        static_job = static_maker.make(structure=relax_job.output.structure, prev_dir=relax_job.output.dir_name)
        nonscf_job = nonscf_maker.make(structure=static_job.output.structure, prev_dir=static_job.output.dir_name)
        flow = Flow([relax_job, static_job, nonscf_job], name=flow_name, output=[nonscf_job.output])
        flow = add_metadata_to_flow(flow, {"flow_name": flow_name})
        
        return flow

    def adsorption_flow(
        self,
        boundary_frac: float,
        substrate_fix_frac: float,
        flow_name: Optional[str] = None,
        slabrelax_set: Optional[RelaxSetGenerator] = None,
        molrelax_set: Optional[RelaxSetGenerator] = None,
        slabstatic_set: Optional[StaticSetGenerator] = None,
        molstatic_set: Optional[StaticSetGenerator] = None,
        cubic_size: int = 20,
        ) -> Flow:
        """
        Generate an adsorption energy workflow. This workflow is suitable for models with clear interface boundaries.

        Parameters:
            boundary_frac (float): Fractional position of the interface between the substrate and the adsorbed molecule (0~1).
            substrate_fix_frac (float): Fraction of atoms in the substrate to be fixed during relaxation.
            flow_name (str, optional): Name of the workflow. If None, a default name will be generated.
            slabrelax_set (RelaxSetGenerator, optional): Input set generator for slab relaxation. If None, a default will be used.
            molrelax_set (RelaxSetGenerator, optional): Input set generator for molecule relaxation. If None, a default will be used.
            slabstatic_set (StaticSetGenerator, optional): Input set generator for slab static calculation. If None, a default will be used.
            molstatic_set (StaticSetGenerator, optional): Input set generator for molecule static calculation. If None, a default will be used.
            cubic_size (int, optional): The size of the cubic box for molecule extraction. Default is 20.

        Returns:
            Flow: The constructed adsorption energy workflow, including relaxation and static calculations for the adsorbed system, the clean slab, and the molecule.
        """
        flow_name = self._get_flow_name("adsorption_flow", flow_name)
        base_structure, mol_structure = SlabModify.split_mol(self.structure, boundary_frac, cubic_size)
        fixed_base_strucutre = SlabModify.fix_slab(base_structure, substrate_fix_frac)
        fixed_abs_structure = SlabModify.fix_slab(self.structure, substrate_fix_frac)
        
        if slabrelax_set is None:
            slabrelax_set = RelaxSetGenerator(user_incar_settings={"ISIF": 2})
            warnings.warn(f'slabrelax_set is not defined, use default MPRelaxSet()')
        if molrelax_set is None:
            molrelax_set = RelaxSetGenerator(user_incar_settings={"ISIF": 2})
            warnings.warn(f'molrelax_set is not defined, use default MPRelaxSet()')
        if slabstatic_set is None:
            slabstatic_set = StaticSetGenerator()
            warnings.warn(f'slabstatic_set is not defined, use default MPStaticSet()')
        if molstatic_set is None:
            molstatic_set = StaticSetGenerator()
            warnings.warn(f'molstatic_set is not defined, use default MPStaticSet()')
        
        # 初始化Maker
        absrelax_maker = RelaxMaker(
            name='abs_relax',
            input_set_generator=slabrelax_set,
            run_vasp_kwargs={"job_type":self.jobtype}
            )
        absstatic_maker = StaticMaker(
            name='abs_static',
            input_set_generator=slabstatic_set,
            run_vasp_kwargs={"job_type":self.jobtype}
            )
        slabrelax_maker = RelaxMaker(
            name='slab_relax',
            input_set_generator=slabrelax_set,
            run_vasp_kwargs={"job_type":self.jobtype}
            )
        slabstatic_maker = StaticMaker(
            name='slab_static',
            input_set_generator=slabstatic_set,
            run_vasp_kwargs={"job_type":self.jobtype}
            )
        molrelax_maker = RelaxMaker(
            name='mol_relax',
            input_set_generator=molrelax_set,
            run_vasp_kwargs={"job_type":self.jobtype}
            )
        molstatic_maker = StaticMaker(
            name='mol_static',
            input_set_generator=molstatic_set,
            run_vasp_kwargs={"job_type":self.jobtype}
            )
        
        # 创建任务链
        abs_relax_job = absrelax_maker.make(fixed_abs_structure)
        abs_static_job = absstatic_maker.make(structure=abs_relax_job.output.structure, prev_dir=abs_relax_job.output.dir_name)
        base_relax_job = slabrelax_maker.make(fixed_base_strucutre)
        base_static_job = slabstatic_maker.make(structure=base_relax_job.output.structure, prev_dir=base_relax_job.output.dir_name)
        mol_relax_job = molrelax_maker.make(mol_structure)
        mol_static_job = molstatic_maker.make(structure=mol_relax_job.output.structure, prev_dir=mol_relax_job.output.dir_name)
        
        flow = Flow([abs_relax_job, abs_static_job, base_relax_job, base_static_job, mol_relax_job, mol_static_job], name=flow_name)
        flow = add_metadata_to_flow(flow, {"flow_name": flow_name})
        
        return flow
    
    def adhesive_work_flow(
        self,
        boundary_frac: float,
        film_fix_frac: float = None,
        substrate_fix_frac: float = None,
        double_fix_heter: bool = False,
        flow_name: Optional[str] = None,
        slabrelax_set: Optional[RelaxSetGenerator] = None,
        slabstatic_set: Optional[StaticSetGenerator] = None,
        ):
        """
        Construct an adhesion energy workflow for heterostructures with a clear interface boundary.

        This workflow splits the input structure into film and substrate parts at the specified interface position,
        applies selective dynamics (fixing atoms) to each part and the combined heterostructure, and then performs
        relaxation and static calculations for each. Finally, it computes the adhesion energy per unit area.

        Parameters:
            boundary_frac (float): Fractional z-coordinate (0~1) indicating the interface between film and substrate.
            film_fix_frac (float, optional): Fractional z-coordinate threshold for fixing atoms in the film during relaxation.
                If None, defaults to half the z-span of the film.
            substrate_fix_frac (float, optional): Fractional z-coordinate threshold for fixing atoms in the substrate during relaxation.
                If None, defaults to half the z-span of the substrate.
            double_fix_heter (bool, optional): If True, fix both film and substrate regions in the heterostructure during relaxation.
                If False, only fix the substrate region. Default is False.
            flow_name (str, optional): Name of the workflow. If None, a default name will be generated.
            slabrelax_set (RelaxSetGenerator, optional): Input set generator for slab relaxation. If None, a default will be used.
            slabstatic_set (StaticSetGenerator, optional): Input set generator for slab static calculation. If None, a default will be used.

        Returns:
            Flow: The constructed jobflow Flow object representing the adhesion energy workflow.
        """
        flow_name = self._get_flow_name("adhesive_work_flow", flow_name)
        film, substrate = SlabModify.split_heter(self.structure, boundary_frac)
        
        if film_fix_frac is None:
            film_fix_frac = 0.5*self._parse_z_frac_coords(film)
        if substrate_fix_frac is None:
            substrate_fix_frac = 0.5*self._parse_z_frac_coords(substrate)
        
        fixed_film_structure = SlabModify.fix_slab(film, max_fix_ratio=film_fix_frac)
        fixed_substrate_structure = SlabModify.fix_slab(substrate, min_fix_ratio=substrate_fix_frac)
        if double_fix_heter:
            fixed_AB_structure = SlabModify.fix_slab(self.structure, min_fix_ratio=substrate_fix_frac, max_fix_ratio=film_fix_frac)
        else:
            fixed_AB_structure = SlabModify.fix_slab(self.structure, min_fix_ratio=substrate_fix_frac)
        
        if slabrelax_set is None:
            slabrelax_set = RelaxSetGenerator(user_incar_settings={"ISIF": 2})
            warnings.warn(f'slabrelax_set is not defined, use default MPRelaxSet()')
        if slabstatic_set is None:
            slabstatic_set = StaticSetGenerator()
            warnings.warn(f'slabstatic_set is not defined, use default MPStaticSet()')
        
        # Initialize Makers
        slabA_relax_maker = RelaxMaker(
            name='slabA_relax',
            input_set_generator=slabrelax_set,
            run_vasp_kwargs={"job_type":self.jobtype}
            )
        slabB_relax_maker = RelaxMaker(
            name='slabB_relax',
            input_set_generator=slabrelax_set,
            run_vasp_kwargs={"job_type":self.jobtype}
            )
        AB_relax_maker = RelaxMaker(
            name='AB_relax',
            input_set_generator=slabrelax_set,
            run_vasp_kwargs={"job_type":self.jobtype}
            )
        slabA_static_maker = StaticMaker(
            name='slabA_static',
            input_set_generator=slabstatic_set,
            run_vasp_kwargs={"job_type":self.jobtype}
            )
        slabB_static_maker = StaticMaker(
            name='slabB_static',
            input_set_generator=slabstatic_set,
            run_vasp_kwargs={"job_type":self.jobtype}
            )
        AB_static_maker = StaticMaker(
            name='AB_static',
            input_set_generator=slabstatic_set,
            run_vasp_kwargs={"job_type":self.jobtype}
            )
        
        
        # Get the cross-sectional area of the AB structure
        slab_area = self.structure.surface_area
        
        # Create job chain
        slabA_relax_job = slabA_relax_maker.make(fixed_film_structure)
        slabB_relax_job = slabB_relax_maker.make(fixed_substrate_structure)
        AB_relax_job = AB_relax_maker.make(fixed_AB_structure)
        slabA_static_job = slabA_static_maker.make(structure=slabA_relax_job.output.structure, prev_dir=slabA_relax_job.output.dir_name)
        slabB_static_job = slabB_static_maker.make(structure=slabB_relax_job.output.structure, prev_dir=slabB_relax_job.output.dir_name)
        AB_static_job = AB_static_maker.make(structure=AB_relax_job.output.structure, prev_dir=AB_relax_job.output.dir_name)
        
        Eadhesive_job = self._calc_Eadhesive(
            AB_static_job.output.energy, 
            slabA_static_job.output.energy, 
            slabB_static_job.output.energy,
            slab_area)
        
        flow = Flow([slabA_relax_job, slabB_relax_job, AB_relax_job, slabA_static_job, slabB_static_job, AB_static_job, Eadhesive_job],
                    name=flow_name,
                    output=[Eadhesive_job.output])
        
        return flow
    
    def _parse_z_frac_coords(
        self,
        slab: Slab | Structure,
        ):
        z_coords = [site.frac_coords[2] for site in slab.sites]
        max_z, min_z = max(z_coords), min(z_coords)
        return (min_z + max_z)
    
    @job
    def _calc_Eadhesive(self, E_AB: float, E_A: float, E_B: float, area: float) -> float:
        """Return adhesion energy per area with the 1/2 factor for two interfaces."""
        return (E_AB - E_A - E_B) / (2.0 * area)
        
class SubWF:
    """
    Workflow submission utility class.

    Provides two submission methods: jobflow_remote and Fireworks.
    """
    @staticmethod
    def sub_jfremote(
        flow:Flow,
        project: str | None = None,
        worker: str | None = None,
        resources: dict | None = None,
        *args, **kwargs):
        """
        Submit the workflow using jobflow_remote.

        Args:
            flow (Flow): Workflow object
            project (str): Project name
            worker (str): Worker name
            resources (dict): Resource configuration
        """
        submit_flow(flow=flow, worker=worker, project=project, resources=resources, *args, **kwargs)
    
    @staticmethod    
    def sub_fireworks(flow:Flow):
        """
        Submit the workflow using Fireworks.

        Args:
            flow (Flow): Workflow object
        """
        wf  = flow_to_workflow(flow)
        lp = LaunchPad.auto_load()
        lp.add_wf(wf)