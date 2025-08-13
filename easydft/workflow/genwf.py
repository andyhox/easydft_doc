"""
This module is used to generate VASP workflows, supporting various tasks such as structure optimization, static calculation, non-self-consistent calculation, adsorption energy, and adhesion energy.
"""
from pymatgen.core import Structure, Lattice
from pymatgen.io.vasp.sets import MPRelaxSet, MPStaticSet, MPNonSCFSet, VaspInputSet
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
import numpy as np
import warnings

class GenWF:
    """
    VASP workflow generator

    Used to generate different types of VASP calculation workflows based on the input structure, including double relaxation, relaxation+static, SCF+NonSCF, adsorption energy, adhesion energy, etc.
    """
    
    def __init__(
        self, 
        structure: Structure,
        use_custodian: bool = True,
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
        date_str = datetime.now().strftime("%Y%m%d")
        return f"{formula}_{date_str}"

    def _get_flow_name(self, suffix: str, flow_name: str = None) -> str:
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
        flow_name: str = None,
        job_name: str = None,
        set1: VaspInputSet = None,
        set2: VaspInputSet = None,
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
            set1 = MPRelaxSet()
            warnings.warn(f'set1 is not defined, use default MPRelaxSet()')
        if set2 is None:
            set2 = MPRelaxSet()
            warnings.warn(f'set2 is not defined, use default MPRelaxSet()')
        relax1_maker = RelaxMaker(input_set_generator=set1,run_vasp_kwargs={"job_type":self.jobtype})
        relax2_maker = RelaxMaker(input_set_generator=set2,run_vasp_kwargs={"job_type":self.jobtype})
        doublerelax_maker = DoubleRelaxMaker(relax1_maker, relax2_maker)
        
        doublerelax_job = doublerelax_maker.make(self.structure)
        flow = Flow(doublerelax_job, name=flow_name)
        flow = add_metadata_to_flow(flow, {"flow_name": flow_name})
        return flow
    
    def relax2scf_flow(
        self,
        flow_name: str = None,
        relaxset: VaspInputSet = None,
        scfset: VaspInputSet = None,
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
            relaxset = MPRelaxSet()
            warnings.warn(f'relaxset is not defined, use default MPRelaxSet()')
        if scfset is None:
            scfset = MPStaticSet()
            warnings.warn(f'staticset is not defined, use default MPStaticSet()')
        relax_maker = RelaxMaker(input_set_generator=relaxset,run_vasp_kwargs={"job_type":self.jobtype})
        static_maker = StaticMaker(input_set_generator=scfset,run_vasp_kwargs={"job_type":self.jobtype})
        
        relax_job = relax_maker.make(self.structure)
        static_job = static_maker.make(structure=relax_job.output.structure, prev_dir=relax_job.output.dir_name)
        
        flow = Flow([relax_job, static_job], name=flow_name, output=static_job.output)
        flow = add_metadata_to_flow(flow, {"flow_name": flow_name})
        
        return flow
    
    def scf2nonscf_flow(
        self, 
        flow_name: str = None,
        scfset: VaspInputSet = None,
        nonscfset: VaspInputSet = None,
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
            scfset = MPStaticSet()
            warnings.warn(f'relaxset is not defined, use default MPRelaxSet()')
        if nonscfset is None:
            nonscfset = MPNonSCFSet()
            warnings.warn(f'staticset is not defined, use default MPStaticSet()')
        static_maker = StaticMaker(input_set_generator=scfset,run_vasp_kwargs={"job_type":self.jobtype})
        nonscf_maker = NonSCFMaker(input_set_generator=nonscfset,run_vasp_kwargs={"job_type":self.jobtype})
        
        static_job = static_maker.make(self.structure)
        nonscf_job = nonscf_maker.make(structure=static_job.output.structure, prev_dir=static_job.output.dir_name)
        flow = Flow([static_job, nonscf_job], name=flow_name, output=[nonscf_job.output])
        flow = add_metadata_to_flow(flow, {"flow_name": flow_name})
        
        return flow
    
    def relax2scf2nonscf_flow(
        self,
        flow_name: str = None,
        relaxset: VaspInputSet = None,
        scfset: VaspInputSet = None,
        nonscfset: VaspInputSet = None,
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
            relaxset = MPRelaxSet()
            warnings.warn(f'relaxset is not defined, use default MPRelaxSet()')
        if scfset is None:
            scfset = MPStaticSet()
            warnings.warn(f'relaxset is not defined, use default MPStaticSet()')
        if nonscfset is None:
            nonscfset = MPNonSCFSet()
            warnings.warn(f'staticset is not defined, use default MPNonSCFSet()')
        relax_maker = RelaxMaker(input_set_generator=relaxset,run_vasp_kwargs={"job_type":self.jobtype})
        static_maker = StaticMaker(input_set_generator=scfset,run_vasp_kwargs={"job_type":self.jobtype})
        nonscf_maker = NonSCFMaker(input_set_generator=nonscfset,run_vasp_kwargs={"job_type":self.jobtype})
        
        relax_job = relax_maker.make(self.structure)
        static_job = static_maker.make(structure=relax_job.output.structure, prev_dir=relax_job.output.dir_name)
        nonscf_job = nonscf_maker.make(structure=static_job.output.structure, prev_dir=static_job.output.dir_name)
        flow = Flow([relax_job, static_job, nonscf_job], name=flow_name, output=[nonscf_job.output])
        flow = add_metadata_to_flow(flow, {"flow_name": flow_name})
        
        return flow

    def adsorption_flow(
        self,
        boundary_frac: float,
        fix_frac: float,
        flow_name: str = None,
        slabrelax_set: VaspInputSet = None,
        molrelax_set: VaspInputSet = None,
        slabstatic_set: VaspInputSet = None,
        molstatic_set: VaspInputSet = None,
        ) -> Flow:
        """
        Adsorption energy workflow. Only supports models with clear interface boundaries.

        Args:
            boundary_frac (float): Interface fraction
            fix_frac (float): Fraction of fixed atoms in the slab
            flow_name (str): Workflow name
            slabrelax_set (VaspInputSet): Input set for slab relaxation
            molrelax_set (VaspInputSet): Input set for molecule relaxation
            slabstatic_set (VaspInputSet): Input set for slab static calculation
            molstatic_set (VaspInputSet): Input set for molecule static calculation
        Returns:
            Flow: Adsorption energy workflow
        """
        flow_name = self._get_flow_name("adsorption_flow", flow_name)
        base_structure, mol_structure = SlabModify.split_slab(self.structure, boundary_frac)
        fixed_base_strucutre = SlabModify.fix_slab(base_structure, fix_frac)
        fixed_abs_structure = SlabModify.fix_slab(self.structure, fix_frac)
        
        if slabrelax_set is None:
            slabrelax_set = MPRelaxSet(user_incar_settings={"ISIF": 2})
            warnings.warn(f'slabrelax_set is not defined, use default MPRelaxSet()')
        if molrelax_set is None:
            molrelax_set = MPRelaxSet(user_incar_settings={"ISIF": 2})
            warnings.warn(f'molrelax_set is not defined, use default MPRelaxSet()')
        if slabstatic_set is None:
            slabstatic_set = MPStaticSet()
            warnings.warn(f'slabstatic_set is not defined, use default MPStaticSet()')
        if molstatic_set is None:
            molstatic_set = MPStaticSet()
            warnings.warn(f'molstatic_set is not defined, use default MPStaticSet()')
        
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
        slabA_structure:Structure,
        slabB_structure:Structure,
        fix_frac: float,
        flow_name: str = None,
        slabrelax_set: VaspInputSet = None,
        slabstatic_set: VaspInputSet = None,
        ):
        """
        Adhesion energy workflow. Only supports models with clear interface boundaries.

        Args:
            slabA_structure (Structure): Structure of slab A
            slabB_structure (Structure): Structure of slab B
            fix_frac (float): Fraction of fixed atoms
            flow_name (str): Workflow name
            slabrelax_set (VaspInputSet): Input set for slab relaxation
            slabstatic_set (VaspInputSet): Input set for slab static calculation
        Returns:
            Flow: Adhesion energy workflow
        """
        flow_name = self._get_flow_name("adhesive_work_flow", flow_name)
        fixed_slabA_structure = SlabModify.fix_slab(slabA_structure, fix_frac)
        fixed_slabB_structure = SlabModify.fix_slab(slabB_structure, fix_frac)
        fixed_AB_structure = SlabModify.fix_slab(self.structure, fix_frac)
        
        if slabrelax_set is None:
            slabrelax_set = MPRelaxSet(user_incar_settings={"ISIF": 2})
            warnings.warn(f'slabrelax_set is not defined, use default MPRelaxSet()')
        if slabstatic_set is None:
            slabstatic_set = MPStaticSet()
            warnings.warn(f'slabstatic_set is not defined, use default MPStaticSet()')
        
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
        
        
        slab_area = SlabModify.slab_area(self.structure)
        
        slabA_relax_job = slabA_relax_maker.make(fixed_slabA_structure)
        slabB_relax_job = slabB_relax_maker.make(fixed_slabB_structure)
        AB_relax_job = AB_relax_maker.make(fixed_AB_structure)
        slabA_static_job = slabA_static_maker.make(structure=slabA_relax_job.output.structure, prev_dir=slabA_relax_job.output.dir_name)
        slabB_static_job = slabB_static_maker.make(structure=slabB_relax_job.output.structure, prev_dir=slabB_relax_job.output.dir_name)
        AB_static_job = AB_static_maker.make(structure=AB_relax_job.output.structure, prev_dir=AB_relax_job.output.dir_name)
        Eadhesive_job = (AB_static_job.output.energy - slabA_static_job.output.energy - slabB_static_job.output.energy) / 2 * slab_area
        
        flow = Flow([slabA_relax_job, slabB_relax_job, AB_relax_job, slabA_static_job, slabB_static_job, AB_static_job, Eadhesive_job],
                    name=flow_name,
                    output=[Eadhesive_job.output])
        
        return flow
        
class SubWF:
    """
    Workflow submission utility class.

    Provides two submission methods: jobflow_remote and Fireworks.
    """
    @staticmethod
    def sub_jfremote(flow:Flow,
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