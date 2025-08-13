from pymatgen.entries.computed_entries import ComputedEntry
from pymatgen.analysis.phase_diagram import *
from pymatgen.io.vasp.outputs import Vasprun
from monty.serialization import dumpfn
import os


class ProcessEntries:
    '''
    Traverse all calculation tasks in the current path and generate computed entries.
    
    Note: Make sure all calculation tasks in the current path are available.
    '''
    def __init__(self, input_path, output_path: None):
        '''
        Initialize the ProcessEntries instance.
        
        Args:
            input_path: Working directory
            output_path: Path to entries.json
        Returns:
            entries.json
        '''
        self.input_path = input_path
        self.output_path = output_path
        self.gen_entries()
        
    def parse_calc_dirs(self):
        '''
        Search all calculation directories and find folders containing vasprun.xml.
        
        Returns:
            calc_dirs: List[str], all directory paths containing vasprun.xml
        '''
        calc_dirs = []
        
        for root, dirs, files in os.walk(self.input_path):
            if 'vasprun.xml' in files:
                calc_dirs.append(root)
            else:
                print(f"No vasprun.xml in {root}, pls check")
        return calc_dirs
    
    def gen_entries(self):
        '''
        Generate computed entries, with entry_id as the calculation directory name by default.
        
        If output_path is None, do not save the json file.
        
        Returns:
            entries: List[ComputedEntry]
            entries.json
        '''
        
        entries = []
        calc_dirs = self.parse_calc_dirs()
        for calc_dir in calc_dirs:
            try:
                vasprun = Vasprun(os.path.join(calc_dir, 'vasprun.xml'))
                entry_id = os.path.basename(calc_dir)
                entry = vasprun.get_computed_entry(entry_id=entry_id)
                entries.append(entry)
                print(f"Successfully created entry for {calc_dir}")
            except Exception as e:
                print(f"Error in {calc_dir}: {e}")
        if self.output_path:
            dumpfn(entries, self.output_path)
            print(f"Entries saved to {self.output_path}")
        
        self.entries = entries
        return self.entries
    
    @classmethod
    def parse_stable_entry_id(
        cls, 
        entries: list[ComputedEntry], 
        mode = 'simple', # 'simple' | 'compound'
        elements: list[str] = None,
        terminal_compositions: list[Composition] = None,
        ):
        '''
        Parse the entry_id of stable entries.
        
        Args:
            entries: List[ComputedEntry]
            mode: 'simple' | 'compound'
                'simple': ordinary phase diagram, elements must be specified
                'compound': compound phase diagram, terminal_compositions must be specified
            elements: List of elements (required if mode='simple')
            terminal_compositions: Terminal compositions (required if mode='compound')
        Returns:
            stable_entry_id: List[str]
        Raises:
            ValueError: If required parameters are missing or mode is not supported
        '''
        
        stable_entry_id = []
        if mode == 'simple':
            if not elements:
                raise ValueError("Please specify elements")
            PD = PhaseDiagram(entries, elements=elements)
        elif mode == 'compound':
            if not terminal_compositions:
                raise ValueError("Please specify terminal_compositions")
            PD = CompoundPhaseDiagram(entries, terminal_compositions=terminal_compositions, normalize_terminal_compositions=False)
        else:
            raise ValueError(f"Unsupported mode: {mode}. Use 'simple' or 'compound'.")
        stable_entries = PD.stable_entries
        entry_map = {(entry.composition.reduced_formula, entry.energy_per_atom): entry for entry in entries}
        for stable_entry in stable_entries:
            if stable_entry.formula in entry_map:
                stable_entry_id.append(entry_map[stable_entry.formula].entry_id)
        
        return stable_entry_id
        