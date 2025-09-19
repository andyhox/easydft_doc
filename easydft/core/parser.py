from pymatgen.io.vasp.outputs import Vasprun, Xdatcar
from pymatgen.io.vasp.inputs import Kpoints
import os
from abc import ABC, abstractmethod
from typing import Any

class Parser(ABC):
    def __init__(self, file_path):
        self.file_path = file_path
        self._validate_file()

    def _validate_file(self):
        if not os.path.exists(self.file_path):
            raise FileNotFoundError(f"File not found: {self.file_path}")

    @abstractmethod
    def parse(self) -> Any:
        pass
    
class VasprunParser(Parser):
    def __init__(self, file_path):
        super().__init__(file_path)
        self.parse()

    def parse(self) -> Vasprun:
        vasprun = Vasprun(self.file_path, parse_projected_eigen=True)
        self.vasprun = vasprun
        return self.vasprun

class XdatcarParser(Parser):
    def __init__(self, file_path):
        super().__init__(file_path)
        self.parse()

    def parse(self) -> Xdatcar:
        xdatcar = Xdatcar(self.file_path)
        self.xdatcar = xdatcar
        return self.xdatcar
    
    @property
    def structures(self):
        return self.xdatcar.structures
    
class KpointsParser(Parser):
    def __init__(self, file_path):
        super().__init__(file_path)
        self.parse()

    def parse(self) -> Kpoints:
        kpoints = Kpoints(self.file_path)
        self.kpoints = kpoints
        return self.kpoints
    
    @property
    def labels(self):
        return self.kpoints.labels