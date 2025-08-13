from typing import List, Dict, Optional, Union
import pandas as pd
import os

class Analyzer():
    def __init__(
        self, 
        work_dir: None,
        file_path: Union[str, list, None], 
        save_data: bool = True, 
        output_path: Optional[str] = None,
        ):
        self.work_dir = work_dir
        self.file_path = file_path
        self.output_path = output_path
        self.save_data = save_data
        self._validate_input()
    def _validate_input(self):
        if not self.work_dir and not (self.poscar_file and self.locpot_file and self.outcar_file):
            raise ValueError("Either 'calc_dir' or all of 'poscar_file', 'locpot_file', 'outcar_file' must be provided.")
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