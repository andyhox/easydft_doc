from pymatgen.io.lobster.inputs import Lobsterin
from typing import Optional, Union
from os import PathLike
from pathlib import Path
import warnings


class Lobster:

    @classmethod
    def gen_lobsterin_from_vasp_inp(
        cls,
        option: Optional[str] = "onlycohpcoop",
        INCAR_file: Optional[PathLike] = None,
        POSCAR_file: Optional[PathLike] = None,
        POTCAR_file: Optional[PathLike] = None,
        calc_dir: Optional[PathLike] = None,
        save_dir: Optional[PathLike] = None,
        dict_for_basis: Optional[dict] = None,
        overwritedict: Optional[dict] = None,
    ):
        """
        Generate a lobsterin file based on VASP input files or calculation directory.

        Parameters
        ----------
        option : str
            Lobsterin generation option, e.g. "onlycohpcoop", "standard".
        INCAR_file, POSCAR_file, POTCAR_file : PathLike
            Paths to VASP INCAR, POSCAR, POTCAR files.
        calc_dir : PathLike
            Directory containing VASP input files.
        save_dir : PathLike
            Directory to save the generated lobsterin.
        """
        INCAR_file = Path(INCAR_file) if INCAR_file else None
        POSCAR_file = Path(POSCAR_file) if POSCAR_file else None
        POTCAR_file = Path(POTCAR_file) if POTCAR_file else None
        calc_dir = Path(calc_dir) if calc_dir else None
        save_dir = Path(save_dir) if save_dir else None

        lobster = None

        # case 1: all three files are provided
        if INCAR_file and POSCAR_file and POTCAR_file:
            lobster = Lobsterin.standard_calculations_from_vasp_files(
                INCAR_input=INCAR_file,
                POSCAR_input=POSCAR_file,
                POTCAR_input=POTCAR_file,
                option=option,
                dict_for_basis=dict_for_basis
            )
        # case 2: try to parse from calc_dir
        elif calc_dir:
            incar_path = calc_dir / "INCAR"
            poscar_path = calc_dir / "POSCAR"
            potcar_path = calc_dir / "POTCAR"
            if not (incar_path.exists() and poscar_path.exists() and potcar_path.exists()):
                raise FileNotFoundError(
                    f"calc_dir does not contain all required files: INCAR, POSCAR, POTCAR"
                )
            lobster = Lobsterin.standard_calculations_from_vasp_files(
                INCAR_input=incar_path,
                POSCAR_input=poscar_path,
                POTCAR_input=potcar_path,
                option=option,
            )
        # case 3: error
        else:
            raise ValueError(
                "You must provide either (INCAR_file, POSCAR_file, POTCAR_file) "
                "or calc_dir containing them."
            )

        # write lobsterin
        if save_dir:
            save_dir.mkdir(parents=True, exist_ok=True)
            lobster.write_lobsterin(path=save_dir / "lobsterin", overwritedict=overwritedict)
        else:
            warnings.warn(
                "save_dir not provided, the lobsterin file will not be saved!"
            )

        return lobster
