#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""


Molecular dynamics data processing toolkit.
Features:
1. Read VASP XDATCAR and OSZICAR files
2. Convert to PDB format trajectory
3. Plot energy and temperature curves
4. Handle periodic boundary conditions
"""

import os
import shutil
import numpy as np
import matplotlib as mpl
import math
from matplotlib import pyplot as plt
from copy import deepcopy

mpl.use('Agg')  # No GUI mode


class EnergyTempReader:
    """
    Read energy and temperature data from OSZICAR file.
    """
    
    def __init__(self, oszicar_path="OSZICAR"):
        """
        Initialize the energy and temperature reader.
        
        Args:
            oszicar_path (str): Path to the OSZICAR file, default "OSZICAR"
        """
        self.temp = []
        self.energy = []
        self.oszicar_path = oszicar_path
        self._read_oszicar()

    def _read_oszicar(self):
        """
        Read the content of the OSZICAR file and extract temperature and energy data.
        
        Raises:
            IOError: Raised if the file does not exist.
        """
        if not os.path.exists(self.oszicar_path):
            raise IOError(f'OSZICAR file does not exist: {self.oszicar_path}')
        
        print('Reading energy and temperature data from OSZICAR...')
        with open(self.oszicar_path, 'r') as f:
            for line in f:
                if "E0=" in line:
                    parts = line.split()
                    self.temp.append(float(parts[2]))  # Temperature
                    self.energy.append(float(parts[6]))  


class XDATCARParser(EnergyTempReader):
    """
    XDATCAR file parser.
    Inherits from EnergyTempReader.
    """
    
    def __init__(self, xdatcar_path="XDATCAR", oszicar_path="OSZICAR"):
        """
        Initialize the XDATCAR parser.
        
        Args:
            xdatcar_path (str): Path to the XDATCAR file, default "XDATCAR"
            oszicar_path (str): Path to the OSZICAR file, default "OSZICAR"
        """
        super().__init__(oszicar_path)
        self.xdatcar_path = xdatcar_path
        self._init_variables()
        self._parse_xdatcar()

    def _init_variables(self):
        """
        Initialize internal variables.
        """
        self.lattice = np.zeros((3, 3))
        self.NPT = False
        self.frames = 0
        self._timestep = 1.0  # Default time step of 1fs
        self.alpha = 0.0
        self.beta = 0.0
        self.gamma = 0.0
        self.current_frame = 0
        self._lattice = np.zeros(3)
        self.lowrange = 0
        self.uprange = 0

    def _parse_xdatcar(self):
        """
        Parse the XDATCAR file to get the number of frames, NPT status, lattice information, etc.
        
        Raises:
            IOError: Raised if the file does not exist.
            ValueError: Raised if the file format is incorrect.
        """
        print('Parsing XDATCAR file...')
        if not os.path.exists(self.xdatcar_path):
            raise IOError(f'XDATCAR file does not exist: {self.xdatcar_path}')

        with open(self.xdatcar_path, 'r') as f:
            # Calculate the total number of frames
            self.frames = sum(1 for line in f if "Direct" in line)
            f.seek(0)
            
            # Check if it's an NPT simulation
            title = f.readline().strip()
            f.seek(0)
            self.NPT = not any("Direct configuration" in line for line in f)
            f.seek(0)
            
            # Read lattice information
            if not self.NPT:
                self._read_lattice_info(f)

        print(f'Total frames: {self.frames}, NPT simulation: {self.NPT}')
        self.uprange = self.frames - 1

    def _read_lattice_info(self, file_handle):
        """
        Read lattice information, including lattice constants, angles, element types, and quantities.
        
        Args:
            file_handle: File handle
        
        Raises:
            ValueError: Raised if the file format does not meet VASP 5.x requirements.
        """
        self.title = file_handle.readline().strip()
        self.scaling_factor = float(file_handle.readline())
        
        # Read lattice vectors
        for i in range(3):
            self.lattice[i] = np.array(
                [float(j) for j in file_handle.readline().split()])
        self.lattice *= self.scaling_factor
        
        # Calculate lattice constants
        self._lattice[0] = np.linalg.norm(self.lattice[0])
        self._lattice[1] = np.linalg.norm(self.lattice[1])
        self._lattice[2] = np.linalg.norm(self.lattice[2])
        
        # Calculate lattice angles
        self.alpha = math.degrees(math.acos(
            np.dot(self.lattice[1], self.lattice[2]) / 
            (self._lattice[1] * self._lattice[2])))
            
        self.beta = math.degrees(math.acos(
            np.dot(self.lattice[0], self.lattice[2]) / 
            (self._lattice[0] * self._lattice[2])))
            
        self.gamma = math.degrees(math.acos(
            np.dot(self.lattice[0], self.lattice[1]) / 
            (self._lattice[0] * self._lattice[1])))

        # Read element information
        self.element_list = file_handle.readline().split()
        try:
            self.element_amount = [int(j) for j in file_handle.readline().split()]
        except ValueError:
            raise ValueError('Requires VASP 5.x format XDATCAR!')
            
        # Generate element type list
        self.total_elements = []
        for elem, count in zip(self.element_list, self.element_amount):
            self.total_elements.extend([elem] * count)
            
        self.total_atom = sum(self.element_amount)
        self.atomic_position = np.zeros((self.total_atom, 3))
        self.cartesian_position = np.zeros((self.total_atom, 3))

    def set_time_range(self, begin=0, end=None):
        """
        Set the analysis time range (frame interval).
        
        Args:
            begin (int): Start frame (0-based)
            end (int): End frame (None means last frame)
        """
        if end is None:
            end = self.frames - 1
            
        self.lowrange = max(0, begin)
        self.uprange = min(self.frames - 1, end)
        print(f'Setting time range: frames {self.lowrange}-{self.uprange}')

    def parse_step(self, step):
        """
        Parse atomic positions for the specified frame.
        
        Args:
            step (int): Frame number (0-based)
        Returns:
            numpy.ndarray: Cartesian coordinates of atoms
        Raises:
            ValueError: Raised if the frame number is out of range.
        """
        if step < 0 or step >= self.frames:
            raise ValueError(f'Frame number {step} is out of range (0-{self.frames-1})')
            
        with open(self.xdatcar_path, 'r') as f:
            # Locate to the specified frame
            current_step = -1
            while current_step < step:
                line = f.readline()
                if "Direct" in line:
                    current_step += 1
                    if current_step == step:
                        break
            
            # Read atomic positions
            positions = []
            for _ in range(self.total_atom):
                line = f.readline()
                positions.append([float(x) for x in line.split()[:3]])
            
        self.atomic_position = np.array(positions)
        self.cartesian_position = np.dot(self.atomic_position, self.lattice)
        return self.cartesian_position

    def write_pdb_frame(self, frame_data, frame_number, output_file):
        """
        Write a single frame of PDB data.
        
        Args:
            frame_data (numpy.ndarray): Atomic coordinates
            frame_number (int): Frame number
            output_file (str): Output file path
        """
        with open(output_file, 'a') as f:
            # Write header information
            f.write(f"MODEL     {frame_number:>5}\n")
            f.write("REMARK   Converted from XDATCAR\n")
            f.write("REMARK   Using VASP toolkit\n")
            f.write(
                f'CRYST1{self._lattice[0]:9.3f}{self._lattice[1]:9.3f}'
                f'{self._lattice[2]:9.3f}{self.alpha:7.2f}'
                f'{self.beta:7.2f}{self.gamma:7.2f}\n')
            
            # Write atomic coordinates
            for i, (elem, pos) in enumerate(zip(self.total_elements, frame_data)):
                f.write(
                    f'ATOM  {i+1:>5} {elem:4} MOL     1    '
                    f'{pos[0]:8.3f}{pos[1]:8.3f}{pos[2]:8.3f}'
                    f'  1.00  0.00          {elem:>2}\n')
            
            f.write('TER\nENDMDL\n')

    def unwrap_pbc(self, current_pos, prev_pos):
        """
        Handle periodic boundary conditions.
        
        Args:
            current_pos (numpy.ndarray): Current frame coordinates
            prev_pos (numpy.ndarray): Previous frame coordinates
        Returns:
            tuple: (Updated previous frame coordinates, displacement vector)
        """
        diff = current_pos - prev_pos
        new_prev = current_pos.copy()
        
        # Handle periodicity in each direction
        for i in range(3):
            mask = diff[:, i] > (self._lattice[i] / 2)
            diff[mask, i] -= self._lattice[i]
            
            mask = diff[:, i] < -(self._lattice[i] / 2)
            diff[mask, i] += self._lattice[i]
        
        return new_prev, diff

    def center_molecule(self, positions, center_atom):
        """
        Center the molecule at the specified atom.
        
        Args:
            positions (numpy.ndarray): Atomic coordinates
            center_atom (int): Index of the center atom (0-based)
        Returns:
            numpy.ndarray: Centered coordinates
        Raises:
            IndexError: Raised if the specified atom does not exist.
        """
        if center_atom >= len(positions):
            raise IndexError("Specified atom does not exist!")
            
        centered = positions.copy()
        center_pos = centered[center_atom].copy()
        
        # Adjust each atom's position
        for i in range(len(centered)):
            displacement = centered[i] - center_pos
            for j in range(3):
                if displacement[j] > self._lattice[j] / 2:
                    centered[i][j] -= self._lattice[j]
                elif displacement[j] < -self._lattice[j] / 2:
                    centered[i][j] += self._lattice[j]
        
        return centered


class MDPlotter:
    """
    Molecular dynamics energy and temperature curve plotting tool.
    """
    
    def __init__(self, parser, output_prefix="", dpi=300, figsize=(5, 4)):
        """
        Initialize the plotting tool.
        
        Args:
            parser (XDATCARParser): Data parser
            output_prefix (str): Output file prefix
            dpi (int): Image resolution
            figsize (tuple): Image size (inches)
        """
        self.parser = parser
        self.output_prefix = output_prefix
        self.dpi = dpi
        self.figsize = figsize

    def plot_profiles(self):
        """
        Plot and save energy and temperature curves.
        """
        # Prepare data
        time_points = np.arange(
            self.parser.lowrange, 
            min(self.parser.uprange + 1, len(self.parser.temp))
            ) * self.parser.timestep
            
        temp_data = self.parser.temp[self.parser.lowrange:self.parser.uprange + 1]
        energy_data = self.parser.energy[self.parser.lowrange:self.parser.uprange + 1]

        label_fontsize = 8
        tick_fontsize = 8
        title_fontsize = 10
        
        # Create figure
        plt.figure(figsize=self.figsize)
        
        # Temperature subplot
        plt.subplot(211)
        plt.plot(time_points, temp_data, 'b-', linewidth=1)
        plt.xlabel('Time (fs)', fontsize=label_fontsize)
        plt.ylabel('Temperature (K)', fontsize=label_fontsize)
        plt.title('Temperature Curve', fontsize=title_fontsize)
        plt.tick_params(axis='both', which='major', labelsize=tick_fontsize)
        
        # Energy subplot
        plt.subplot(212)
        plt.plot(time_points, energy_data, 'r-', linewidth=1)
        plt.xlabel('Time (fs)', fontsize=label_fontsize)
        plt.ylabel('Energy (eV)', fontsize=label_fontsize)
        plt.title('Energy Curve', fontsize=title_fontsize)
        plt.tick_params(axis='both', which='major', labelsize=tick_fontsize)
        
        # Save figure
        plt.tight_layout()
        output_path = f'{self.output_prefix}ENERGY.png'
        plt.savefig(output_path, dpi=self.dpi, bbox_inches='tight')
        plt.close()
        print(f'Energy curve plot saved: {output_path}')
        
        # Save data files
        self._save_data(time_points, temp_data, 'Temperature.dat')
        self._save_data(time_points, energy_data, 'Energy.dat')

    def _save_data(self, time_points, data, filename):
        """
        Save data to file.
        
        Args:
            time_points (array): Time points
            data (array): Data values
            filename (str): File name
        """
        output_path = f'{self.output_prefix}{filename}'
        with open(output_path, 'w') as f:
            f.write("Time(fs)   Value\n")
            for t, d in zip(time_points, data):
                f.write(f"{t:.1f}   {d:.6f}\n")
        print(f'Data file saved: {output_path}')


class MDProcessor:
    """
    Molecular dynamics data processing controller.
    """
    
    def __init__(self, xdatcar_path="XDATCAR", oszicar_path="OSZICAR"):
        """
        Initialize the processor.
        
        Args:
            xdatcar_path (str): Path to the XDATCAR file, default "XDATCAR"
            oszicar_path (str): Path to the OSZICAR file, default "OSZICAR"
        """
        self.xdatcar_path = xdatcar_path
        self.oszicar_path = oszicar_path
        self.parser = None
        
    def process(self, **kwargs):
        """
        Main processing function.
        
        Args:
            begin (int): Start frame, default 0
            
            end (int): End frame, default None (all frames)
            
            timestep (float): Time step (fs), default 1.0
            
            output_pdb (str/bool):
                - str: Specify output path
                - True: Use default path "XDATCAR.pdb"
                - False/None: Do not output
            interval (int): Output interval, default 1
            
            unwrap_pbc (bool): Whether to handle periodic boundary, default False
            
            center_atom (int): Center atom (1-based), default None
            
            plot_profiles (bool): Whether to plot curves, default False
            
            output_prefix (str): Output prefix, default ""
        """
        # Initialize parser
        self.parser = XDATCARParser(self.xdatcar_path, self.oszicar_path)
        
        # Set parameters
        begin = kwargs.get('begin', 0)
        end = kwargs.get('end', None)
        timestep = kwargs.get('timestep', 1.0)
        self.parser.timestep = timestep
        self.parser.set_time_range(begin, end)
        
        # PDB conversion
        output_pdb = kwargs.get('output_pdb', False)
        if output_pdb:
            pdb_path = "XDATCAR.pdb" if output_pdb is True else output_pdb
            self._convert_to_pdb(
                pdb_path=pdb_path,
                interval=kwargs.get('interval', 1),
                unwrap_pbc=kwargs.get('unwrap_pbc', False),
                center_atom=kwargs.get('center_atom', None)
            )
        
        # Plot curves
        if kwargs.get('plot_profiles', False):
            self._plot_profiles(
                output_prefix=kwargs.get('output_prefix', "")
            )
    
    def _convert_to_pdb(self, pdb_path, interval, unwrap_pbc, center_atom):
        """
        Execute PDB conversion.
        
        Args:
            pdb_path (str): Output PDB file path
            interval (int): Frame output interval
            unwrap_pbc (bool): Whether to handle periodic boundary
            center_atom (int): Center atom (1-based)
        """
        if os.path.exists(pdb_path):
            shutil.copyfile(pdb_path, f"{pdb_path}.bak")
            os.remove(pdb_path)
        
        prev_pos = None
        real_pos = None
        frame_count = 0
        
        print(f'Starting conversion to PDB format, output path: {pdb_path}')
        print(f'Frame interval: {interval}, Handle PBC: {unwrap_pbc}, Center atom: {center_atom}')
        
        for frame in range(self.parser.lowrange, self.parser.uprange + 1):
            if frame % interval != 0:
                continue
                    
            # Parse the current frame
            positions = self.parser.parse_step(frame)
            
            # Handle periodic boundary conditions
            if unwrap_pbc:
                if frame == self.parser.lowrange:
                    real_pos = positions.copy()
                    if center_atom is not None:
                        real_pos = self.parser.center_molecule(
                            real_pos, center_atom - 1)
                    prev_pos = positions.copy()
                else:
                    prev_pos, disp = self.parser.unwrap_pbc(positions, prev_pos)
                    real_pos += disp
                    if center_atom is not None:
                        real_pos = self.parser.center_molecule(
                            real_pos, center_atom - 1)
                self.parser.cartesian_position = real_pos
            
            # Write the PDB frame
            self.parser.write_pdb_frame(
                self.parser.cartesian_position, 
                frame_count + 1, 
                pdb_path
            )
            frame_count += 1
        
        print(f'Conversion completed, total {frame_count} frames converted')
        print(f'Time interval: {self.parser.timestep * interval}fs')
    
    def _plot_profiles(self, output_prefix=""):
        """
        Plot energy and temperature curves.
        
        Args:
            output_prefix (str): Output file prefix
        """
        print(f'Starting to plot energy and temperature curves, prefix: {output_prefix}')
        plotter = MDPlotter(self.parser, output_prefix)
        plotter.plot_profiles()