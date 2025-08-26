"""
This module is used to query VASP calculation results via jobflow_remote, including energy, structure, band structure, density of states, etc.
"""
from pymatgen.core.structure import Structure
from pymatgen.electronic_structure.dos import CompleteDos
from pymatgen.electronic_structure.bandstructure import BandStructureSymmLine
from jobflow_remote import get_jobstore
from pprint import pprint
import pandas as pd
import numpy as np

class querydata:
    """
    Utility class for querying the Jobflow remote database.
    Supports querying energy, structure, band structure, density of states, and more.
    """

    def __init__(self,jf_project):
        """
        Initialize the query utility.
        
        Args:
            jf_project (str): Jobflow project name.
        """
        self.jf_project = jf_project
        self.store = None
        self._set_JOBSTORE()
        
    def _set_JOBSTORE(self):
        """
        Set up the JobStore connection.
        
        Returns:
            store: JobStore object.
        """
        self.store = get_jobstore(project_name=self.jf_project)
        return self.store
    
    def _parse_results(self, results, properties: list):
        """
        Dynamically parse query results.

        Args:
            results (dict): Result dictionary returned by the query.
            properties (list): List of property paths, e.g. ['output.output.energy'].

        Returns:
            dict: Dictionary containing property paths and their corresponding values.

        Raises:
            KeyError: If a key in the path does not exist.
        """
        parsed_results = {}
        for prop in properties:
            try:
                keys = prop.split('.')
                value = results
                for key in keys:
                    value = value[key]
                parsed_results[prop] = value
            except KeyError as e:
                raise KeyError(f"KeyError in parsing property '{prop}': {e}. Results: {results}")
        return parsed_results

    def query(self, 
              criteria: dict = None,
              properties: dict | list = None, 
              load: bool = False):
        """
        Custom data query.

        Args:
            criteria (dict): Query criteria.
            properties (dict | list): Properties to query.
            load (bool): Whether to load all objects.
        Returns:
            dict: Query result.
        Raises:
            ConnectionError: If JobStore is not connected.
            ValueError: If no result is found.
        """
        if not self.store:
            raise ConnectionError("JobStore is not set")
        self.store.connect(force_reset=True)
        results = self.store.query_one(
            criteria=criteria,
            properties=properties,
            load=load,
        )
        if not results:
            raise ValueError(f"No document found for criteria: {criteria}")
        
        return results

    def query_final_energy(self, 
                     criteria: dict = None,
                     properties: dict | list = ['output.output.energy'],):
        """
        Query the final energy.

        Args:
            criteria (dict): Query criteria.
            properties (dict | list): Properties to query, default ['output.output.energy'].
        Returns:
            float: Energy value.
        """
        results = self.query(criteria=criteria, properties=properties)
        parsed_results = self._parse_results(results, properties)
        return parsed_results[properties[0]]

    def query_final_structure(self, 
                              criteria: dict = None,
                              properties: dict | list = ['output.output.structure'],) -> Structure:
        """
        Query the final structure.

        Args:
            criteria (dict): Query criteria.
            properties (dict | list): Properties to query, default ['output.output.structure'].
        Returns:
            Structure: pymatgen Structure object.
        """
        results = self.query(criteria=criteria, properties=properties)
        parsed_results = self._parse_results(results, properties)
        return Structure.from_dict(parsed_results[properties[0]])
        
    def query_dos(self,
                  criteria: dict = {"name": "non-scf uniform"},
                  properties: dict | list = ['output.vasp_objects.dos'],
                  load: bool = True) -> CompleteDos:
        """
        Query the density of states (DOS).

        Args:
            criteria (dict): Query criteria, default {"name": "non-scf uniform"}.
            properties (dict | list): Properties to query, default ['output.vasp_objects.dos'].
            load (bool): Whether to load all objects.
        Returns:
            CompleteDos: pymatgen CompleteDos object.
        """
        results = self.query(criteria=criteria, properties=properties, load=load)
        parsed_results = self._parse_results(results, properties)
        return CompleteDos.from_dict(parsed_results[properties[0]])
        
    def query_bandstructure(self,
                            criteria: dict = {"name": "non-scf line"},
                            properties: dict | list = ['output.vasp_objects.bandstructure'],
                            load: bool = True) -> BandStructureSymmLine:
        """
        Query the band structure.

        Args:
            criteria (dict): Query criteria, default {"name": "non-scf line"}.
            properties (dict | list): Properties to query, default ['output.vasp_objects.bandstructure'].
            load (bool): Whether to load all objects.
        Returns:
            BandStructureSymmLine: pymatgen BandStructureSymmLine object.
        """
        results = self.query(criteria=criteria, properties=properties, load=load)
        results = self.query(criteria=criteria, properties=properties, load=load)
        parsed_results = self._parse_results(results, properties)
        return BandStructureSymmLine.from_dict(parsed_results[properties[0]])
        
    # def query_elastic_tensor(self,
    #                          criteria: dict = {"name":"fit_elastic_tensor"},
    #                          properties: dict | list = ['output.elastic_tensor.raw'],
    #                          load: bool = False):
    #     results = self.query(criteria=criteria, properties=properties, load=load)
    #     parsed_results = self._parse_results(results, properties)
    #     return pd.DataFrame(parsed_results[properties[0]])
        
    # def query_elastic_properties(self,
    #                              criteria: dict = {"name":"fit_elastic_tensor"},
    #                              properties: dict | list = ['output.derived_properties'],
    #                              load: bool = False):
    #     results = self.query(criteria=criteria, properties=properties, load=load)
    #     parsed_results = self._parse_results(results, properties)
    #     elastic_properties = parsed_results[properties[0]]
    #     df_elastic_properties = pd.DataFrame(list(elastic_properties.items()), columns=["Property", "Value"])
    #     return df_elastic_properties