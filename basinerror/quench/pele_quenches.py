"""
Interface for pele quenches to find the minima.
"""

from .base_minima_finder import BaseMinimaFinder




class PeleMinimaFinder(BaseMinimaFinder):
    """
    Interface to use pele quenches to find minima
    """
    def __init__(self, base_folder, paramstring, quench_function):
        """
        Initialize the PeleMinimaFinder
        """
        BaseMinimaFinder.__init__(self, base_folder, paramstring)
        self.quench_function = quench_function
    def find_minimum_production(self, coords, **kwargs):
        """         Run the quench function on the given coordinates
        Production version treats raised exceptions as failures.

        Parameters
        ----------
        coords : np.ndarray
            starting coordinates for which we want to find the minimum

        Returns
        -------
        ret : Enum
            Contains information about the completed run.
        """        
        try:
            ret = self.quench_function(coords, self.potential, **kwargs)
        except:
            ret = None
        return ret
    def find_minimum(self, coords, **kwargs):
        """
        Run the quench function on the given coordinates

        Parameters
        ----------
        coords : np.ndarray
            starting coordinates for which we want to find the minimum

        Returns
        -------
        ret : Enum
            Contains information about the completed run.
        """
        ret = self.quench_function(coords, self.potential, **kwargs)
        return ret
