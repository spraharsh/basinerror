"""
Base class for Minima finders. Provides a common interface for all quenches.
to be used in the quench.py module.
"""

from .utils.params import get_system_parameters_from_folder, setup_inverse_power


class BaseMinimaFinder(object):
    """
    Base class for Minima finders. Provides a common interface for all Minima finding algorithms.

    """

    def __init__(self, base_folder, paramstring):
        """
        Initialize the Minima finder.

        Parameters
        ----------
        config : dict
            the quench configuration
        """
        self.base_folder = base_folder
        self.SysPars = get_system_parameters_from_folder(self.base_folder)
        # this is currently setup for finding minima for inverse power potentials
        self.potential = setup_inverse_power(self.SysPars)


    def find_minimum_production(self, initial_coords):
        """
        Find the minimum of the system. 
        This is a template for a production run, assuming
        assuming that any errors are failed runs.

        Parameters
        ----------
        initial_coords : np.ndarray
            starting coordinates

        Raises
        ------
        NotImplementedError
        """
        raise NotImplementedError('quench method not implemented')

    def get_potential(self):
        """ 
        Return the potential object

        Returns
        -------
        potential : potential.Potential
            the potential object
        """
        return self.potential
