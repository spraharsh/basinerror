"""
Interface for the Julia version of Mixed Descent.
"""


from .base_minima_finder import BaseMinimaFinder
from julia import Main
# julia imports
Main.include("../../basins.jl/src/optimizer/newton.jl")
Main.include("../../basins.jl/src/minimumassign/mxopt.jl")


class MixedDescentJulia(BaseMinimaFinder):
    def __init__(self, basefolder, paramstring):
        super().__init__(basefolder, paramstring)
        self.interfaced_potential = Main.PythonPotential(self.potential)

    def find_minimum_production(self, coords, opt_param_dict):
        """
        Find the minimum corresponding to the coords.

        Parameters
        ----------
        coords : np.ndarray
            Coordinates for which we want to find the corrsponding minima
        """
        nls = Main.NewtonLinesearch(
            self.interfaced_potential, coords, opt_param_dict['tol'])
        mxd = Main.Mixed_Descent(self.interfaced_potential, Main.CVODE_BDF(
        ), nls, coords, opt_param_dict['T'], opt_param_dict['rtol'], opt_param_dict['conv_tol'], opt_param_dict['tol'])
        try:
            mxd.run_b(mxd, opt_param_dict['n_steps'])
        except:
            return None
        return (mxd.optimizer.x0, mxd.converged, mxd.n_g_evals, mxd.iter_number, mxd.n_h_evals)