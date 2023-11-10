import numpy as np
import matplotlib.pyplot as plt

import FiniteDifferenceMethod as FDM
import FiniteElementMethod as FEM

def Run_all(Source, Dirichlet, Neumann, N):
    driver1 = FDM.FiniteDifferenceMethod(N)
    driver2 = FDM.FiniteDifferenceMethod(N)
    driver3 = FEM.MixedFiniteElementMethod(N)
    driver4 = FEM.ExoticFiniteElementMethod(N)
    data = {}
    results = {}
    driver1.solve(Source, Dirichlet, Neumann, True)
    driver2.solve(Source, Dirichlet, Neumann, False)
    driver3.solve(Source, Dirichlet, Neumann)
    driver4.solve(Source, Dirichlet, Neumann)
    data["fdm_fast"] = driver1.get_result()
    data["fdm_accurate"] = driver2.get_result()
    data["fem_mixed"] = driver3.get_result()
    data["fem_exotic"] = driver4.get_result()
    
    
    positions = np.array([subset[0][0] for subset in data["fdm_fast"]])
    for key in ["fdm_fast","fdm_accurate"]:
        results[key] = np.array([subset[1] for subset in data[key]])
    key = "fem_mixed"
    results[key] = np.array([subset[1][1] for subset in data[key]])
    key = "fem_exotic"
    results[key] = np.array([subset[1][0] for subset in data[key]])
    return data, positions, results




if __name__ == "main":
    print("Hello world from python wrapper function")