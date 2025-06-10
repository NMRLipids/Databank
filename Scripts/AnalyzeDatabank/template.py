""" File demonstrating how to start using DatabankLib API.
TODO: supplement!
"""
from DatabankLib.core import initialize_databank
import DatabankLib.databankLibrary as dbl

if __name__ == "__main__":
    systems = initialize_databank()

    for system in systems:
        print(system)
        print("Fraction of POPC: %.2f"
              % (dbl.calcLipidFraction(system, "POPC")*100))
        print("--")
