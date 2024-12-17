""" File demonstrating how to start using DatabankLib API.
TODO: supplement!
"""
import DatabankLib.databankLibrary as dbl

if __name__ == "__main__":
    systems = dbl.initialize_databank()

    for system in systems:
        print(system)

        print("Fraction: ", dbl.calcLipidFraction(system))
