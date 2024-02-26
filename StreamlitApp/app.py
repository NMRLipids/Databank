# my_streamlit_app.py
import streamlit as st



import streamlit as st
import os
import sys
import numpy as np
import json
import MDAnalysis
import urllib.request
import yaml
import random
import matplotlib.pyplot as plt

import pandas as pd
import math




# Set page configuration to wide layout
st.set_page_config(layout="wide")


# This defines the path for the NMRlipids databank on your computer.
# Default is that this repository and the NMRlipids databank repository are cloned to the same folder.
# If this is not the case, change this to the folder where the NMRlipids databank repository is located.
databankPath = '../'

# Path for the NMRlipids databank
databankPath = '../'

# Access to functions defined in the NMRlipids databank
sys.path.insert(1, databankPath + '/Scripts/BuildDatabank/')
from databankLibrary import *

# Initialize the databank
systems = initialize_databank(databankPath)


# Sidebar title
st.sidebar.title("NMRlipids")


main_menu_selection = st.sidebar.radio("Main Menu", ["Ranking", "Plot Simulations"])




# Display secondary menu and content based on main menu selection
if main_menu_selection == "Ranking":
    ranking_selection = st.sidebar.selectbox("Ranking", ["Simulation Ranking", "Lipid Ranking", "POPC-CHOL Ranking"])

    # Content for Ranking Tables page
    if ranking_selection == "Simulation Ranking":
        st.title("Ranking tables of simulations against experimental data")


        # Selectbox for choosing sorting criteria
        Fragments = ['total', 'tails', 'headgroup', 'FormFactor']
        SortBasedOn = st.selectbox('Select Sorting Criteria', Fragments)

        # Button to show sorted table
        if st.button('Show Sorted Table'):
            st.write(f'Sorted based on {SortBasedOn} quality')

            # Path to the ranking file
            FFrankingPath = databankPath + '/Data/Ranking/SYSTEM_' + SortBasedOn + '_Ranking.json'
            
            # Load ranking data
            with open(FFrankingPath) as json_file:
                FFranking = json.load(json_file)

            # Display the table
            ShowTableGui(FFranking, 'TotalQuality')

    # Content for other pages...
    elif ranking_selection == "Lipid Ranking":
        st.title("Show ranking separately for each lipid")


        # Selectbox for choosing sorting criteria
        Fragments = ['total','sn-1','sn-2','headgroup']
        SortBasedOn = st.selectbox('Select Sorting Criteria', Fragments)

        # Button to show sorted table
        if st.button('Show Sorted Table'):
            st.write(f'Sorted based on {SortBasedOn} quality')

            # Path to the ranking file
            for lipid in lipids_dict:
                FFrankingPath = databankPath +  '/Data/Ranking/' + lipid + '_' + SortBasedOn + '_Ranking.json'
            
                # Load ranking data
                with open(FFrankingPath) as json_file:
                    FFranking = json.load(json_file)

                # Display the table
                ShowTableGui(FFranking, lipid)


    # Content for other pages...
    elif ranking_selection == "POPC-CHOL Ranking":
        st.title("Show sn-1 ranking only for systems with POPC and CHOLESTEROL")

            
        ### This is showing the ranking only for systems containing POPC and cholesterol
        
        FFrankingPath = databankPath + '/Data/Ranking/POPC_sn-1_Ranking.json'
        lipid = 'CHOL'
        
        
        with open(FFrankingPath) as json_file:
            FFranking = json.load(json_file)
        json_file.close()
        
        NewRank = []
        for i in FFranking:
            #print(i)
            #for tst in i:
            #    print(tst)
            if lipid in i['system']['COMPOSITION']:
                NewRank.append(i)
           
        ShowTableGui(NewRank,'POPC')


elif main_menu_selection == "Plot Simulations":

    plot_simulations_selection = st.sidebar.selectbox("Choose a Simulation Plot", ["ID"])


    # Content for new page...
    if plot_simulations_selection == "ID":
        st.title("Detailed Simulation Analysis")

        # Text input for ID
        input_id = st.text_input("Enter the Simulation ID", "")

        # Button to trigger analysis
        if st.button('Analyze Simulation'):
            if input_id:
                try:
                    # Convert input ID to integer
                    ID = int(input_id)
                    st.write(f'Analyzing simulation with ID: {ID}')

                    # Find the system with the given ID
                    selected_system = None
                    for system in systems:
                        if system['ID'] == ID:
                            selected_system = system
                            break

                    if selected_system:
                        # Display information and plots for the selected system
                        APL = CalcAreaPerMolecule(selected_system)
                        st.write('Membrane area per lipid:', APL)

                        thickness = GetThickness(selected_system)
                        st.write('Membrane thickness:', thickness)

                        # Display equilibration times and plots
                        st.write('Relative equilibration time for each lipid in the simulation:')
                        ShowEquilibrationTimesGui(selected_system)

                        st.write('Plot form factor and C-H bond order parameters from the simulation together with experimental data if available')
                        for lipid in selected_system['COMPOSITION']:
                            if lipid not in lipids_dict:
                                continue
                            plotSimulationGui(selected_system['ID'], lipid)  # Implement the plotSimulation function accordingly
                    else:
                        st.error("No simulation found with the given ID.")
                except ValueError:
                    st.error("Please enter a valid integer ID.")
            else:
                st.error("Please enter a Simulation ID.")





