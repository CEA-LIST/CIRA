# ----------------------------------------------------------------------------
#                        LIST / DSCIN / LSTA
# ----------------------------------------------------------------------------
# 
# 
#  File        : CIRA
# 
#  Description : This program can be used to analyze the repairability of
#                a die-to-die interface. See the file README.md for
#                user instructions.
#  
#  Copyright (C) 2025 CEA-LIST
#  Author      : Théo BERMOND (theo.bermond@cea.fr)
#
#  Maintainer  : Adrian Evans (adrian.evans@cea.fr)
#
# Licensed under the LGPL-3.0 (the "License");
# You may not use this file except in compliance with the License.
# You may obtain a copy of the License at:
# https://www.gnu.org/licenses/lgpl-3.0.fr.html#license-text
#
# THE SOFTWARE IS PROVIDED “AS IS” AND THE AUTHOR DISCLAIMS ALL WARRANTIES WITH
# REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY AND
# FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY SPECIAL, DIRECT,
# INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER RESULTING FROM
# LOSS OF USE, DATA OR PROFITS, WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE OR
# OTHER TORTIOUS ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE OR
# PERFORMANCE OF THIS SOFTWARE.
# 
# ----------------------------------------------------------------------------

import argparse
import os
import pandas as pd
import yaml
import json
import drawsvg as dw
import numpy as np 
from math import sqrt
from collections import defaultdict, deque
from itertools import combinations
import random
import time
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter


#Parser initialisation.
parser = argparse.ArgumentParser() 

#Arguments for SVG displaying.
parser.add_argument('--BumpMap_file_name', type = str, help = 'Name of the Bumpmap file.', default = r'DEMO\DEMO_BumpMap_V3.yaml')
parser.add_argument('--Create_SVG', action = 'store_true', help = 'Create the bumpmap in SVG format.') 
parser.add_argument('--Aspect_file_name', type = str, help = 'Filename for the colors and shapes of bumps (must be a csv).', default = r'UserFiles\colors_shapes_dict.csv')
parser.add_argument('--Open_SVG', action = 'store_true', help = 'Directly open the SVG image that was generated using system default application.')
parser.add_argument('--Bump_Diameter', type = float, help = 'Bump size, in µm.', default = 1)
parser.add_argument('--Pitch', type = float, help = 'Define the pitch (or the minimal distance between two connections) of the interface. Will be used to scale the interface', default = 25)
parser.add_argument('--Input_X_scale', type = float, help = 'Scaling factor applying to X-axis in bump-map file.', default = 1)
parser.add_argument('--Input_Y_scale', type = float, help = 'Scaling factor applying to Y-axis in bump-map file.', default = 1)
parser.add_argument('--Legend', action = 'store_true', help = 'Flag to display the legend of the bumpmap.')
parser.add_argument('--Margin', type = int, help = 'Indicate the margin around the bumpmap, is a multiple of the Pitch.', default = 1)
parser.add_argument('--Bump_Name', action = 'store_true', help = 'Flag to display the bump names on the map.')
parser.add_argument('--Stroke_color', type = str, help = 'Define the stroke color.', default = 'black')
parser.add_argument('--Font', type = str, help = 'Choose the font.', default = 'Arial')
parser.add_argument('--Font_Size', type = float, help = 'Choose a scaling factor for the font.', default = 1)
parser.add_argument('--BumpMap_SVG_image_file_name', type = str, help = 'Filename for outputting the resulting image.', default = r'OutputFiles\BumpMap.svg' ) 
parser.add_argument('--Display_Reparability_SVG', action = 'store_true', help = 'Flag to display the reparability in the SVG image of the choosen interface.')

#Arguments for Reparability Stats.
parser.add_argument('--Reparability_Statistics', action = 'store_true', help = 'Flag to output the reparability statistics of the choosen interface.')
parser.add_argument('--Repair_Solutions', action = 'store_true', help = 'Flag to output the repair solution of every faults of the choosen interface.')
parser.add_argument('--IRL_file_name', type = str, help = 'Name of the IRL file that contains the repair informations of the chiplet.', default = r'DEMO\IRL_DEMO_V3.yaml')
parser.add_argument('--Fault_Table_file_name', type = str, help = 'The file that is written containing the repair informations for the interface.', default = r'OutputFiles\Fault_Table.yaml')
parser.add_argument('--Reparability_Table_file_name', type = str, help = 'The file that is written containing the repair informations for the interface.', default = r'OutputFiles\Repair_Table.yaml')
parser.add_argument('--Repair_Solutions_Table_file_name', type = str, help = 'The file that is written containing the repair informations for the interface.', default = r'OutputFiles\Repair_Solutions_Table.yaml')
parser.add_argument('--Print_Fault', action = 'store_true', help = 'Flag to print each fault.')

#Arguments for Fault Model.
parser.add_argument('--Fault_Type', type = str, help = 'Choose the fault type to analyze [Short, Open].', default = 'Short')
parser.add_argument('--Faults_Number', type = int, help = 'Choose the fault multiplicity for the given fault type. For example : Fault_Type = "Open" , Faults_Number = 2 corresponds to two arbitrary open faults anywhere on the interface.', default = 1)
parser.add_argument('--Short_Distance', type = float, help = 'Choose the upper threshold for the short distance in µm.', default = 26)
parser.add_argument('--Shorted_Bumps_Number', type = int, help = 'Choose the number of bumps affected by the fault, only works with short.For example : Shorted_Bumps_Number = 3 corresponds to three bumps shorted together within the Short_Distance.', default = 2)

#Arguments for MetaCIRA.
parser.add_argument('--Meta_Analysis', action = 'store_true', help = 'Flag to enable analysis using MetaCIRA.')
parser.add_argument('--System_description_file_name', type = str, help = 'The file that will contains the informations on the system.', default = r'MetaCIRA_Benchmark\System_1.yaml')
parser.add_argument('--System_Analysis', action = 'store_true', help = 'Flag to enable system analysis.')
parser.add_argument('--Min_Yield', type = float, help = 'Minimum electrical yield considered.', default = 0.95)
parser.add_argument('--Max_Yield', type = float, help = 'Maximum electrical yield considered.', default = 1)
parser.add_argument('--Number_of_faults_tested', type = int, help = 'Define the number of randomly generated faults for analyzing an interface.', default = 100)
parser.add_argument('--Number_of_electrical_yield_tested', type = int, help = 'Define the number of tested electrical yield.', default = 10)
parser.add_argument('--Log_Scale', action = 'store_true', help = 'Flag to set the electrical yield in log scale.') 

#Arguments for Bundle Repair Mechanisms (BRM).
parser.add_argument('--Bundle_Flag', action = 'store_true', help = 'A boolean to indicate if the interface and the repair mechanism is at the bundle level.')

#Arguments.
args = parser.parse_args()

BumpMap_file_name = args.BumpMap_file_name
Create_SVG = args.Create_SVG
Aspect_file_name = args.Aspect_file_name 
Open_SVG = args.Open_SVG
Bump_Diameter = args.Bump_Diameter
Pitch = args.Pitch
Input_X_scale = args.Input_X_scale
Input_Y_scale = args.Input_Y_scale
Legend = args.Legend
Margin = args.Margin
Bump_Name = args.Bump_Name
Stroke_Color = args.Stroke_color
Font = args.Font
Font_Size = args.Font_Size
BumpMap_SVG_image_file_name = args.BumpMap_SVG_image_file_name
Display_Reparability_SVG = args.Display_Reparability_SVG

Reparability_Statistics = args.Reparability_Statistics
Repair_Solutions = args.Repair_Solutions
Interface_IRL_file_name = args.IRL_file_name
Fault_Table_file_name = args.Fault_Table_file_name
Reparability_Table_file_name = args.Reparability_Table_file_name
Repair_Solutions_Table_file_name = args.Repair_Solutions_Table_file_name
Print_Fault = args.Print_Fault

Fault_Type = args.Fault_Type
Faults_Number = args.Faults_Number
Short_Distance = args.Short_Distance
Shorted_Bumps_Number = args.Shorted_Bumps_Number

Meta_Analysis = args.Meta_Analysis
System_description_file_name = args.System_description_file_name
System_Analysis = args.System_Analysis
Min_Yield = args.Min_Yield
Max_Yield = args.Max_Yield
Number_of_faults_tested = args.Number_of_faults_tested
Number_of_electrical_yield_tested = args.Number_of_electrical_yield_tested
Log_Scale = args.Log_Scale

Bundle_Flag = args.Bundle_Flag

start = time.time()

# Section 1 : Data Loading and Preparation.
def file_loading_as_a_DataFrame(file_name):
    """
    Loads a specified file into a pandas DataFrame based on its extension.

    Parameters:
    - file_name (str): Path to the file to be loaded.

    Returns:
    - pd.DataFrame: A DataFrame containing the loaded data if successful, otherwise None.
    """

    # Extract the file name and extension
    name, ext = os.path.splitext(file_name)
    
    # Determine the file type and load the data accordingly
    if ext[1:] == 'csv':  # Check for CSV files
        # Load CSV file into a DataFrame
        df = pd.read_csv(file_name)  # Read CSV file directly
        
    elif ext[1:] in ['yaml', 'yml']:  # Check for YAML or YML files
        # Load YAML or YML file into a DataFrame
        with open(file_name, 'r') as file:
            data = yaml.safe_load(file)
            df = pd.DataFrame(data)  # Create a DataFrame from the loaded data
        
    elif ext[1:] == 'json':  # Check for JSON files
        # Load JSON file into a DataFrame
        with open(file_name, 'r') as file:
            data = json.load(file)
            df = pd.DataFrame(data)  # Create a DataFrame from the loaded data
        
    elif ext[1:] == 'xml':  # Check for XML files
        # Load XML file into a DataFrame
        with open(file_name, 'r') as file:
            df = pd.read_xml(file)  # Read XML file directly
    
    else:  # Handle unsupported file formats
        # Print an error message for unsupported file formats
        print("Format de fichier non supporté")  
        return None  # Return None to indicate failure

    # Return the DataFrame containing the loaded data
    return df

def Repair_yaml_file_loading_into_a_dataframe(Interface_IRL_file_name):
    """
    This function loads a YAML file containing repair information into a pandas DataFrame.
    The YAML file is expected to have a specific structure with keys representing repair chains,
    and nested dictionaries containing information about functional and physical ports.

    Parameters:
    - Interface_IRL_file_name (str): Path to the YAML file containing repair information.

    Returns:
    - pd.DataFrame: A DataFrame containing the repair capabilities extracted from the YAML file.
    """
    # This function loads a YAML file containing repair information into a pandas DataFrame.
    # The YAML file is expected to have a specific structure with keys representing repair chains,
    # and nested dictionaries containing information about functional and physical ports.

    # Open the YAML file and load its contents into a dictionary.
    with open(Interface_IRL_file_name, 'r') as file:
        data = yaml.safe_load(file)
        
    # Initialize an empty DataFrame with the desired columns.
    df_Repair_Capabilities = pd.DataFrame(columns=['Signal', 'Connection', 'Mux', 'Sel', 'Status', 'RepairChain'])

    # Extract information from the dictionary.
    # The keys of the dictionary represent repair chains.
    RepairChains = list(data.keys())
    for RepairChain in RepairChains:
    # Each repair chain contains a dictionary with functional ports as keys.
        RepairChain_dict = data[RepairChain]
        for functional_port, functional_port_info in RepairChain_dict.items():
        # Extract the name of the functional port.
            functional_port_name = functional_port_info['Name']
            for physical_port, physical_port_info in functional_port_info.items():
            # Skip the 'Name' key as it is already processed.
                if physical_port != 'Name':
                # Extract the name of the physical port.
                    physical_port_name = physical_port_info.get('To')
                # Extract multiplexer (Mux) and selector (Sel) information.
                    mux_info = physical_port_info.get('Control', {}).get('Mux', '')
                    sel_info = physical_port_info.get('Control', {}).get('Sel', '')

                # Create a new row with the extracted information.
                    new_row = {
                        'Signal': [functional_port_name],
                        'Connection': [physical_port_name],
                        'Mux': [mux_info],
                        'Sel': [sel_info],
                        'Status': [physical_port],
                        'RepairChain': [RepairChain],
                    }

                # Append the new row to the DataFrame.
                    df_Repair_Capabilities = pd.concat([df_Repair_Capabilities, pd.DataFrame(new_row)], ignore_index=True)

    # Return the DataFrame containing the repair capabilities.
        return df_Repair_Capabilities

def Repair_IRL_file_loading_into_a_dataframe(Interface_IRL_file_name):
    """
    Standalone function that reads an IRL file and returns the same DataFrame
    as the original YAML function.
    Args:
        Interface_IRL_file_name (str): Path to the .irl file.
        
    Returns:
        pd.DataFrame: DataFrame with the columns ['Signal', 'Connection', 'Mux', 'Sel', 'Status', 'RepairChain'].
    """
    try:
        # Read the IRL file
        with open(Interface_IRL_file_name, 'r', encoding='utf-8') as file:
            content = file.read()
        
        # Remove header comments to get pure YAML
        lines = content.split('\n')
        content_lines = [line for line in lines if not line.strip().startswith('#')]
        clean_content = '\n'.join(content_lines).strip()
        
        # Parse as YAML (IRL format is compatible)
        data = yaml.safe_load(clean_content)
        
        # Initialize an empty DataFrame with the desired columns.
        df_Repair_Capabilities = pd.DataFrame(columns=['Signal', 'Connection', 'Mux', 'Sel', 'Status', 'RepairChain'])

        # Extract information from the dictionary (same logic as the YAML function)
        RepairChains = list(data.keys())
        for RepairChain in RepairChains:
            RepairChain_dict = data[RepairChain]
            for functional_port, functional_port_info in RepairChain_dict.items():
                functional_port_name = functional_port_info['Name']
                for physical_port, physical_port_info in functional_port_info.items():
                    if physical_port != 'Name':
                        physical_port_name = physical_port_info.get('To')
                        mux_info = physical_port_info.get('Control', {}).get('Mux', '')
                        sel_info = physical_port_info.get('Control', {}).get('Sel', '')

                        new_row = {
                            'Signal': [functional_port_name],
                            'Connection': [physical_port_name],
                            'Mux': [mux_info],
                            'Sel': [sel_info],
                            'Status': [physical_port],
                            'RepairChain': [RepairChain],
                        }

                        # Append the new row to the DataFrame
                        df_Repair_Capabilities = pd.concat([df_Repair_Capabilities, pd.DataFrame(new_row)], ignore_index=True)

        return df_Repair_Capabilities
        
    except Exception as e:
        print(f"Error reading the IRL file: {e}")
        return pd.DataFrame(columns=['Signal', 'Connection', 'Mux', 'Sel', 'Status', 'RepairChain'])

def Avoid_bump_name_iteration(file_name):
    """
    This function processes a DataFrame loaded from a specified file to ensure that all bump names are unique.
    It iterates through the 'Name' column of the DataFrame, appending a suffix to any duplicate names to make them unique.

    Parameters:
    - file_name (str): Path to the file to be loaded. The file should contain a 'Name' column.

    Returns:
    - pd.DataFrame: A DataFrame with unique bump names.
    """
    # Load the DataFrame from the specified file
    df = file_loading_as_a_DataFrame(file_name)

    # Dictionary to keep track of name counts
    name_counts = defaultdict(int)

    # List to store the modified names
    modified_names = []

    # Iterate through each name in the 'Name' column
    for name in df['Name']:
        # Get the current count of the name
        count = name_counts[name]
        # Increment the count for the name
        name_counts[name] += 1
        
        # If the name has not been seen before, append it as is
        if count == 0:
            modified_names.append(name)
        else:
            # If the name has been seen before, append it with a suffix indicating the count
            modified_names.append(f"{name}_{count}")

    # Update the 'Name' column with the modified names
    df['Name'] = modified_names
    
    # Return the DataFrame with unique bump names
    return df

# Section 2 : SVG Generation.    
def Display_SVG(BumpMap_file_name, Aspect_file_name, BumpMap_SVG_image_file_name, Open_SVG,
                Bump_Diameter, Pitch, Input_X_scale, Input_Y_scale,
                Legend, Margin, Bump_Name, Stroke_Color,
                Font, Font_Size, Display_Reparability_SVG):
    """
    This function generates an SVG image of the bump map based on the provided parameters.
    It reads the bump map data from a specified file, applies scaling factors, and creates
    an SVG representation of the bumps with optional features such as legend, bump names,
    and reparability display.

    Parameters:
    - BumpMap_file_name (str): Path to the bump map file.
    - Aspect_file_name (str): Path to the file containing colors and shapes of bumps.
    - BumpMap_SVG_image_file_name (str): Path to save the generated SVG image.
    - Open_SVG (bool): Flag to open the generated SVG image using the system default application.
    - Bump_Diameter (float): Diameter of the bumps in micrometers.
    - Pitch (float): Pitch of the interface, defining the minimal distance between two connections.
    - Input_X_scale (float): Scaling factor for the X-axis in the bump map file.
    - Input_Y_scale (float): Scaling factor for the Y-axis in the bump map file.
    - Legend (bool): Flag to display the legend of the bump map.
    - Margin (int): Margin around the bump map, specified as a multiple of the pitch.
    - Bump_Name (bool): Flag to display the bump names on the map.
    - Stroke_Color (str): Color of the stroke around the bumps.
    - Font (str): Font to be used for text in the SVG image.
    - Font_Size (float): Scaling factor for the font size.
    - Display_Reparability_SVG (bool): Flag to display reparability information in the SVG image.
    """
    # Code to generate the SVG image
    # (The code remains the same as provided in the original script)

    # Load bump map data into a DataFrame
    df = Avoid_bump_name_iteration(BumpMap_file_name)

    # Load aspect data into a DataFrame
    aspect = file_loading_as_a_DataFrame(Aspect_file_name)

    # Scale the x and y coordinates of the bumps
    df['X'] = df['X'] * Input_X_scale 
    df['Y'] = df['Y'] * Input_Y_scale 

    # Determine the maximum and minimum distances in the x and y directions
    max_X = max(df['X']) 
    min_X = min(df['X'])
    max_Y = max(df['Y'])
    min_Y = min(df['Y']) 

    n_row = len(df['Y'].unique())
    n_columns = len(df['X'].unique())
    
    X_Pitch = (max_X-min_X)/(n_columns/2)
    Y_Pitch = (max_Y-min_Y)/(n_row/2)

    # If the pitch is not defined by the user, it will be extracted by taking the mean between X_Pitch and Y_Pitch.  
    if Pitch == 0:
        Pitch = (X_Pitch+Y_Pitch)/2
        print('Warning : The pitch is not defined by the user. To define it, please use : --Pitch int (in µm)')

    # Define the margin that scales with the pitch
    Margin = Margin * 0.7 * Pitch 

    # Define the canvas dimensions with a margin
    width = ((max_X - min_X) + 2 * Margin) 
    height = ((max_Y - min_Y) + 2 * Margin) 
    
    # Scale the bump size (s) and the margin with the bump diameter
    s = 0.2 * Bump_Diameter * Pitch
    Margin = Margin * Bump_Diameter 

    # Define a list of bump types for the legend if the legend is enabled
    if Legend:
        legend_list = df['Type'].unique().tolist()
        if 'SPARE' not in legend_list:
            legend_list.append('SPARE')
        legend_margin = 2.5 * Margin
    else:
        legend_list = []
        legend_margin = 0

    # Define functions to create different shapes (circle, triangle, square)
    def Circle(x, y, s, color, Stroke_Color, a):
        return dw.Circle(0, 0, s, fill=color, stroke=Stroke_Color, stroke_width = s / 10, transform=f'translate({x},{y})', fill_opacity = a, stroke_opacity = a)

    def Triangle(x, y, s, color, Stroke_Color, a):
        return dw.Lines(-s, s, s, s, 0, -s, fill=color, close='true', stroke=Stroke_Color, stroke_width = s / 10, transform = f'translate({x},{y})',  fill_opacity = a, stroke_opacity = a)

    def Square(x, y, s, color, Stroke_Color, a):
        return dw.Rectangle(-s, -s, 1.7 * s, 1.7 * s, fill=color, stroke=Stroke_Color, stroke_width = s / 10, transform = f'translate({x},{y})',  fill_opacity = a, stroke_opacity = a)

    def Text(x, y, bump_type, s, Font_Size, Font, BaseLine, anchor, angle):
        return dw.Text(f'{bump_type}', font_size=Font_Size * 0.8 * s, x = x, y = y, dominant_baseline = BaseLine, font_family = Font, text_anchor = anchor, transform = f'rotate ({angle}, {x}, {y})')

    def Line(X1, Y1, X2, Y2, color, s):
        return dw.Line(X1, Y1, X2, Y2, stroke = color, stroke_width = s / 6)
    
    # Dictionary mapping bump shapes to their corresponding functions
    shape_function_dict = {
        'Circle': Circle,
        'Triangle': Triangle,
        'Square': Square,
    }
    
    # Define the SVG canvas, including a margin for the legend if enabled
    bumpmap = dw.Drawing(width + legend_margin + 2 * Margin, height + 2 * Margin, origin=(-2.5 * Margin + min_X, -2.5 * Margin + min_Y))

    # If Bump_Name is called then the alpha parameters (opacity) will be set at 0.7. 
    if Bump_Name:
        a = 0.7
    else:
        a = 1

    # If the Display_Reparability_SVG flag is set, generate the repair solutions table
    if Display_Reparability_SVG:   
        # Call the function to generate repair statistics using a logic solver
        Repair_Solutions_Table = Repair_Statistics_using_LogicSolver(BumpMap_file_name, Fault_Type, Shorted_Bumps_Number, Short_Distance, Faults_Number, Interface_IRL_file_name, Reparability_Table_file_name, Fault_Table_file_name, Print_Fault) 

        if Fault_Type == 'Short' and Shorted_Bumps_Number == 2:
            
            for index, fault in Repair_Solutions_Table.iterrows():
                # Initialize a list to store the coordinates of the bumps involved in the fault
                fault_coordinates = []
                # Extract the list of bumps involved in the fault
                bumps = fault['Fault']
                # Iterate through each bump in the fault
                for bump in bumps:
                    # Initialize a list to store the coordinates of the current bump
                    bump_coordinates = []
                    # Get the x-coordinate of the bump from the DataFrame
                    X_coord = df[df['Name'] == bump]['X'].values[0]
                    # Get the y-coordinate of the bump from the DataFrame
                    Y_coord = df[df['Name'] == bump]['Y'].values[0]
                    # Append the coordinates to the bump_coordinates list
                    bump_coordinates.append(X_coord)
                    bump_coordinates.append(Y_coord)
                    # Append the bump_coordinates list to the fault_coordinates list
                    fault_coordinates.append(bump_coordinates)
                
                # Get the repair type for the current fault
                Reparability = fault['Repair_Type']
                # Get the color corresponding to the repair type from the aspect DataFrame
                color = aspect.loc[aspect['Type'] == Reparability, 'Color'].values[0]

                # If the repair type is 'Catastrophic', draw a thick, solid line between the bumps
                if Reparability == 'Catastrophic':
                    bumpmap.append(Line(fault_coordinates[0][0], fault_coordinates[0][1], fault_coordinates[1][0], fault_coordinates[1][1], color, 2*s))
                
                # If the repair type is 'Benign', draw a thin, dashed line between the bumps
                if Reparability == 'Benign':
                    bumpmap.append(dw.Line(fault_coordinates[0][0], fault_coordinates[0][1], fault_coordinates[1][0], fault_coordinates[1][1], stroke = color, stroke_width = s/6, stroke_dasharray = '2,2'))

                # For other repair types, draw a regular line between the bumps
                else: 
                    bumpmap.append(Line(fault_coordinates[0][0], fault_coordinates[0][1], fault_coordinates[1][0], fault_coordinates[1][1], color, s))

            # If the legend is enabled, add the repair types to the legend list
            if Legend:    
                for i in Repair_Solutions_Table['Repair_Type'].unique().tolist():  
                    legend_list.append(i)
        else: 
            print('Warning : Display_Reparability_SVG does not work with others fault model than 2-bumps short. Please specify the correct fault model with : --Fault_Type and --Shorted_Bumps_Number')
               
    # Iterate over each bump in the DataFrame
    for i in range(df.shape[0]):

        # Define the x, y coordinates using the columns from the DataFrame
        x = df['X'][i]
        y = df['Y'][i]

        bump_type = df['Type'][i]

        # Recover colors and shapes from the aspect DataFrame corresponding to the bump type
        color = aspect.loc[aspect['Type'] == bump_type, 'Color'].values[0]
        shape = aspect.loc[aspect['Type'] == bump_type, 'Shape'].values[0]

        # If a bump is a spare, get the shape from the aspect file
        if df['Spare'][i] == True:
            shape = aspect.loc[aspect['Type'] == 'SPARE', 'Shape'].values[0]

        # Generate a shape for each bump at the position x, y with the shape, color, and stroke color required
        bumpmap.append(shape_function_dict[shape](x, y, s, 'white', Stroke_Color, 1))
        bumpmap.append(shape_function_dict[shape](x, y, s, color, Stroke_Color, a))

        # Show the bump name if required and if the bump is not a POWER or GND connection  
        if Bump_Name:
            if df['Type'][i] != 'GND' and df['Type'][i] != 'POWER': 
                bump_name = df['Name'][i].replace('_phy', '')
                bumpmap.append(Text(x, y, bump_name, 1.2*s, Font_Size, Font, 'middle', 'start', -15))

    # Show the legend if required
    if Legend:
        
        X_edge = width + legend_margin - Margin 
        for j in legend_list:
            index = legend_list.index(j)
            color = aspect.loc[aspect['Type'] == j, 'Color'].values[0]
            shape = aspect.loc[aspect['Type'] == j, 'Shape'].values[0]

            X_Shape = X_edge - 0.65 * legend_margin
            Y_Shape = (index * Pitch)

            if shape == 'Line':
                if j == 'Catastrophic':
                    bumpmap.append(Line(X_Shape - s, (Y_Shape + s) + min_Y, X_Shape + s, (Y_Shape - s) + min_Y, color, 2*s))
                elif j == 'Benign':
                    bumpmap.append(dw.Line(X_Shape - s, (Y_Shape + s) + min_Y, X_Shape + s, (Y_Shape - s) + min_Y, stroke = color, stroke_width = s/6, stroke_dasharray = '5.5'))   
                else:
                    bumpmap.append(Line(X_Shape - s, (Y_Shape + s) + min_Y, X_Shape + s, (Y_Shape - s) + min_Y, color, s))
            else:
                bumpmap.append(shape_function_dict[shape](X_Shape, Y_Shape + min_Y, s, color, Stroke_Color, 1))

            bumpmap.append(Text(X_edge - 0.53 * legend_margin, Y_Shape + min_Y, j, 1*s, Font_Size, Font, 'middle', 'start', 0))

        bumpmap.append(dw.Rectangle(X_edge - 0.8 * legend_margin, (-0.5 * Margin) + min_Y, legend_margin, len(legend_list) * Pitch, fill = 'none', stroke = 'black', stroke_width = s/6))

    n_step_Y = n_row 
    n_step_X =  n_columns 

    # Y Axis
    # Define the starting x-coordinate for the Y-axis line
    Yx = -0.75 * Margin + min_X
    # Define the starting y-coordinate for the Y-axis line
    Y1 = - Margin + min_Y
    # Define the ending y-coordinate for the Y-axis line
    Y2 = 0.5 * Margin + max_Y
    # Draw the Y-axis line
    bumpmap.append(Line(Yx, Y1, Yx, Y2, 'black', s))

    # Iterate over the range of y-coordinates to draw tick marks and labels
    for y in np.linspace(0, max_Y - min_Y, n_step_Y):
        # Round the y-coordinate to 2 decimal places
        y = round(y, 2)
        # Draw the tick mark label for the y-coordinate
        bumpmap.append(Text(Yx - 0.2 * Margin, y + min_Y, f'{y} µm', 1.2*s, Font_Size, Font, 'middle', 'end', 0))
        # Draw the tick mark line for the y-coordinate
        bumpmap.append(Line(Yx - 0.1 * Margin, y + min_Y, Yx + 0.1 * Margin, y + min_Y, 'black', s))

    # X Axis
    # Define the starting y-coordinate for the X-axis line
    Xy = -0.75 * Margin + min_Y
    # Define the starting x-coordinate for the X-axis line
    X1 = - Margin + min_X
    # Define the ending x-coordinate for the X-axis line
    X2 = 0.5 * Margin + max_X
    # Draw the X-axis line
    bumpmap.append(Line(X1, Xy, X2, Xy, 'black', s))

    # Iterate over the range of x-coordinates to draw tick marks and labels
    for x in np.linspace(0, max_X - min_X, n_step_X):
        # Round the x-coordinate to 2 decimal places
        x = round(x, 2)
        # Draw the tick mark label for the x-coordinate
        bumpmap.append(Text(x + min_X, Xy - 0.2 * Margin, f'{x} µm', 1.2*s, Font_Size, Font, 'middle', 'end', 35))
        # Draw the tick mark line for the x-coordinate
        bumpmap.append(Line(x + min_X, Xy - 0.1 * Margin, x + min_X, Xy + 0.1 * Margin, 'black', s))


    # Save the SVG image and open it directly if required
    bumpmap.save_svg(BumpMap_SVG_image_file_name)
    if Open_SVG:
        os.system(f'inkscape {BumpMap_SVG_image_file_name}') 
 
# Section 3 : Raw Solvers
def LogicSolver(Chain_list, Route_Table, df_bump, fault):
    """
    This function determines the reparability of a fault based on the available repair chains and route table.
    The function will evaluate how many spares are available per affected repair chain. 
    It will compare this number to the number of connection that need to be repaired.
    The fault can be repaired if there is enough spares in the affected repair chains.

    Parameters:
    - Chain_list (list): List of repair chains.
    - Route_Table (pd.DataFrame): DataFrame containing the route table with repair information.
    - df_bump (pd.DataFrame): DataFrame containing bump information.
    - fault (list): List of faulty connections.

    Returns:
    - str: 'Repairable' if the fault can be repaired, 'Unrepairable' otherwise.
    """

    # Initialize the unrepairable flag
    UnrepairableFlag = False

    # Iterate over each repair chain of the Chain_list
    for Chain in Chain_list:
        # Initialize lists to keep track of spare connections and faulty connections
        Spare_list = []
        Spare_Count = 0

        # If the fault is not marked as Unrepairable yet.
        if UnrepairableFlag == False: 
            Faulty_Connections = []

            # Get the signal set (a list of the signal in the affected repair chain) for the current repair chain
            SignalSet = Route_Table[(Route_Table['RepairChain'] == Chain)]
      
            # Iterate over each signal in the signal set
            for index, signal_row in SignalSet.iterrows():

                # Check if the connection is a spare or if it only exists for repair purpose in the Route_Table, the second condition is related to the special case of HBM2 mode 1, use of DBI connections for repair. 
                if (df_bump.loc[df_bump['Name'] == signal_row['Connection'], 'Spare'].values[0] == True) or Route_Table[(Route_Table['Connection'] == signal_row['Connection']) & (Route_Table['Status'] == 'Default')].empty:
                    if signal_row['Connection'] not in Spare_list:
                        Spare_list.append(signal_row['Connection']) 
                        # Increment the spare count per 1
                        Spare_Count += 1

            # If the connection is found in the SignalSet's 'Connection' values, it is added to the Faulty_Connections list
            for connection in fault:

                if connection in SignalSet['Connection'].values: # HBM2 Mode 1 DBI specifics
                    Faulty_Connections.append(connection)

                signal = connection.replace('_phy', '')
                if signal in SignalSet['Signal'].values:
                    # If the connection does not have any repair option. In other words, if it only has a default connection, then the fault is marked as Unrepairable
                    if Route_Table[(Route_Table['Signal'] == signal) & (Route_Table['Status'] != 'Default')].empty:
                        UnrepairableFlag = True
                        
            # If the number of faulty connections exceeds the spare count, mark as Unrepairable
            if len(Faulty_Connections) > Spare_Count:
                UnrepairableFlag = True

    # Return the reparability status
    if UnrepairableFlag == True:
        return 'Unrepairable'
    else:
        return 'Repairable'

def RecursiveSolver(Used_IS_list, Routed_PFS_list, PFS_to_route_list, PFS_counter, Solution, Solutions_dict_per_RepairChain, Route_Table, Repair_Type):
    """
    This function recursively solves the routing problem for a given set of Physical Functionnal Sources (PFS).
    It attempts to route each PFS to an available Interconnect Source (IS) and keeps track of the solutions found.
        
    Parameters:
    - Used_IS_list (list): List of Interconnect Sources (IS) that have already been used.
    - Routed_PFS_list (list): List of Physical Functional Sources (PFS) that have already been routed.
    - PFS_to_route_list (list): List of PFS that need to be routed.
    - PFS_counter (int): Counter to keep track of the current PFS being routed.
    - Solution (pd.DataFrame): Current solution being constructed.
    - Solutions_dict_per_RepairChain (dict): Dictionary to store all possible solutions for the repair chain.
    - Route_Table (pd.DataFrame): DataFrame containing the route table with repair information.

    Returns:
    - None: The function modifies the Solutions_dict_per_RepairChain dictionary in place.
    """
    
    # Base case: If all PFS are routed, store the solution
    if Routed_PFS_list == PFS_to_route_list: #If we have routed every PFS 
        Solutions_dict_per_RepairChain[f'Solution_{len(Solutions_dict_per_RepairChain)}'] = Solution.set_index('Mux')['Sel'].to_dict()
        return

    # if there exists no solution for the fault yet
    if len(Solutions_dict_per_RepairChain) == 0:

        # Get the current PFS to route
        PFS_to_route = PFS_to_route_list[PFS_counter]
        
        # Find possible routes for the current PFS in Route_Table
        possible_PFS_route = Route_Table[(Route_Table['Signal'] == PFS_to_route)]

        for index, route in possible_PFS_route.iterrows():
            # We marked the connection corresponding to the current route as used
            Used_IS = route['Connection']

            # Check if the IS and PFS are not already used or routed
            if Used_IS not in Used_IS_list and PFS_to_route not in Routed_PFS_list:

                # Update the lists and the current solution
                Used_IS_list.append(Used_IS)
                Routed_PFS_list.append(PFS_to_route)
                Solution = Solution._append(route, ignore_index=True)

                # Recursively solve for the next PFS
                RecursiveSolver(Used_IS_list, Routed_PFS_list, PFS_to_route_list, PFS_counter + 1, Solution, Solutions_dict_per_RepairChain, Route_Table, Repair_Type)

                if len(Solutions_dict_per_RepairChain) == 0:
                    # Backtrack to explore other routes
                    Routed_PFS_list.pop()
                    Used_IS_list.pop()
                    Solution = Solution.drop(Solution.index[-1])

def BundleSolver(df_bump, fault, Route_Table):
    """
    This function determines the reparability of a fault based on the available bundles and route table.
    It only works with interface made of bundles. 
    See the README for more informations.

    Parameters:
    - df_bump (pd.DataFrame): DataFrame containing bump information.
    - fault (list): List of faulty connections.
    - Route_Table (pd.DataFrame): DataFrame containing the route table with repair information.

    Returns:
    - str: 'Repairable' if the fault can be repaired, 'Unrepairable' otherwise.
    """
    # Initialize an empty list to store unique bundles associated with the faulty connections
    Bundle_list = []

    # Initialize a flag to determine if the fault is unrepairable
    UnrepairableFlag = False

    # Iterate through each connection in the fault list
    for connection in fault:
        # Get the bundle associated with the current connection
        Bundle = df_bump.loc[df_bump['Name'] == connection, 'Bundle'].values[0]
        
        # If the bundle is not already in the Bundle_list, add it
        if Bundle not in Bundle_list:
            Bundle_list.append(Bundle)

    # Iterate through each unique bundle in the Bundle_list
    for Bundle in Bundle_list:
        # If the fault is not already marked as unrepairable and the bundle is not None
        if UnrepairableFlag == False and Bundle != None:
            # Check if there is a default route for the bundle in the Route_Table, in other words : check if the bundle is for repair purpose only 
            if not Route_Table[(Route_Table['Connection'] == Bundle) & (Route_Table['Status'] == 'Default')].empty:
                # Get the signal name associated with the bundle
                Bundle_Signal = Bundle.replace('_phy', '')

                # Get the repair connection for the bundle signal
                Repair_Bundle = Route_Table[(Route_Table['Signal'] == Bundle_Signal) & (Route_Table['Status'] == 'Repair')]['Connection'].values[0]

                # If the repair bundle assigned to the current bundle is in the Bundle_list, mark the fault as unrepairable
                if Repair_Bundle in Bundle_list:
                    UnrepairableFlag = True

    # Return the reparability status based on the UnrepairableFlag
    if UnrepairableFlag == True: 
        return 'Unrepairable'
    else:
        return 'Repairable'
    
# Section 4 : Reparability Statistics
def Fault_Table_Generator(Interface_IRL_file_name, BumpMap_file_name, Faults_Number, Shorted_Bumps_Number, Short_Distance, Fault_Type, Fault_Table_file_name):
     

    # Get the IRL and Bumpmap file 
    Route_Table = Repair_IRL_file_loading_into_a_dataframe(Interface_IRL_file_name)
    df_bump = Avoid_bump_name_iteration(BumpMap_file_name)

    # Define the type of each connection regarding the reparability
    fault_reparation = {
        'POWER': 'Benign',
        'GND': 'Benign',
        'DATA': 'Repair', 
        'CLK': 'Repair', 
        'ADDR': 'Repair', 
        'SIDEBAND': 'Repair', 
        'SPARE': 'Benign', 
        'NONE': 'Benign'}
    
    # Initiate the fault table as an empty dataframe
    Fault_Table = pd.DataFrame()

    def euclidean_distance(point1, point2):

        """
        Calculate the Euclidean distance between two bumps.
        
        Args:
            point1: First point with coordinates accessible via dictionary-like access (e.g., point1['X'], point1['Y'])
            point2: Second point with coordinates accessible via dictionary-like access
        
        Returns:
            float: Euclidean distance between the bumps
        """
        # Assumes both bumps have the same dimensions
        # Extract all columns that might represent coordinates (e.g., X, Y, Z...)

        coords = [col for col in point1.index if col in point2.index and col.upper() in ['X', 'Y', 'Z']]
        
        # Calculate sum of squared differences for each coordinate
        sum_squared_diff = sum((point1[coord] - point2[coord])**2 for coord in coords)
        
        # Return the square root of the sum
        return np.sqrt(sum_squared_diff)

    def is_short(bumps, threshold):
        """
        Determine if a set of bumps forms a short.
        A short is defined as:
        1. Each point is connected to at least one other point (distance < threshold)
        2. All bumps form a single connected component (can reach any point from any other)
        
        Args:
            bumps: List of bumps, where each point is a pandas Series or similar with X, Y, etc. coordinates
            threshold: Maximum distance for two bumps to be considered connected
        
        Returns:
            bool: True if bumps form a short, False otherwise
        """

        if len(bumps) == 0:
            # Empty list can't be short
            return False
        
        elif len(bumps) == 1:
            # A single point is an Open
            return True
        
        else:
            # Build an adjacency list representing the graph
            # Each point is a node, and edges connect bumps with distance < threshold
            graph = defaultdict(list)
            
            # Check if each point has at least one connection
            has_connection = [False] * len(bumps)
            
            # Build the graph
            for i in range(len(bumps)):
                for j in range(i+1, len(bumps)):  # Only check each pair once
                    if euclidean_distance(bumps[i], bumps[j]) < threshold:
                        graph[i].append(j)
                        graph[j].append(i)
                        has_connection[i] = True
                        has_connection[j] = True
            
            # Check if any point has no connections
            if not all(has_connection):
                return False
            
            # Use BFS to check if the graph is connected
            # Start from the first point (index 0)
            visited = [False] * len(bumps)
            queue = deque([0])
            visited[0] = True
            
            while queue:
                node = queue.popleft()
                for neighbor in graph[node]:
                    if not visited[neighbor]:
                        visited[neighbor] = True
                        queue.append(neighbor)
            
            # If all nodes are visited, the graph is connected
            return all(visited)

    # Convert the df_bump as a list of Series. Each containing the informations of a connections=
    all_bumps = [df_bump.iloc[i] for i in range(len(df_bump))]

    # Check if n_bumps is valid
    if Shorted_Bumps_Number < 1 or Shorted_Bumps_Number > len(all_bumps):
        raise ValueError(f"Shorted_Bumps_Number must be between 1 and {len(all_bumps)}")
    
    # Set the number of bumps that will be affected by the fault.  
    if Fault_Type == 'Open':
        Bumps_Number = 1

    if Fault_Type == 'Short':
        Bumps_Number = Shorted_Bumps_Number

    # Test each combination
    # The first use of combinations is for the short faults. For example, for the 3-bump short, we will check every combination of three connections (Bumps_Number)
    # The second use of combinations is for the multiple fault scenario. For example, for the two 2-bump short, we will check every combination of two combinations (Faults_Number) of two connections
    for combo_of_combo in combinations(combinations(range(len(all_bumps)), Bumps_Number), Faults_Number):  

        index_list = []
        for combo_index in combo_of_combo:
         
            # FIXME : Does not work for scenarios with multiple short faults. The code puts the indices of each connection in the same list.
            # A double short with two connections will be treated as a short with 4 connections.
            # Furthermore, in two different short (happening at the same time): the same connection may appear twice.
            for index in combo_index:
                index_list.append(index)
            
        combo_bumps = [all_bumps[i] for i in index_list]
     
        # Initialize a flag to determine if the fault is valid
        FaultFlag = True
        # If the fault type is 'Short', check if the bumps form a valid short
        if Fault_Type == 'Short':
            FaultFlag = is_short(combo_bumps, Short_Distance)

        # If the FaultFlag is True, proceed with the fault analysis
        if FaultFlag == True:
            # Initialize lists and variables to store fault information
            new_row = []
            Chain_list = []
            Route_Table_copy = Route_Table.copy()
            fault = []
            Repair_Type = 'Benign'
            GNDFlag = False
            POWERFlag = False

            # Iterate through each bump in the combination of bumps
            for bump in combo_bumps:  
                # Add the bump name to the fault list
                fault.append(bump['Name'])  

                # Check if the fault is catastrophic
                if bump['Type'] == 'GND':
                    GNDFlag = True
                if bump['Type'] == 'POWER':
                    POWERFlag = True
                if POWERFlag == True and GNDFlag == True and Fault_Type == 'Short':
                    Repair_Type = 'Catastrophic'

                # Generate a list containing the involved repair chains 
                if not Route_Table[(Route_Table['Connection'] == bump['Name'])]['RepairChain'].empty:
                    RepairChain = Route_Table[(Route_Table['Connection'] == bump['Name'])]['RepairChain'].values[0]
                    Chain_list.append(RepairChain)

                # If the bump is a spare, update its type to 'SPARE'
                if bump['Spare'] == True:
                    bump = bump.copy()
                    bump['Type'] = 'SPARE'

                # If the bump does not have a default route and is not of type 'DATA', update its type to 'NONE' # HBM2 Mode 1 DBI specifics
                if Route_Table_copy[(Route_Table_copy['Connection'] == bump['Name']) & (Route_Table_copy['Status'] == 'Default')].empty and bump['Type'] != 'DATA':
                    bump = bump.copy()
                    bump['Type'] = 'NONE'

                # Check if the fault needs a repair action
                if fault_reparation[bump['Type']] == 'Repair' and Fault_Type != 'Catastrophic':
                    Repair_Type = 'Repair'

            # Add the fault information to the new row
            new_row.insert(0, set(Chain_list))
            new_row.insert(0, Repair_Type)
            new_row.insert(0, fault)

            # Append the new row to the Fault_Table
            Fault_Table = Fault_Table._append([new_row], ignore_index=True)

    Fault_Table = Fault_Table.rename(columns={0: 'Fault', 1: 'Repair_Type', 2: 'Chain_list'})
    Fault_Table.set_index('Fault')
    Fault_Table.to_csv(Fault_Table_file_name, index=True)

    return Fault_Table
            
def Repair_Statistics_using_LogicSolver(BumpMap_file_name, Fault_Type, Shorted_Bumps_Number, Short_Distance, Faults_Number, Interface_IRL_file_name, Reparability_Table_file_name, Fault_Table_file_name, Print_Fault):
    """
    This function generates repair statistics using a logic solver.
    It first generates a fault table using the Fault_Table_Generator function.
    Then, it iterates over each fault in the fault table, determines the reparability of the fault using the LogicSolver function,
    and updates the repair type in the fault table.
    Finally, it calculates and prints the repair statistics and saves the repair table to a CSV file.
    """
    # Generate the fault table using the Fault_Table_Generator function
    Fault_Table = Fault_Table_Generator(Interface_IRL_file_name, BumpMap_file_name, Faults_Number, Shorted_Bumps_Number, Short_Distance, Fault_Type, Fault_Table_file_name)

    # Load the route table and bumpmap 
    Route_Table = Repair_IRL_file_loading_into_a_dataframe(Interface_IRL_file_name)
    df_bump = Avoid_bump_name_iteration(BumpMap_file_name)  

    # Iterate over each fault in the fault table
    for index, Fault_row in Fault_Table.iterrows():
        # Extract the fault information
        fault = Fault_row['Fault']
        Chain_list = list(set(Fault_row['Chain_list']))
        Repair_Type = Fault_row['Repair_Type']

        if Print_Fault: 
            print(fault)

        # If the repair type is 'Repair', determine the reparability using the LogicSolver function
        if Repair_Type == 'Repair':
            Repair_Type = LogicSolver(Chain_list, Route_Table, df_bump, fault)

        # Update the repair type in the fault table
        Fault_row['Repair_Type'] = Repair_Type

    # Create a copy of the fault table to use as the repair table
    Repair_Table = Fault_Table.copy()

    # Calculate the number of repairable, benign, catastrophic, and unrepairable faults
    Repairable_fault = len(Repair_Table[(Repair_Table['Repair_Type'] == 'Repairable')])
    Benign_fault = len(Repair_Table[(Repair_Table['Repair_Type'] == 'Benign')])
    Catastrophic_fault = len(Repair_Table[(Repair_Table['Repair_Type'] == 'Catastrophic')])
    Unrepairable_fault = len(Repair_Table[(Repair_Table['Repair_Type'] == 'Unrepairable')])
    Total_fault = len(Repair_Table)

    # Calculate the reparability percentage
    Reparability_percentage = (Repairable_fault + Benign_fault) / Total_fault * 100
    print(f'Repair Statistics using LogicSolver : Total faults : {Total_fault} , Repairable faults : {Repairable_fault}, Benign faults :  {Benign_fault}, Catastrophic faults : {Catastrophic_fault}, Unrepairable faults : {Unrepairable_fault}, {Reparability_percentage}%')

    # Save the Repair_Table to a CSV file
    Repair_Table.to_csv(Reparability_Table_file_name, index=True)

    # Return the Repair_Table DataFrame
    return Repair_Table

def Repair_Solutions_using_RecursiveSolver(BumpMap_file_name, Interface_IRL_file_name, Repair_Solutions_Table_file_name, Faults_Number, Shorted_Bumps_Number, Short_Distance, Fault_Type, Fault_Table_file_name, Print_Fault):
    """
    This function generates repair solutions using a recursive solver. It first generates a fault table using the Fault_Table_Generator function.
    Then, it iterates over each fault in the fault table, determines the reparability of the fault using the RecursiveSolver function,
    and updates the repair type in the fault table. Finally, it calculates and prints the repair statistics and saves the repair solutions table to a CSV file.
    """

    # Generate the fault table using the Fault_Table_Generator function
    Fault_Table = Fault_Table_Generator(Interface_IRL_file_name, BumpMap_file_name, Faults_Number, Shorted_Bumps_Number, Short_Distance, Fault_Type, Fault_Table_file_name)

    # Load the route table and bumpmap
    Route_Table = Repair_IRL_file_loading_into_a_dataframe(Interface_IRL_file_name)
    df_bump = Avoid_bump_name_iteration(BumpMap_file_name) 
    
    # Initialize an empty DataFrame to store repair solutions
    Repair_Solutions_Table = pd.DataFrame()

    # Define a list of functional types that require routing
    Functionnal_type_list = ['DATA', 'ADDR', 'SIDEBAND', 'CLK']

    # Iterate over each row in the fault table
    for Fault_index, Fault_row in Fault_Table.iterrows():
        
        Route_Table_copy = Route_Table.copy()
        new_row = []
        Chain_list = list(set(Fault_row['Chain_list']))
        fault = Fault_row['Fault']
        Repair_Type = Fault_row['Repair_Type']

        if Print_Fault:     
            print(fault)

        # Check if a repair action is needed
        if Repair_Type == 'Repair':
        # if fault == ['rxdata32_phy', 'rxvldRD_phy'] or fault == ['VSS_phy_18', 'rxdata0_phy']:

            # Initialize a dictionary to map Physical Functional Sources (PFS) to their possible routes
            PFS_to_route_dict = defaultdict(list)

            # Iterate over each connection in the fault list
            for connection in fault:

                # Get the bump details from the DataFrame
                bump = df_bump[df_bump['Name'] == connection]
                # Extract the signal name by removing the '_phy' suffix
                bump_signal = bump['Name'].values[0].replace('_phy','')
                # Initialize a list to store the PFS that need to be routed
                PFS_to_route_list = []

                # Check if the bump type is in the list of functional types and is not a spare
                if bump['Type'].values[0] in Functionnal_type_list and bump['Spare'].values[0] != True:
                    # Get the repair chain associated with the faulty bump
                    faulty_bump_RepairChain = Route_Table[(Route_Table['Connection'] == bump['Name'].values[0])]['RepairChain'].values[0]
                    # Get the PFS to reroute for the repair chain
                    PFS_to_route = Route_Table[(Route_Table['RepairChain'] == faulty_bump_RepairChain)]

                    # Iterate over each PFS to route
                    for index, PFS in PFS_to_route.iterrows():
                        # Check if the connection is not a spare and is not already in the PFS_to_route_list
                        if df_bump.loc[df_bump['Name'] == PFS['Connection'], 'Spare'].values[0] != True and PFS['Signal'] not in PFS_to_route_list:
                            # Add the signal to the PFS_to_route_list
                            PFS_to_route_list.append(PFS['Signal'])

                    # If the bump signal is in the PFS_to_route_list and is in the second half of the list, reverse the list, to accelerate the solver. 
                    if bump_signal in PFS_to_route_list and PFS_to_route_list.index(bump_signal) > len(PFS_to_route_list)/2:
                        PFS_to_route_list = list(reversed(PFS_to_route_list))

                    # Add the PFS_to_route_list to the PFS_to_route_dict with the repair chain as the key
                    PFS_to_route_dict[faulty_bump_RepairChain] = PFS_to_route_list

                # Filter the Route_Table to exclude the current bump connection
                Route_Table_copy = Route_Table_copy[(Route_Table_copy['Connection'] != bump['Name'].values[0])]

            # Initialize a list to store unique PFS that need to be routed
            PFS_filter_list = []
    
            # Iterate over each repair chain and its associated PFS list in the PFS_to_route_dict
            for RepairChain, PFS_list in PFS_to_route_dict.items():
                for PFS in PFS_list:
                    # If the PFS is not already in the PFS_filter_list, add it
                    if PFS not in PFS_filter_list:
                        PFS_filter_list.append(PFS)

            # Filter Route_Table with the PFS we want to repair
            Route_Table_copy = Route_Table_copy[(Route_Table_copy['Signal'].isin(PFS_filter_list))]

            # Initialize a list to store repair chains
            Solution_Total = []
            Repair_Flag = True
            Repair_Type = 'Unrepairable'

            # Iterate over each repair chain and its associated PFS list in the PFS_to_route_dict
            for RepairChain, PFS_to_route_list in PFS_to_route_dict.items():
           
                if Repair_Flag == True:
                    
                    # Initialize a dictionary to store all possible solutions for the repair chain
                    Solutions_dict_per_RepairChain = defaultdict(list)

                    # Initialize a DataFrame to store the current solution
                    Solution = pd.DataFrame(columns=['Signal', 'Connection', 'Mux', 'Sel', 'Status'])

                    # Initialize lists to keep track of used Interconnect Sources (IS) and routed Physical Functional Sources (PFS)
                    Used_IS_list = []
                    Routed_PFS_list = []

                    # Start the recursive solver with an empty list                       
                    # This will initiate the recursive process to find all possible routing solutions for the given PFS list.
                    RecursiveSolver(Used_IS_list, Routed_PFS_list, PFS_to_route_list, 0, Solution, Solutions_dict_per_RepairChain, Route_Table_copy, Repair_Type)
         
                    # Check if any solutions were found
                    # If no solutions are found, mark the fault as unrepairable.
                    if len(Solutions_dict_per_RepairChain) == 0:
                        Repair_Flag = False
                        Repair_Type = 'Unrepairable'
                    else:
                        Repair_Type = 'Repairable'
                        Repair_Solution = [] 

                        # Iterate over each solution in the Solutions_dict_per_RepairChain
                        # This loop will extract the mux and sel values from each solution.
                        for Solution_key, Solution_value in Solutions_dict_per_RepairChain.items():
                            # Iterate over each mux and sel in the solution value
                            # This nested loop will append each mux and sel pair to the Repair_Solution list.
                            for mux, sel in Solution_value.items():
                                # Append the mux and sel to the new solution list
                                Repair_Solution.append([mux, sel])
                            # Insert the repair chain at the beginning of the new solution list
                            # This will ensure that the repair chain is associated with the solution.
                            Repair_Solution = [RepairChain, Repair_Solution]
                        Solution_Total.insert(0, Repair_Solution)

            # Insert the Solution for all the repair chain in the new row
            new_row.insert(0, Solution_Total)
        
        # Insert fault informations in the new row 
        new_row.insert(0, set(Chain_list))
        new_row.insert(0, Repair_Type)
        new_row.insert(0, fault)
       
        # Build the new table row per row
        Repair_Solutions_Table = Repair_Solutions_Table._append([new_row], ignore_index=True)
    
    # Rename and set the index to the 'Fault' column
    Repair_Solutions_Table = Repair_Solutions_Table.rename(columns={0: 'Fault', 1: 'Repair_Type', 2: 'Chain_list', 3: 'Repair_Solutions'})
    Repair_Solutions_Table.set_index('Fault')

    # Calculate the total number of faults
    Total_fault = len(Repair_Solutions_Table['Repair_Type'])

    # Calculate the number of repairable faults
    Repairable_fault = len(Repair_Solutions_Table[(Repair_Solutions_Table['Repair_Type'] == 'Repairable')])

    # Calculate the number of benign faults
    Benign_fault = len(Repair_Solutions_Table[(Repair_Solutions_Table['Repair_Type'] == 'Benign')])

    # Calculate the number of catastrophic faults
    Catastrophic_fault = len(Repair_Solutions_Table[(Repair_Solutions_Table['Repair_Type'] == 'Catastrophic')])

    # Calculate the number of unrepairable faults
    Unrepairable_fault = len(Repair_Solutions_Table[(Repair_Solutions_Table['Repair_Type'] == 'Unrepairable')])

    # Calculate the reparability percentage
    Reparability_percentage = (Repairable_fault + Benign_fault) / Total_fault * 100

    # Save the Repair_Solutions_Table to a CSV file
    Repair_Solutions_Table.to_csv(Repair_Solutions_Table_file_name, index=True)

    # Print the repair statistics
    print(f'RecursiveSolver : Total faults : {Total_fault} , Repairable faults : {Repairable_fault}, Benign faults :  {Benign_fault}, Catastrophic faults : {Catastrophic_fault}, Unrepairable faults : {Unrepairable_fault}, {Reparability_percentage}%')
    
    # Return the Repair_Solutions_Table DataFrame
    return Repair_Solutions_Table

# Section 5 : Yield and Cost Analysis
def MetaCIRA(BumpMap_file_name, Interface_IRL_file_name, System_description_file_name, System_Analysis, Min_Yield, Max_Yield, Number_of_faults_tested, Number_of_electrical_yield_tested, Bundle_Flag, Log_Scale, seed=None):

    fault_reparation = {
        'POWER': 'Benign',
        'GND': 'Benign',
        'DATA': 'Repair', 
        'CLK': 'Repair', 
        'ADDR': 'Repair', 
        'SIDEBAND': 'Repair', 
        'SPARE': 'Benign', 
        'NONE': 'Benign'}

    if seed is not None:
        random.seed(seed)

    # This function classifies faults based on the number of faults tested and the electrical yield.
    # It returns the number of benign and repairable faults.
    def Fault_Classifier(N, Number_of_faults_tested, Electrical_Yield, Route_Table, df_bump):
        # Initialize counters for benign and repairable faults.
        RepairCounter = 0
        BenignCounter = 0
        
        # Calculate the number of connections per fault.
        Nc = round((1-Electrical_Yield) * N, 8)
        # Integrer part of Nc
        A = int(Nc)
        # Decimal part of Nc
        a = round((Nc - A), 8)
      
        # Calculate the number of faults with A and A+1 connections.
        Nsup = int(Number_of_faults_tested * a) 
        Ninf = Number_of_faults_tested - Nsup 

        # Generate combinations of faulty connections.
        Faulty_Combinations = []

        # Generate combinations with A connections.
        for _ in range(Ninf):
            Faulty_Index = random.sample(list(range(N)), A) 
            Faulty_Combinations.append(Faulty_Index) 

        # Generate combinations with A+1 connections.
        for _ in range(Nsup):
            Faulty_Index = random.sample(list(range(N)), A + 1)
            Faulty_Combinations.append(Faulty_Index)  

        # Shuffle the combinations to ensure randomness.
        random.shuffle(Faulty_Combinations)

        # Iterate over each combination of faulty connections.
        for Combination in Faulty_Combinations:
            # Get the faulty connections from the DataFrame.
            Faulty_Connections = df_bump.iloc[Combination].copy()
            Chain_list = []
            Route_Table_copy = Route_Table.copy()
            fault = []
            Repair_Type = 'Benign'
            RepairFlag = False 

            # If Bundle_Flag is True, use BundleSolver to determine reparability.
            if Bundle_Flag:
                # Iterate through every connections in the fault.
                for index, bump in Faulty_Connections.iterrows(): 
                    # Add the bump name to the fault list.
                    fault.append(bump['Name'])
                    # Check if the fault type requires repair.
                    if fault_reparation[bump['Type']] == 'Repair':
                        RepairFlag = True
                # If any fault requires repair, determine the reparability using BundleSolver.
                if RepairFlag == True:
                    Repair_Type = BundleSolver(df_bump, fault, Route_Table)
                    # If the fault is repairable, increment the repair counter.
                    if Repair_Type == 'Repairable':
                        RepairCounter += 1 
                # If the fault is benign, increment the benign counter.
                if Repair_Type == 'Benign':
                    BenignCounter += 1

            # If Bundle_Flag is False, use LogicSolver to determine reparability.
            else:
                # Iterate through each connection in the fault.
                for index, bump in Faulty_Connections.iterrows(): 
                    # Append the bump name to the fault list.
                    fault.append(bump['Name'])
                    # Generate a list containing the repair chains involved.
                    if not Route_Table[(Route_Table['Connection'] == bump['Name'])]['RepairChain'].empty:
                        RepairChain = Route_Table[(Route_Table['Connection'] == bump['Name'])]['RepairChain'].values[0]
                        Chain_list.append(RepairChain)

                    # If the bump is a spare, update its type to 'SPARE'.
                    if bump['Spare'] == True:
                        bump = bump.copy()
                        bump['Type'] = 'SPARE'

                    # If the bump does not have a default route and is not of type 'DATA', update its type to 'NONE'.
                    if Route_Table_copy[(Route_Table_copy['Connection'] == bump['Name']) & (Route_Table_copy['Status'] == 'Default')].empty and bump['Type'] != 'DATA':
                        bump = bump.copy()
                        bump['Type'] = 'NONE'

                    # If the bump type requires repair, set the RepairFlag to True.
                    if fault_reparation[bump['Type']] == 'Repair':
                        RepairFlag = True

                # If the RepairFlag is True, determine the reparability using the LogicSolver function.
                if RepairFlag == True:  
                    # Remove duplicate repair chains from the Chain_list.
                    Chain_list = list(set(Chain_list))
                    Repair_Type = LogicSolver(Chain_list, Route_Table, df_bump, fault)

                    # If the fault is repairable, increment the repair counter.
                    if Repair_Type == 'Repairable':
                        RepairCounter += 1 

                # If the fault is benign, increment the benign counter.
                if Repair_Type == 'Benign':
                    BenignCounter += 1

        # Return the number of benign and repairable faults.
        return BenignCounter, RepairCounter

    if Log_Scale:
        yield_range = [1- 10 ** (-exp) for exp in range(1, Number_of_electrical_yield_tested + 1)]
    else:
        yield_range = np.linspace(Min_Yield, Max_Yield, num = Number_of_electrical_yield_tested + 1)

    print(f'Yield range : {yield_range}')    
    yield_without_repair_list = []
    yield_with_repair_list = []
    Surface_ratio_list = []
    Total_Surface = 0
    Total_Surface_repair = 0 

    # If System_Analysis is True, read the system description from the file and process each die in the system.
    if System_Analysis: 
        # Open the system description file and load its contents.
        with open(System_description_file_name, 'r') as file:
            System_description = yaml.safe_load(file)

        # Iterate over each die in the system description.
        for Die in list(System_description.keys()):
            # Extract the number of dies, interfaces, and resources for the current die.
            Die_Number = System_description[Die]['Die_Number']
            Interface_Number = System_description[Die]['Interface_Number']
            Ressources = System_description[Die]['Ressources']
            Die_Surface = Ressources['Surface']
            Interface_BumpMap_file_name = System_description[Die]['BumpMap_file_name']

            # Load the bump map data for the current interface.
            Interface_df_bump = Avoid_bump_name_iteration(Interface_BumpMap_file_name) 

            # Calculate the maximum and minimum X and Y coordinates of the bumps.
            max_X = max(Interface_df_bump['X']) # In µm²
            min_X = min(Interface_df_bump['X'])
            max_Y = max(Interface_df_bump['Y']) # In µm²
            min_Y = min(Interface_df_bump['Y'])

            # Calculate the surface area of the interface in mm².
            Interface_Surface = (max_X - min_X) * (max_Y - min_Y) * 10**-6

            # Update the total surface area and the total surface area for repair.
            Total_Surface += Die_Surface * Die_Number
            Number_Spares = len(Interface_df_bump[Interface_df_bump['Spare'] == True])
            Surface_repair = (Number_Spares/len(Interface_df_bump)) * Interface_Surface * Interface_Number
            Total_Surface_repair += Surface_repair

    for Electrical_Yield in yield_range:  # Iterate over each electrical yield value in the yield range

        print(Electrical_Yield)  # Print the current electrical yield value for debugging purposes

        if System_Analysis:  # If system analysis is enabled

            System_yield_with_repair = 1  # Initialize the system yield with repair to 1
            System_yield_without_repair = 1  # Initialize the system yield without repair to 1

            for interface in list(System_description.keys()):  # Iterate over each interface in the system description

                Interface_BumpMap_file_name = System_description[interface]['BumpMap_file_name']  # Get the bump map file name for the current interface
                Interface_IRL_file_name = System_description[interface]['IRL_file_name']  # Get the IRL file name for the current interface

                Interface_df_bump = Avoid_bump_name_iteration(Interface_BumpMap_file_name)  # Load the bump information from the bump map file into a DataFrame
                Interface_Route_Table = Repair_IRL_file_loading_into_a_dataframe(Interface_IRL_file_name)  # Load the route table from the IRL file into a DataFrame

                N = len(Interface_df_bump)  # Get the number of bumps in the current interface

                BenignCounter, RepairCounter = Fault_Classifier(N, Number_of_faults_tested, Electrical_Yield, Interface_Route_Table, Interface_df_bump)  # Classify the faults and get the number of benign and repairable faults

                yield_without_repair = BenignCounter/Number_of_faults_tested  # Calculate the yield without repair for the current interface
                yield_with_repair = (RepairCounter + BenignCounter)/Number_of_faults_tested  # Calculate the yield with repair for the current interface

                System_yield_without_repair = System_yield_without_repair * yield_without_repair  # Update the system yield without repair
                System_yield_with_repair = System_yield_with_repair * yield_with_repair  # Update the system yield with repair

            yield_without_repair_list.append(System_yield_without_repair)  # Append the system yield without repair to the list
            yield_with_repair_list.append(System_yield_with_repair)  # Append the system yield with repair to the list

            Wasted_Surface = (1 - System_yield_without_repair) * Total_Surface  # Calculate the wasted surface due to faults without repair
            Surface_ratio = Wasted_Surface / Surface_repair  # Calculate the ratio of wasted surface to repair surface
            Surface_ratio_list.append(Surface_ratio)  # Append the surface ratio to the list

        else:
            Interface_df_bump = Avoid_bump_name_iteration(BumpMap_file_name)  
            Interface_Route_Table = Repair_IRL_file_loading_into_a_dataframe(Interface_IRL_file_name)
            Number_Spares = len(Interface_df_bump[Interface_df_bump['Spare'] == True])

            max_X = max(Interface_df_bump['X']) # In µm²
            min_X = min(Interface_df_bump['X'])
            max_Y = max(Interface_df_bump['Y']) # In µm²
            min_Y = min(Interface_df_bump['Y'])
            Interface_Surface = (max_X - min_X) * (max_Y - min_Y) * 10**-6 # In mm²

            # Surface_repair = (Number_Spares/len(Interface_df_bump)) * Interface_Surface

            N = len(Interface_df_bump)
    
            BenignCounter, RepairCounter = Fault_Classifier(N, Number_of_faults_tested, Electrical_Yield, Interface_Route_Table, Interface_df_bump)
            yield_without_repair = BenignCounter / Number_of_faults_tested
            yield_with_repair = (RepairCounter + BenignCounter) / Number_of_faults_tested
            yield_without_repair_list.append(yield_without_repair)
            yield_with_repair_list.append(yield_with_repair)

            # Wasted_Surface = Interface_Surface * (1 - yield_without_repair)
            # Surface_ratio = Wasted_Surface / Surface_repair
            # Surface_ratio_list.append(Surface_ratio) 

    # Print the yield without repair and yield with repair lists for debugging purposes
    print(f'Interface yield without repair action : {yield_without_repair_list}')
    print(f'Interface yield with repair action : {yield_with_repair_list}')

    # Create a new figure for plotting
    plt.figure()

    # Scatter plot for yield with repair
    plt.scatter(yield_range, yield_with_repair_list, label='With repair', color='blue', marker='+')

    # Scatter plot for yield without repair
    plt.scatter(yield_range, yield_without_repair_list, label='Without repair', color='red', marker='x')

    # Set the x-axis to use scientific notation without offset
    formatter = ScalarFormatter(useOffset=False)
    formatter.set_scientific(True)
    plt.gca().xaxis.set_major_formatter(formatter)

    # Set the x-axis label
    plt.xlabel('Electrical Yield')

    # Set the y-axis label
    plt.ylabel('System Yield')

    # Add a legend to the plot
    plt.legend()

    # Add a grid to the plot for better readability
    plt.grid()

    # Set the title of the plot
    plt.title('System yield vs Electrical yield, with and without repair')
        
    # plt.figure()
    # plt.scatter(yield_range, Surface_ratio_list, label = 'Surface Ratio', color = 'green', marker = 'o')

    # formatter = ScalarFormatter(useOffset=False)
    # formatter.set_scientific(True)
    # plt.gca().xaxis.set_major_formatter(formatter)

    # plt.xlabel('Electrical yield')
    # plt.ylabel('Surface wasted without repair on surface used for repair')

    # plt.legend()
    # plt.grid()
    # plt.title('Interest of reparability vs electrical yield')

    plt.savefig("plot.pdf", format="pdf")
    print('Wrote graph to plot.pdf')
    plt.savefig("plot.svg", format="svg")
    print('Wrote graph to plot.svg')
    plt.show()


if Create_SVG:
    Display_SVG(BumpMap_file_name, Aspect_file_name, BumpMap_SVG_image_file_name, Open_SVG, Bump_Diameter, Pitch, 
    Input_X_scale, Input_Y_scale, Legend, Margin, Bump_Name, Stroke_Color, Font, Font_Size, Display_Reparability_SVG)
 
if Reparability_Statistics:
    Repair_Statistics_using_LogicSolver(BumpMap_file_name, Fault_Type, Shorted_Bumps_Number, Short_Distance, 
    Faults_Number, Interface_IRL_file_name, Reparability_Table_file_name, Fault_Table_file_name, Print_Fault)

if Repair_Solutions:
    Repair_Solutions_using_RecursiveSolver(BumpMap_file_name, Interface_IRL_file_name, Repair_Solutions_Table_file_name, 
    Faults_Number, Shorted_Bumps_Number, Short_Distance, Fault_Type, Fault_Table_file_name, Print_Fault)

if Meta_Analysis:
    MetaCIRA(BumpMap_file_name, Interface_IRL_file_name, System_description_file_name, System_Analysis, Min_Yield, 
    Max_Yield, Number_of_faults_tested, Number_of_electrical_yield_tested, Bundle_Flag, Log_Scale, seed = None)

end = time.time()
print(f'Execution time = {end - start} s')    
