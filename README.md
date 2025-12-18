# CIRA

CIRA (Chiplet Interface Repair Analyzer) is a tool which analyze the reparability of a chiplet interface. 
It uses bump map (physical description of each connection) and IRL (Interconnect Repair Langage, describe the repair option of each connection) files. 
The user can specify a fault model (open, short etc) and CIRA will ouput reparability statistics, repair solution for every possible faults, SVG representation of the interface etc. 

## Requirements

### Python Version 
This script requires Python 3.11 or higher.

### Dependencies
Install the required packages using: 

```bash
pip install -r requirements.txt
```

pandas: Data manipulation and analysis.
numpy: Numerical computing.
matplotlib: Plotting and visualization.
PyYAML: YAML file parsing for data extraction.
drawsvg: SVG drawing for chiplet interface visualization.

## Usage

In this section, several examples will by detailled. Please refer to the file DEMO\MyChipletInterface\MCI_description.txt for more details about the interfaces available in the folder DEMO. 
At the end, another interface will be used, called HYDRA. Again, please refer to the file DEMO\HYDRA\HYDRA_description.txt.

### Examples 


#### Display
To simply display an interface in an SVG image, please run :

```bash
python CIRA.py --Create_SVG --Open_SVG --Legend --Bump_Name --BumpMap_file_name .\DEMO\MyChipletInterface\MCI_1_BumpMap.yaml --BumpMap_SVG_image_file_name .\OutputFiles\BumpMap.svg --Aspect_file_name .\DEMO\colors_shapes_dict.csv --Pitch 12
```

CIRA will create the SVG in the folder OutputFiles. It will also open it. On the SVG you will find the interface, but also a legend and the bumps name. 
The aspect file is used to specify the shapes and colors of the connections according to their type. 
Finally, the --Pitch argument is used to specify the pitch of the interface, always in µm. 
If it's not specified by the user, CIRA will estimate it but it may not be precise enough.

#### Reparability Statistics
To analyze an interface for a specific fault models, please run : 

```bash
python CIRA.py  --BumpMap_file_name .\DEMO\MyChipletInterface\MCI_1_BumpMap.yaml --IRL_file_name .\DEMO\MyChipletInterface\MCI_1.irl  --Reparability_Statistics --Fault_Type 'Short' --Shorted_Bumps_Number 2 --Faults_Number 1 --Short_Distance 12 --Fault_Table_file_name .\OutputFiles\Fault_Table.yaml --Reparability_Table_file_name .\OutputFiles\Repair_Table.yaml --Print_Fault
``` 

CIRA will read the bump map and IRL file. It will analyze the interface and find every fault corresponding to one short, affecting two connections, separated by less than 12 µm. 
It will create a first file named Fault_Table.yaml, which will contains informations on all the faults. 
The affected connections and repair chain and the type of fault (Benign, Catastrophic, Repair). 
It will then create a second file named Repair_Table.yaml, containing the same information but instead of Repair in the fault type column, the user will find Repairable or Unrepairable. 
Finally, CIRA will also print the reparability statistics. 

The argument --Print_Fault is a flag, that if called, will enable CIRA to print every fault it analyzes in the terminal. It is particularly useful for debugging purpose or just to follow the progression of CIRA. 

For the moment, CIRA can analyze any numbers of open (double-open, triple-open etc).
It can also analyze 2-bump short, 3-bump short etc.
But CIRA cannot analyze combination of two (or more) shorts. 

#### Repair Solutions 
To go a step further and generate the repair solutions for each fault (ie : the state of each MUX in the affected repair chain), please run : 

```bash
python CIRA.py  --BumpMap_file_name .\DEMO\MyChipletInterface\MCI_1_BumpMap.yaml --IRL_file_name .\DEMO\MyChipletInterface\MCI_1.irl  --Repair_Solutions --Fault_Type 'Short' --Shorted_Bumps_Number 2 --Faults_Number 1 --Short_Distance 12 --Fault_Table_file_name .\OutputFiles\Fault_Table.yaml --Repair_Solutions_Table_file_name .\OutputFiles\Repair_Solutions_Table.yaml --Print_Fault
``` 
In this example we replace Reparability_Statistics by Repair_Solutions. 
Instead of having a file that only contains the reparability of each fault, the generated file (named Repair_Solutions_Table.yaml) will contains the repair solution for every repairable faults. 
The other arguments stay unchanged. 

#### Display Reparability
CIRA can also display 2-bumps shorts on the SVG representation of the interface. 
Repairable, Unrepairable, Catastrophic and Benign short will appears. 
It only works with 2-bumps short.
To display reparability, please run : 

```bash
python CIRA.py --Create_SVG --Open_SVG --Legend --Bump_Name --BumpMap_file_name .\DEMO\MyChipletInterface\MCI_1_BumpMap.yaml --IRL_file_name .\DEMO\MyChipletInterface\MCI_1.irl --BumpMap_SVG_image_file_name .\OutputFiles\BumpMap.svg --Aspect_file_name .\DEMO\colors_shapes_dict.csv --Pitch 12 --Shorted_Bumps_Number 2 --Faults_Number 1 --Short_Distance 12 --Fault_Table_file_name .\OutputFiles\Fault_Table.yaml --Reparability_Table_file_name .\OutputFiles\Repair_Table.yaml --Display_Reparability_SVG 
```
We just need to call the argument --Display_Reparability_SVG in addition to the argument need to display SVG and the argument for reparability statistics. 

#### MetaCIRA
The MetaCIRA function is able to calculate the yield of an interface for a specific electrical yield. 
For example, if we have an electrical yield of 0.999, it means 1 in 1000 connections is not able to transmit current.
From this, MetaCIRA can calculate the interface yield, for example, for 1000 interfaces, how many will be functionnal, with and without repair action. 
This function aims to demonstrate the interest of reparability as a function of electrical yield. 
To run MetaCIRA, please run : 

```bash
python CIRA.py --BumpMap_file_name .\DEMO\MyChipletInterface\MCI_1_BumpMap.yaml --IRL_file_name .\DEMO\MyChipletInterface\MCI_1.irl --Meta_Analysis --Min_Yield 0.95 --Max_Yield 1 --Number_of_faults_tested 1000 --Number_of_electrical_yield_tested 10
```
CIRA will load the bump and IRL file and proceed to the "Meta Analysis". 
You need to specify the yield range (min and max yield), Max_Yield has be superior to Min_Yield but inferior or equal to 1. Min_Yield has to be superior to 0. 
The yield range represent the electrical yield range that will be analyzed. 
Number_of_electrical_yield_tested specify the number of yield in the yield range. 
Finally, Number_of_faults_tested specify the number of faults that will be analyzed. 
For example, in an interface composed of 100 connections, with an electrical yield of 0.99, MetaCIRA will analyze the reparability of 1000 faults which affect only 1 connections in the interface. 
This function will ouput the interface yield with and without repair, in a list format but also on a plot. 

The user can add the argument --Log_Scale to represent the electrical yield in a logarithmic scale. For example : [0.9, 0.99, 0.999, 0.9999, 0.99999, 0.999999, 0.9999999]. 
To analyze the same electrical yield, please run : 

```bash
python CIRA.py --BumpMap_file_name .\DEMO\MyChipletInterface\MCI_1_BumpMap.yaml --IRL_file_name .\DEMO\MyChipletInterface\MCI_1.irl --Meta_Analysis --Number_of_faults_tested 1000 --Number_of_electrical_yield_tested 7 --Log_Scale
```

#### Bundle Repair Mechanisms (BRM)

WARNING : this only works with files made with the same architecture as the one available in the HYDRA folder. 
It won't work with the MCI files. 

One new repair strategy consist in the use of bundle. 
A bundle can be defined as a group of connection that will be shifted at the same time, instead of a repair chain where a shift is done at the connection scale.
For example, an interface composed of 100 connections could be composed of 5 bundles, each composed of 20 connections. 
If a fault affect one bundle, then the 20 corresponding connection will be shifted to a repair bundle. 
The last function of CIRA enables the user to analyze the efficiency of this new mechanism.   
This only works with MetaCIRA, on a log scale or not. 

The example in this case will be an interface called HYDRA, for HYbrid bonDing Repair Architecture. 
Please run : 

```bash
python CIRA.py --BumpMap_file_name .\DEMO\HYDRA\HYDRA_4-1_BumpMap.yaml --IRL_file_name .\DEMO\HYDRA\HYDRA_4-1_2RB.irl --Meta_Analysis --Bundle_Flag --Number_of_faults_tested 1000 --Number_of_electrical_yield_tested 7 --Log_Scale
```

CIRA will load the files and proceed to the "Meta Analysis", with 1000 faults per electrical yield and 7 electrical yield on a log scale. 
We need to add the argument Bundle_Flag to tell CIRA to treat the interface as a set of bundles.
Again, please refer to the file DEMO\HYDRA\HYDRA_description.txt

## What's next ? 
