# Models, data, and graphical results
## for draft publication titled A data assimilation framework for improved pumping-induced subsidence estimation in data-scarce settings."
## Authors: Jenny T. Soonthornrangsan<sup>1</sup>, Femke C. Vossepoel<sup>1</sup>, Mark Bakker<sup>2</sup>
##### <sup>1</sup>Department of Geoscience & Engineering, Faculty of Civil Engineering and Geoscience, Delft University of Technology, Stevinweg 1, 2628 CN Delft, The Netherlands
##### <sup>2</sup>Department of Water Management, Faculty of Civil Engineering and Geoscience, Delft University of Technology, Stevinweg 1, 2628 CN Delft, The Netherlands
<br />
<br />
<br />

Current release:  [![DOI](url)

Libraries:
`Pastas`: 1.3.0
`Pandas`: 2.0.3
`Numpy`: 1.26.3

Various python scripts are provided that create different graphical results.

- `Figures.py`: Produces the figures shown in the main text of the paper

- `SubsidenceModelResults_1978-2020.py`: Computes land subsidence at each well nest location. Plots annual subsidence rates (cm) bar graphs during 1978-2020 (Shown in the main text and supplemental information), error maps (subsidence RMSE for each well nest), forecast maps (2020-2060 in cm), sensitivity analysis (LCBKK013 for parameter sets Sske (clay), Sske (sand), Sskv (clay), K (clay), thickness; long run time so only calculating for one well nest at a time) (Shown in Supplemental Information)
  - Notebook version: https://githubtocolab.com/jsoontho/BKKSubPastasModels/blob/main/JupyterNotebooks/SubsidenceModelResults_1978-2020.ipynb

- `Pastas_ModelGraphs_1950-2020.py`: Creates Pastas models with the option to save and import the model as well as produces graphical results shown in the paper and supplemental information. Models simulate groundwater for each well nest
  - Notebook version: https://githubtocolab.com/jsoontho/BKKSubPastasModels/blob/main/JupyterNotebooks/Pastas_ModelGraphs_1950-2020.ipynb
  
- `Pastas_ResultsMaps_1950-2020.py`: Creates spatial maps that show the groundwater RMSE and t<sub>90</sub> results for each well in each well nest. Imports previously created Pastas models
  - Notebook version: https://githubtocolab.com/jsoontho/BKKSubPastasModels/blob/main/JupyterNotebooks/Pastas_ResultsMaps_1950-2020.ipynb


<br />
<br />
<br />

`figures\` : figures created from the scripts

`figures_test\`: figures used in publication to be compared to output figures from script. Should match.

`models\`: Pastas models created from Pastas_ModelGraphs_1950-2020.py script (Model files end with .pas)

`inputs\`: inputs needed for scripts 

- Groundwater observations for each well nest. Note that in the paper, each well nest is referred to without "LC" (`LC******.xlsx`)

- Subsidence observations for each well nest (`SurveyingLevels.xlsx`)

- Parameters required for subsidence calcultion for each well nest e.g. Sske, Sskv, K, thickness (`SUBParameters.xlsx`)

- Initial estimate of steady state heads for each well nest (`SS_Head_GWWellLocs.xlsx`)

- Land surface elevation for each well nest (`LandSurfElev_GWWellLocs.xlsx`)
   - Coastal DEM 2.1: `https://www.climatecentral.org/coastaldem-v2.1`
   - Kulp, S. A., and B. H. Strauss. 2021. CoastalDEM v2. 1: A high-accuracy and high-resolution global coastal elevation model trained on ICESat-2 satellite lidar.

- Location, well depth, screen depth, available observation years for each well for each well nest (`GroundwaterWellLocs.xls`)

- Basin-wide pumping data 1954-2060 (`BasinPumping.xlsx`)
  - 500,000 m<sup>3</sup>/day scenario: EstTotalPump_54-60_Int50 sheet
  - 250,000 m<sup>3</sup>/day scenario: EstTotalPump_54-60_IntF25 sheet
  - 1,000,000 m<sup>3</sup>/day scenario: EstTotalPump_54-60_IntF100 sheet
  - 500,000 to 250,000 m<sup>3</sup>/day in 2040 scenario: EstTotalPump_54-60_IntF50_25 sheet
  - No pumping scenario: EstTotalPump_54-60_IntF0 sheet



