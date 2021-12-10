
### Overview

This repository contains code associated with the publication 
**Deep Phenotyping of Alzheimer’s Disease Leveraging Electronic Medical Records Identifies Sex-Specific Associations**

Citation: TBD

---
###  Data Preparation
In order to run this code, electronic medical record data will need to be extracted and put into formats listed below. A skeleton of the data format is shown in the `\Data` folder.

**Patient identification**
* Inclusion criteria:
    * Patients diagnosed with ICD-10-CM codes G30.1, G30.8, G30.9
    * Patients age > 64

* Propensity-score matched controls:
    * All other patients in the database excluding cohort   
    * Matched on: Race, Age, Sex, Death_Status
    * Sentivity analysis controls matched on: Race, Age, Sex, Death_Status, Total_Number_of_Encounters, Total_Duration_of_Records_in_EMR
        * Controls filtered to patients with >10 encounters/visits, and >1 year of records
    * The code to identify controls is in `0-controls.R`

**Table Format**

**x_demographics.csv**

Columns: PatientID, Sex, Ethnicity, Race, BirthDate, DeathStatus, Age

**x_diagnosis.csv**

Columns: PatientID, FullDiagnosisName, ICD10_Code, Level2_Category, Level3_Category

**x_medications.csv**

Columns: PatientID, MedicationName, MedicationGenericName(dosage_route_removed)

**x_labresults.csv**

Columns: PatientID, TestName, TestResult

### Running the Code

1. Create python environment and install dependencies
`conda create --name ehr_phenotype python=3.8 pip`
`conda activate ehr_phenotype`
`pip install -r requirements.txt`

2. Identify cohort and controls. An example of propensity score matching using R for control identification is in `0-controls.R`.

3. `1-LowDimRepresent.ipynb` contains code for creating low dimensional representations for patient cohorts. 

4. `2-CompareADControls.ipynb` contains code for performing association analysis.

5. `3-ComorbidityNetworks.ipynb` contains code for creating comorbidity networks and saving as graphml files. These files can be imported into Cytoscape for visualization and analysis. 
