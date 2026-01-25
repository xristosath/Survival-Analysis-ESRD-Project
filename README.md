## Project Overview
This project analyzes the survival outcomes of 8,796 elderly patients with End-Stage Renal Disease (ESRD), comparing those who initiated Haemodialysis against those who chose conservative management. The primary goal was to estimate the relative effectiveness of treatment while adjusting for confounding factors like age, sex, and comorbidities.

## Methodology
The analysis follows a rigorous statistical pipeline implemented in R:


Descriptive Stats: Crude rates and incidence rate ratios (Poisson approach).


Non-Parametric Analysis: Kaplan-Meier estimation and Log-Rank/Peto-Peto tests to compare survival curves.


Cox Proportional Hazards (PH) Model: Univariate and Multivariable modeling.

Assumption Checking:


Linearity: Martingale residuals (confirmed linear fit for Age).


Proportional Hazards: Scaled Schoenfeld residuals and Log-Log plots.

Addressing PH Violations: The analysis revealed that the PH assumption was violated for treatment. This was resolved using a Piecewise Cox Model and Time-Dependent Coefficients, allowing the Hazard Ratio to vary between the early phase (0-1 year) and late phase (1+ years).


Parametric Modeling: Fitted a Weibull model and compared it against the semi-parametric results.

## Key Findings

Crude estimates were misleading: Simple relative risk calculations underestimated the benefit due to differing follow-up times.

Time-Dependent Effect: The effectiveness of Haemodialysis is not constant. It offers a massive protective effect in the first year (Adjusted HR: 0.15, ~85% risk reduction). However, this benefit diminishes significantly after year 1, likely due to survivor bias in the conservative group.

## Libraries Used
survival

tidyverse 
