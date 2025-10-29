---
editor_options: 
  markdown: 
    wrap: 72
---

# Agricultural Loss Estimation 

## Overview

This repository contains a statistical analysis of agricultural
production losses from climate shocks (drought and flood) as a function
of farm size. The analysis builds upon a
[dataset](https://github.com/Better-Planet-Laboratory/microshock)
developed by Zia Mehrabi, Julie Fortin and Navin Ramankutty, which
harmonized household survey data across 17 countries to examine the
relationship between farm size and climate shock impacts.

## Data Source

The analysis uses the `microdata_spei_merge.RDS` dataset from this prior
project, which combines:

-   **Farm size data**: Extracted from household surveys and
    agricultural censuses across 17 countries (Albania, Bangladesh,
    Burkina Faso, Colombia, Ethiopia, Guatemala, India, Mali, Malawi,
    Mexico, Niger, Nigeria, Pakistan, Tanzania, Timor-Leste, Uganda)
-   **Climate shock data**: Binary and continuous measures of production
    losses from droughts and floods reported by each household
-   **SPEI data**: Standardized Precipitation Evapotranspiration Index
    data for validation of climate anomalies during survey periods and
    for use in thresholding effects
-   **Geographic information**: Administrative units at multiple
    hierarchical levels for spatial analysis

### Countries and Surveys Included

The dataset encompasses: - 17 countries across Africa, Asia, and Latin
America - Multiple household survey types (LSMS, agricultural censuses,
specialized surveys) - Survey years ranging from 2000-2015 - Over
100,000 households with agricultural production data (\>1M households
are available in entirely but 1M are from one country alone, Colombia).

## Analysis Structure (Loss_calc.R)

### Data Preparation

#### Farm Size Distribution Analysis

The script begins with exploratory analysis of farm size distributions
across countries from the surveys in order to inspect data quality,
using the `getdistros()` function, which calculates: - Weighted median
and mean farm sizes by country - 10th and 90th percentiles - Min/max
farm size. Several data quality issues were identified and corrected:

1.  **Guatemala outlier correction**: Removed farms ≥1000 ha due to
    clear data errors
2.  **Farm size scaling for specific countries**: Applied scaling
    corrections for Bangladesh, Ethiopia, and Burkina Faso based on FAO
    World Census of Agriculture data:
    -   Bangladesh: Scaled by ratio to 0.52 ha mean
    -   Ethiopia: Scaled by ratio to 1.03 ha mean
    -   Burkina Faso: Scaled by ratio to 4.18 ha mean

#### Variable Transformations

Data were transformed prior to modelling.

-   **Farm size**: Log-transformed (LogFSIZE) to handle skewed
    distributions, with zero values replaced by 0.0001
-   **Household weights**: Rounded to integers and NA/zero values
    replaced with 1.
-   **Geographic hierarchy**: Smallest Sampling Unit (SSU) determined
    using `afn.findssu()` function. Note this is different from the
    Smallest Administrative Unit (SRAU) for which surveys were designed.

### Country Coverage by Analysis Type

#### Binary Loss Analysis (Production Loss Occurrence)

**Countries included**: Albania, Bangladesh, Burkina Faso, Ethiopia,
Guatemala, India, Mali, Malawi, Mexico, Niger, Nigeria, Pakistan,
Tanzania, Timor-Leste, Uganda *Note: Colombia excluded due to sample
size - alternative sampling approach available*

#### Continuous Loss Analysis (Revenue Loss Percentage)

**Countries with drought data**: Burkina Faso, India, Mali, Niger,
Pakistan **Countries with flood data**: Niger, Pakistan *Note:
Continuous loss analysis requires specific survey questions about
harvest amounts and losses, available in fewer countries*

### Statistical Models

The analysis implements multiple statistical approaches to model the
relationship between farm size and production losses:

#### 1. Binary Loss Models

**Model Specification:**

```         
Loss ~ LogFSIZE * Event + (1|LVL0/SRAU)
```

**Mathematical Form:**

```         
logit(P(Loss = 1)) = β₀ + β₁·LogFSIZE + β₂·Event + β₃·(LogFSIZE × Event) + u₀ + u₁
```

Where: - `Loss`: Binary indicator (0/1) for production loss reporting
due to a flood or drought and - `LogFSIZE`: Natural logarithm of farm
size in hectares - `Event`: Factor variable (Drought vs. Flood) -
`u₀ ~ N(0, σ²₀)`: Country-level random effect - `u₁ ~ N(0, σ²₁)`:
SRAU-level random effect nested within country

**Implementation:** - **Frequentist approach**: `glmer()` with binomial
family and logit link - **Bayesian quantile based approach**: `brm()`
with Bernoulli family

#### 2. Continuous Loss Models

**Model Specification:**

```         
Lost.Revenue ~ LogFSIZE * Event + (1|LVL0/SRAU/SSU)
```

**Mathematical Form:**

```         
logit(E[Lost.Revenue]) = β₀ + β₁·LogFSIZE + β₂·Event + β₃·(LogFSIZE × Event) + u₀ + u₁ + u₂
```

Where: - `Lost.Revenue`: Beta-distributed revenue loss proportion
(0,1) - Three-level nested random effects: Country/SRAU -
`u₀, u₁ ~ N(0, σ²ᵢ)`

**Revenue Loss Calculation:** Original loss proportions (0, 1) were
transformed to avoid boundary issues: - 0 → 0.00001 - 1 → 0.9999

**Implementation:** - **Primary model**: `glmmTMB()` with beta family
and logit link - **Alternative approach**: `qgamV()` for quantile
regression at median (0.5 quantile)

### Model Diagnostics

Diagnostics are performed on the non-quantile based models
`simulateResiduals()` for distributional assumption checking; outlier
detection via boostrapping, model validation using inbuild functions.
Generally while some issues with models appear to be present, quantile
based approaches produce very similar results, and so any issues
unlikely to be problematic for use.

### SPEI thresholding

Models are re-run after filtering for extreme SPEI min (\< -1.49) and
max (\> 1.49) for drought and flood, respectively. These models can be
used for downstream modelling attempting to address the farm size loss
relationships under extreme wet and dry conditions. Min and max spei
values are computed relative to the growing season and survey window of
each input dataset (they are the monthly min/max over the full
concatenated set of growing seasons referred to in the survey question,
for example if the survey asked "Did you experience a loss due to
drought in the last 12 months" the max spei would be the max monthly
value in the growing season for that period).

## Key Findings and Model Outputs

### Model Results Storage:

-   `model_binary_con.RDS`: Binary loss GLM results
-   `model_con.RDS`: Continuous loss GLM results
-   `brms_bin.RDS`: Bayesian binary model results
-   `cont_qgam.RDS`: Median regression results
-   `model_binary_th.RDS`: Binary loss GLM results thresholded
-   `model_con_th.RDS`: Continuous loss GLM results thresholded

### Visualization Outputs:

-   `model_binary_con.png`: Predicted probability of loss by farm size
    and event type

-   `model_con.png`: Predicted revenue loss proportion by farm size and
    event type

-   `model_bin_brms.png`: Bayesian model predictions

-   `model_con_gam.png`: Quantile regression predictions

-   `model_binary_th.png`: Predicted probability of loss by farm size
    and event type thresholded

-   `model_con_th.png`: Predicted revenue loss proportion by farm size
    and event type thresholded

## Key Results and Visualizations

### Binary Loss Models: Probability of Production Loss

![Binary Loss Model - Frequentist Approach](figs/model_binary_con.png)

*Figure 1: Predicted probability of experiencing production loss by farm
size (log scale) for drought and flood events. Based on mixed-effects
logistic regression*

![](figs/model_bin_brms.png)

*Figure 2: Bayesian median based estimates of loss probability by farm
size and event type, providing robust estimates less sensitive to
outliers.*

### Continuous Loss Models: Revenue Loss Percentage

![](figs/model_con.png)

*Figure 3: Predicted percentage of revenue lost due to climate shocks by
farm size. Based on beta regression for countries with detailed
harvest/loss data.*

![](figs/model_con_gam.png)

*Figure 4: Median regression estimates of revenue loss, providing robust
estimates less sensitive to outliers.*

## Citation

Mehrabi, Zia. Agricultural Loss Estimation. Better Planet Laboratory.
2025.

### 
