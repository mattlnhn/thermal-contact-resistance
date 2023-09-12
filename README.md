# Thermal Contact Resistance

## IHCP solver based on the function specification method to calculate unknown heat transfer coefficients

Explanation of the maths behind this code can be found in the .pdf.

prepdata.m is a function to deal with errors in data gathering. It removes duplicate rows which would appear to give zero heat flux.

Requirements:

0. Clone repo

1. Add column headings to data files: time, pos, load; T_oven, T_Cu1, T_Cu2, T_Inco1, T_H25, T_Inco2, T_Cu3, T_Cu4

2. Run ihcp_exported.m

    a) Select file

    b) Set parameters, settings

    c) Run

