# internship

## IHCP solver based on the function specification method to calculate unknown heat transfer coefficients

Explanation of the maths behind this code can be found in the .pdf.

To use this code, place the data file (.dat, .csv, .txt, etc.) in this directory, change the filenames in upstreamdownstream.m and main.m, and run upstreamdownstream.m then main.m.

temp1I_PC_QBC.m (one interface, perfect conduction, heat flux boundary conditions) and temp2I_CR_QBC (two interfaces, contact resistance, heat flux boundary conditions) are functions implementing a tri-diagonal matrix method to solve the direct problem in the specified scenarios.

Column headings in the data file should be: pos, load, time, T_oven, T_Cu1, T_Cu2, T_Inco1, T_H25, T_Inco2, T_Cu3, T_Cu4.
