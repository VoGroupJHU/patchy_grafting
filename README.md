Here, we provide a set of sample scripts for the simulation component

All codes run on MATLAB version R2018 or higher. MATLAB installation instructions can be obtained from the software provider.

The main script is called "patchy_grafting.m"

Running "patchy_grafting.m" perform the implimented scaling theory and place chain on the prism surface based on the computed grafting probabilities as described in the main text and save the data to a MATLAB mat binary file. Parameters such as the defined grafting density (sigma), chain length (N), or temperature (kT) can be changed in the script. The code should complete within at most 5 minutes

An example out "Prism3_sigma_0.1_N_50_kT_1.mat" is provided. This file can be loaded into matlab for plotting/analysis. 

There are a couple other included convenient functions defined as followed:

Compute the intertial tensor of a given shape
calc_inertial_tensor.m 			

Compute a surface mesh for a given input shape
gen_poly.m 							    

Compute an object that evaluate the kernel shape parameter
parameterize_shape.m 				

Converts from cartesian to spherical coordinates 
xyz_to_kernel.m 					  

Convenient data structure of pre-computed surface mesh
core_data.mat 						  

Coorindates defining vertices of the prism
Prism3.txt 							    
