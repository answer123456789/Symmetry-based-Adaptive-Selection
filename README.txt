The following instructions are provided to reproduce analyses in the manuscript.

Step 1: Load the following R-packages.
knockoff 0.3.6                  MASS 7.3-56.                glmnet 4.1-8.              
ggplot2 3.5.0                    latex2exp 0.9.6.             patchwork 1.2.0
ggthemes 5.1.0.               doSNOW 1.0.20.           doParallel 1.0.17
foreach 1.5.2
 
Step 2: Run the simulations of Section 5 and Appendixes B, C, and E. Codes are in Folder “Simulation”.
(1). Check the path of .Rmd file and set it as a working directory before you run the codes.
(2). Set the file “Rmarkdown/Functions.R” under the path of .Rmd file.
(3). Run the file “Figures_1_2.Rmd” for checking Figures 1-2.
(4). Run the file “Figures_3_4.Rmd” for checking Figures 3-4.
(5). Run the files “DS_1234.Rmd” for checking Figure 5, “DS_5.Rmd” for checking Figure 6, “Knock_1234.Rmd” for checking Figure 7, “Knock_5.Rmd” for checking Figure 8 in turn.
(6). Run the files “DS_small_np.Rmd” for checking Figure B1, “Knockoff_small_np.Rmd” for checking Figure B2.
(7). Run the file “Consistency_plot.Rmd” for checking Figure C7.
(8). Run the files “DERDS_1234.Rmd” for checking Figure E8, “DERDS_5.Rmd” for checking Figure E9, “DERKnock_1234.Rmd” for checking Figure E10, “DERKnock_5.Rmd” for checking Figure E11 in turn.
 
Step 3: Run the real data analysis of Section 6 and Appendixes B and E. Codes are in Folder “Real_data”.
(1). Check the path of .Rmd file and set it as a working directory before you run the codes.
(2). Set the file “Rmarkdown/Functions.R” under the path of .Rmd file.
(3). Run the file “DER_HIV.Rmd” for checking Figures 10, B3, E12, and E13.
(4). Run the file “HIV_scatter.Rmd” for checking Figure 9 and Table B1.
(5). Run the file “HIV_stability.Rmd” for checking Figure E14.
(6). Run the file “Supermarket.Rmd” for checking Figures B4-B6 and E15-E16.
 
All the results of .Rmd files can be founded in Folder “Html”.

