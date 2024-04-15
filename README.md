The following instructions are provided to reproduce analyses in the manuscript.

Step 1: Load the following R-packages.
knockoff 0.3.6                  MASS 7.3-56.                glmnet 4.1-8.              
ggplot2 3.5.0                    latex2exp 0.9.6.             patchwork 1.2.0
ggthemes 5.1.0.               doSNOW 1.0.20.           doParallel 1.0.17
foreach 1.5.2
 
Step 2: Run the simulations of Section 5 and Appendix B-C. Codes are in Folder “Simulation”.
(1). Check the path of .Rmd file and set it as a working directory before you run the codes.
(2). Set the file “Rmarkdown/Functions.R” under the path of .Rmd file.
(3). Run the file “Figures_1_2.Rmd” for checking Figures 1-2.
(4). Run the file “Figures_1_2.Rmd” for checking Figures 3-4.
(5). Run the files “DS_1234_final.Rmd” for checking Figures 5 and C6, “DS_5_final.Rmd” for checking Figures 6 and C7, “Knock_1234_final.Rmd” for checking Figures 7 and C8, “Knock_5_final.Rmd” for checking Figures 8 and C9 in turn.
(6). Run the file “Consistency_plot.Rmd” for checking Figure B1.
 
Step 3: Run the real data analysis of Section 6 and Appendix C. Codes are in Folder “Real_data”.
(1). Check the path of .Rmd file and set it as a working directory before you run the codes.
(2). Set the file “Rmarkdown/Functions.R” under the path of .Rmd file.
(3). Run the file “HIV_bar_final.Rmd” for checking Figures 10, C4, C10, and C11.
(4). Run the file “HIV_scatter_final.Rmd” for checking Figure 9.
(5). Run the file “Supermarket_final.Rmd” for checking Figures 11, C5, and C12.
 
All the results of .Rmd files can be founded in Folder “Html”.

