This code was created to calculate the clustering properties of the stars, 
specifically, the Herbig Ae/Be stars, using data from Gaia DR2. 

The clustering properties are found analyzing the distributions of the parallax and proper motions in RA and DEC. It is important, but not necessary, to know the parallax and proper motions in RA and DEC of the Herbig Ae/Be star because those ones would be used to make the selection of the data.  

The input file with the general information about the sample should contain, in the following order: star_name, starRA, starDEC, spt, para_know, parae_know, pmra_know, pmrae_know, pmdec_know, pmdece_know, g_mag_know, bp_rp_mag_know, dist_exp. 

The input file of each star can be in a VOTABLe format if the data come directly from the GAIA archive or in a text file format, following order: IDGaia, ID, ra, ra_error, dec, dec_error, para, para_error, pmra, pmra_error, pmdec, pmdec_error, g_mag, bp_rp_mag, bp_g_mag, g_rp_mag.  One of these formats needs to be select before running the code. 

The path also needs to be defined at the beginning of the code.

When the code start running would print on the screen information about the stars and the selection and also questions which answers have to be yes or no. Use 'y' for yes or 'n' for no, using the simple quotation marks.   

Plots would be created for each star after each selection in the different parameters (parallax, pmra, pmdec), which would be saved in a specific folder under the name of the star. 

At the end of each run, a file by stars would be created with the information selected. 

Also, after running each run, a final question would be print on the screen asking if "This star is in a cluster? (Yes=1; No=0; Maybe=2)" which should be answered with a number from 0 to 2. This information would be saved in the end, on a file which header would be: Name  RA DEC  PARALLAX PARALLAX_Error PMRA PMRA_Error PMDEC PMDEC_Error ClusterClassNum NumObjSelect. 


           