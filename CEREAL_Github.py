#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 11 19:22:25 2020

@author: Alice Perez-Blanco


"""


#---------------------------------------------
### Packages ###

import numpy as np
import scipy
import matplotlib.pyplot as plt
import warnings
import time
from astropy.io.votable import parse_single_table
import os

#---------------------------------------------
### Path to save the data and plots ###
# Add the directorios where the figures and data would be save

pathfig  = '/PATH/'             
pathdata = '/PATH/'         
path_inputdata= '/PATH/' # This should contain the diretory where the VOTABLE (Gaia data) is save

#---------------------------------------------
### Data input ###

# This file have to contain the general information for each star. This would be the input for the known value of each target.
ClusterInfo   = scipy.genfromtxt(path_inputdata + 'FILE.txt', dtype='string',comments='#')   
                                 
# This file should be a list with the input data of each target for GAIA DR2, or other.                              
Clusterdata   = scipy.genfromtxt(path_inputdata + 'FILE.txt', dtype='string',comments='#')        

#---------------------------------------------                                
### Variables ###                              

then = time.time()                          # Time before the operations start 
warnings.filterwarnings('ignore')           # Do not print the warnings on scream 
NumStar = 1                                 # Print on scream the number of stars on the list and which is runnig at the moment
ClusterNum,Star,A1,A2,A3,A4,A5,A6,A7,A8,A9,B1 = [],[],[],[],[],[],[],[],[],[],[],[]   # Empty arrays to created a final file from  the clusters clasification

#---------------------------------------------#---------------------------------------------
### Code body ###


for z in range(len(ClusterInfo)):
    
    print '\033[1;31m'
    print 'Star '+ '\t'  + str(NumStar) + '\t' + 'of' + '\t' + str(len(ClusterInfo))
    print '\033[1;30m \n'


    # This should be the basic information you should know about the target 
    # and should be on the file save in 'ClusterInfo'
    info = scipy.genfromtxt(ClusterInfo[z], dtype='string',comments='#')

    
    star_name       = info[0]       
    starRA          = info[1].astype(float)
    starDEC         = info[2].astype(float)
    spt             = info[3]
    para_know       = info[5].astype(float)
    parae_know      = info[6].astype(float)
    pmra_know       = info[7].astype(float)
    pmrae_know      = info[8].astype(float)
    pmdec_know      = info[9].astype(float)
    pmdece_know     = info[10].astype(float)
    g_mag_know      = info[11].astype(float)
    bp_rp_mag_know  = info[12].astype(float)
    dist_exp        = info[13].astype(float)    
    
    Mgstars = g_mag_know + 5 * np.log10(para_know) - 10 
    
    '''
    ### The following inputs would depends from the data that would be analysis ###
    #From GAIA: this should be  VOTABLE files
    VOTGaia      = parse_single_table(path_inputdata + Clusterdata[z]) #, 'r'
    vodat_gaia   = VOTGaia.array
    IDGaia       = np.array(VOTGaia.array['designation'])
    ID           = VOTGaia.array['source_id']
    ra           = VOTGaia.array['ra']   
    ra_error     = VOTGaia.array['ra_error']
    dec          = VOTGaia.array['dec']   
    dec_error    = VOTGaia.array['dec_error']
    para         = VOTGaia.array['parallax']   
    para_error   = VOTGaia.array['parallax_error']
    pmra         = VOTGaia.array['pmra']   
    pmra_error   = VOTGaia.array['pmra_error']
    pmdec        = VOTGaia.array['pmdec']   
    pmdec_error  = VOTGaia.array['pmdec_error']
    g_mag        = VOTGaia.array['phot_g_mean_mag']   
    bp_rp_mag    = VOTGaia.array['bp_rp']
    bp_g_mag     = VOTGaia.array['bp_g']
    g_rp_mag     = VOTGaia.array['g_rp']
    '''
    
    # From a txt file: This should be an ASCII file 
    data         = scipy.genfromtxt(path_inputdata + Clusterdata[z], dtype='string',comments='#') #
    ID           = data[:,0].astype(float) 
    ra           = data[:,1].astype(float)   
    ra_error     = data[:,2].astype(float)  
    dec          = data[:,3].astype(float) 
    dec_error    = data[:,4].astype(float) 
    para         = data[:,5].astype(float) 
    para_error   = data[:,6].astype(float) 
    pmra         = data[:,7].astype(float) 
    pmra_error   = data[:,8].astype(float) 
    pmdec        = data[:,9].astype(float)  
    pmdec_error  = data[:,10].astype(float) 
    g_mag        = data[:,11].astype(float) 
    bp_rp_mag    = data[:,12].astype(float) 
    bp_g_mag     = data[:,13].astype(float) 
    g_rp_mag     = data[:,14].astype(float) 
    
       
    MG = g_mag + 5 * np.log10(para) - 10                  
      
        
    print '\033[1;34m'
    print("Object name = {0}".format(star_name))
    print("Initial number of objects = {0}".format(len(para)))
    print("Known value of PARALLAX and error = {:.3f} +/- {:.3f}[mas]".format(para_know,parae_know))
    print("Known value of PMRA and error = {:.3f} +/- {:.3f}[mas/yr]".format(pmra_know,pmrae_know))
    print("Known value of PMDEC and error = {:.3f} +/- {:.3f}[mas/yr]".format(pmdec_know,pmdece_know))
    print '\033[1;30m \n'
    
    
    # Creating individual folders for each star
    filename = str(star_name)                          # This would be the name for the new folder to save the data for each star
    results_dir = os.path.join(pathfig, filename)      # This join the original path, to created the new one
    if not os.path.isdir(results_dir):                 # This line and next one would check if the folder have been alreary created or not
        os.makedirs(results_dir)
    
    
    ## Building the Density profile 
    starRA_list   = np.linspace(starRA, starRA, len(ra)) 
    starDEC_list  = np.linspace(starDEC,starDEC, len(dec))
    DeltaRA       = np.cos(dec/(180/np.pi)) * (ra - (starRA_list))
    DeltaDEC      = dec - (starDEC_list)
    Distance      = np.sqrt(np.power(DeltaRA,2) + np.power(DeltaDEC,2))
    
    count, bins, bars = plt.hist(Distance,30, color='white',histtype=u'barstacked', rwidth=0.8, lw=2, alpha=0.8) 
    plt.close()
     
    newbin          = bins - np.min(Distance) 
    dist_lower      = np.linspace(np.min(newbin), np.max(newbin), len(count))
    diff            = np.diff(dist_lower)
    high            = np.linspace(np.min(diff), np.max(diff), len(count))
    dist_high       = dist_lower + high
    ringarea        = np.pi * np.power(dist_high,2) - np.pi * np.power(dist_lower,2)
    density         = count / ringarea
    count_err       = np.sqrt(count)
    density_err     = count_err / ringarea
    
    # Figures
    plt.figure(figsize=(12.5,6.5))   
    
    ## Parallax
    plt.subplot(2,4,1)
    count_para, bins_para, bars_para = plt.hist(para, 40, color='greenyellow',histtype=u'barstacked', rwidth=0.8, lw=2, alpha=0.8)
    plt.axvline(x=para_know, color='k', linestyle='--', lw=1)
    plt.xlabel(r'$\varpi $',fontsize=14)
    plt.ylabel('Number objects',fontsize=11)
    plt.tick_params(axis='both', which='major', labelsize=11)
    plt.subplots_adjust(hspace=0.5,wspace=0.6)

    ## Proper motion RA
    plt.subplot(2,4,2)
    count_pmra, bins_pmra, bars_pmra =plt.hist(pmra, 40, color='aqua',histtype=u'barstacked', rwidth=0.8, lw=2, alpha=0.8)
    plt.axvline(x=pmra_know, color='k', linestyle='--', lw=1)
    plt.xlabel(r'$\mu_{\alpha*} $',fontsize=14)
    plt.tick_params(axis='both', which='major', labelsize=11)
    plt.subplots_adjust(hspace=0.5,wspace=0.6)

    ## Proper motion DEC
    plt.subplot(2,4,3)
    count_pmdec, bins_pmdec, bars_pmdec = plt.hist(pmdec,40, color='royalblue',histtype=u'barstacked', rwidth=0.8, lw=2, alpha=0.8)
    plt.axvline(x=pmdec_know, color='k', linestyle='--', lw=1)
    plt.xlabel(r'$\mu_{\delta} $',fontsize=14)
    plt.tick_params(axis='both', which='major', labelsize=11)
    plt.subplots_adjust(hspace=0.5,wspace=0.6)

    #Density profile
    plt.subplot(2,4,4)
    plt.plot(dist_high,density, '--.', color='midnightblue', alpha=0.8)
    plt.errorbar(dist_high,density, yerr=density_err, fmt='.', color='midnightblue', alpha=0.8, lw=.8)
    plt.xlabel('Distance (deg) ',fontsize=11) 
    plt.ylabel('Density (deg$^{-2}$)',fontsize=11) 
    plt.tick_params(axis='both', which='major', labelsize=11)
    plt.subplots_adjust(hspace=0.5,wspace=0.6)
    
    ## Spatial distribution
    plt.subplot(2,4,5)
    sc=plt.scatter(ra,dec, c=para, marker='o', s=(para*4.5), cmap= 'YlGnBu', norm=None, vmin=np.min(para), vmax=np.max(para), alpha=None, linewidths=None, verts=None, edgecolors=None, hold=None, data=None)
    plt.plot(starRA,starDEC,'*',ms=6,color='darkviolet',alpha=0.65)
    plt.xlabel('RA',fontsize=11)
    plt.ylabel('DEC',fontsize=11)
    plt.gca().invert_xaxis()
    plt.colorbar(sc).remove()
    plt.tick_params(axis='both', which='major', labelsize=11)
    plt.subplots_adjust(wspace=0.5)
    
    ## Proper motion  distribution 
    plt.subplot(2,4,6)
    sc=plt.scatter(pmra,pmdec, c=para, marker='o',s=(para*4.5), cmap= 'YlGnBu', zorder=  35, norm=None, vmin=np.min(para), vmax=np.max(para), alpha=None, linewidths=None, verts=None, edgecolors=None, hold=None, data=None)
    plt.plot(pmra_know,pmdec_know,'*',ms=6, color='darkviolet',alpha=0.65)
    plt.colorbar(sc).remove()
    plt.xlabel(r'$\mu_{\alpha*} $',fontsize=14)
    plt.ylabel(r'$\mu_{\delta} $',fontsize=14)
    plt.tick_params(axis='both', which='major', labelsize=11)
    plt.subplots_adjust(wspace=0.5)
    
    ## CMD
    plt.subplot(2,4,7)
    sc=plt.scatter(bp_rp_mag,MG, c=para, marker='o',s=(para*4.5), cmap= 'YlGnBu', zorder=  35,norm=None, vmin=np.min(para), vmax=np.max(para), alpha=None, linewidths=None, verts=None, edgecolors=None, hold=None, data=None)
    plt.plot(bp_rp_mag_know,Mgstars,'*',ms=6,color='darkviolet',alpha=0.65)
    plt.colorbar(sc, aspect=25).set_label(label=r'$\varpi $', size=12)
    plt.xlabel('G$_{Bp}$ - G$_{Rp}$',fontsize=11)
    plt.ylabel('M$_{G}$',fontsize=11)
    plt.gca().invert_yaxis()
    plt.tick_params(axis='both', which='major', labelsize=11)
    plt.subplots_adjust(wspace=0.5)

    # Information of the target at this stage
    plt.subplot(2,4,8)
    plt.axis([0, 10, 0, 10])
    plt.text(0.5, 9, "Object= {0}".format(star_name), wrap=True, fontsize=10)
    plt.text(0.5, 8, "Number of objects= {0}".format(len(para)), wrap=True, fontsize=10)
    plt.text(0.5, 6, "Known PARALLAX= {:.2f} +/- {:.2f} [mas]".format(para_know,parae_know), wrap=True, fontsize=10)
    plt.text(0.5, 5, "Known PMRA= {:.2f} +/- {:.2f} [mas/yr]".format(pmra_know,pmrae_know), wrap=True, fontsize=10)
    plt.text(0.5, 4, "Known PMDEC= {:.2f} +/- {:.2f} [mas/yr]".format(pmdec_know,pmdece_know), wrap=True, fontsize=10)
    plt.text(0.5, 2, "Mean PARALLAX= {:.2f} [mas]".format(np.mean(para)) , wrap=True, fontsize=10)
    plt.text(0.5, 1, "Mean PMRA= {:.2f} [mas/yr]".format(np.mean(pmra)), wrap=True, fontsize=10)
    plt.text(0.5, 0, "Mean PMDEC= {:.2f} [mas/yr]".format(np.mean(pmdec)), wrap=True, fontsize=10)
    plt.setp(plt.gca(), frame_on=False, xticks=(), yticks=())
    plt.subplots_adjust(wspace=0.6)

    plt.savefig(results_dir + '/InputData_' + str(star_name) + '.pdf')
    plt.tight_layout(h_pad=1.2, w_pad=1.2)
    plt.show()
    plt.close()
    
    ## Started the selection of the data. First parameter: Parallax
    again = 'y'                
    while again == 'y':
        
        IDGaia_l,ID_l,ra_l,rae_l,dec_l,dece_l,para_l,parae_l,pmra_l,pmrae_l,pmdec_l,pmdece_l,gmag_l,BpRp_l,BpG_l, GRp_l=[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]

        print '\033[1;35m'
        print ("Make your selection!\n")
        print '\033[1;30m'
    

        # Define the range to select the data
        para_lower = input("PARALLAX RANGE>=")
        para_upper = input("PARALLAX RANGE<=")

        for i in range(len(ID)):
            if para[i]>=para_lower and para[i]<=para_upper: 
                IDGaia_l.append(IDGaia[i])
                ID_l.append(ID[i])
                ra_l.append(ra[i])
                rae_l.append(ra_error[i])
                dec_l.append(dec[i])
                dece_l.append(dec_error[i])
                para_l.append(para[i])
                parae_l.append(para_error[i])
                pmra_l.append(pmra[i])            
                pmrae_l.append(pmra_error[i])
                pmdec_l.append(pmdec[i]) 
                pmdece_l.append(pmdec_error[i])
                gmag_l.append(g_mag[i])
                BpRp_l.append(bp_rp_mag[i])
                BpG_l.append(bp_g_mag[i])
                GRp_l.append(g_rp_mag[i])

        Mg_l = gmag_l + 5*np.log10(para_l) - 10
        
        mu_para,mu_pmra,mu_pmdec = np.mean(para_l),np.mean(pmra_l),np.mean(pmdec_l)
    
        # General information before the selection
        print '\033[1;34m'
        print("Number of objects in this selection = {0}".format(len(para_l)))   
        print '\033[1;30m'
        
        print '\033[1;34m'
        print("Known value of PARALLAX = {0:.3f}[mas]".format(para_know))
        print("Known value of PMRA =  {0:.3f}[mas/yr]".format(pmra_know))
        print("Known value of PMDEC = {0:.3f}[mas/yr]".format(pmdec_know))
        print '\033[1;30m'
        
        print '\033[1;34m'
        print("Mean value of PARALLAX = {0:.3f} [mas]".format(mu_para))
        print("Mean value of PMRA = {0:.3f} [mas/yr]".format(mu_pmra))
        print("Mean value of PMDEC = {0:.3f} [mas/yr]".format(mu_pmdec))
        print '\033[1;30m'
    
        starRA_list0   = np.linspace(starRA, starRA, len(ra_l)) 
        starDEC_list0  = np.linspace(starDEC,starDEC, len(dec_l))  
        DECresult0     = np.asarray(dec_l)                     
        DeltaRA0       = np.cos(DECresult0/(180/np.pi)) * (ra_l - (starRA_list0))
        DeltaDEC0      = dec_l - (starDEC_list0)
        Distance0      = np.sqrt(np.power(DeltaRA0,2) + np.power(DeltaDEC0,2))
    
        count0, bins0, bars0 = plt.hist(Distance0, 30, color='white',histtype=u'barstacked', rwidth=0.8, lw=2, alpha=0.8) 
        plt.close()
        
        newbin0          = bins0 - np.min(Distance0) 
        dist_lower0      = np.linspace(np.min(newbin0), np.max(newbin0), len(count0))
        diff0            = np.diff(dist_lower0)
        high0            = np.linspace(np.min(diff0), np.max(diff0), len(count0))
        dist_high0       = dist_lower0 + high0
        ringarea0        = np.pi * np.power(dist_high0,2) - np.pi * np.power(dist_lower0,2)
        density0         = count0 / ringarea0
        count_err0       = np.sqrt(count0)
        density_err0     = count_err0 / ringarea0
    
        plt.figure(figsize=(12.5,6.5))    

        ## Parallax 
        plt.subplot(2,4,1)
        count_para0, bins_para0, bars_para0 = plt.hist(para_l, 35,  color='greenyellow',histtype=u'barstacked', rwidth=0.8, lw=2, alpha=0.8)
        plt.axvline(x=para_know, color='k', linestyle='--', lw=1)
        plt.xlabel(r'$\varpi $',fontsize=14)
        plt.ylabel('Number objects',fontsize=11)
        plt.tick_params(axis='both', which='major', labelsize=11)
        plt.subplots_adjust(hspace=0.5,wspace=0.6)

        ## Proper motion RA
        plt.subplot(2,4,2)
        count_pmra0, bins_pmra0, bars_pmra0 = plt.hist(pmra_l,   35,  color='aqua',histtype=u'barstacked', rwidth=0.8, lw=2, alpha=0.8)
        plt.axvline(x=pmra_know, color='k', linestyle='--', lw=1)
        plt.xlabel(r'$\mu_{\alpha*} $',fontsize=14)
        plt.tick_params(axis='both', which='major', labelsize=11)
        plt.subplots_adjust(hspace=0.5,wspace=0.6)

        ## Proper motion DEC
        plt.subplot(2,4,3)
        count_pmdec0, bins_pmdec0, bars_pmdec0 = plt.hist(pmdec_l, 35,  color='royalblue',histtype=u'barstacked', rwidth=0.8, lw=2, alpha=0.8)
        plt.axvline(x=pmdec_know, color='k', linestyle='--', lw=1)
        plt.xlabel(r'$\mu_{\delta} $',fontsize=14)
        plt.tick_params(axis='both', which='major', labelsize=11)
        plt.subplots_adjust(hspace=0.5,wspace=0.6)
        
        ## Density profile
        plt.subplot(2,4,4)
        plt.plot(dist_high0,density0, '--.', color='midnightblue', alpha=0.8)
        plt.errorbar(dist_high0,density0, yerr=density_err0, fmt='.', color='midnightblue', alpha=0.8, lw=.8)
        plt.xlabel('Distance (deg) ',fontsize=11) 
        plt.ylabel('Density (deg$^{-2}$)',fontsize=11) 
        plt.tick_params(axis='both', which='major', labelsize=11)
        plt.subplots_adjust(hspace=0.5,wspace=0.6)
                
        ## Spatial distribution
        plt.subplot(2,4,5)
        sc=plt.scatter(ra_l,dec_l, c=para_l, marker='o', s=(np.array(para_l)*4.5) ,cmap= 'YlGnBu', norm=None, vmin=np.min(para_l), vmax=np.max(para_l), alpha=None, linewidths=None, verts=None, edgecolors=None, hold=None, data=None)
        plt.plot(starRA,starDEC,'*',ms=6,color='darkviolet',alpha=0.65)
        plt.xlabel('RA',fontsize=11)
        plt.ylabel('DEC',fontsize=11)
        plt.gca().invert_xaxis()
        plt.colorbar(sc).remove()
        plt.tick_params(axis='both', which='major', labelsize=11)
        plt.subplots_adjust(wspace=0.5)
  
    
        ## Proper motion  distribution
        plt.subplot(2,4,6)
        sc=plt.scatter(pmra_l,pmdec_l, c=para_l, marker='o',s=(np.array(para_l)*4.5), cmap= 'YlGnBu', zorder=  35, norm=None, vmin=np.min(para_l), vmax=np.max(para_l), alpha=None, linewidths=None, verts=None, edgecolors=None, hold=None, data=None)
        plt.plot(pmra_know,pmdec_know,'*',ms=6, color='darkviolet',alpha=0.65)
        plt.colorbar(sc).remove()
        plt.xlabel(r'$\mu_{\alpha} $',fontsize=14)
        plt.ylabel(r'$\mu_{\delta} $',fontsize=14)
        plt.tick_params(axis='both', which='major', labelsize=11)
        plt.subplots_adjust(wspace=0.5)    

        ## CMD
        plt.subplot(2,4,7)
        sc=plt.scatter(BpRp_l,Mg_l, c=para_l, marker='o',s=(np.array(para_l)*4.5), cmap= 'YlGnBu', zorder=  35,norm=None, vmin=np.min(para_l), vmax=np.max(para_l), alpha=None, linewidths=None, verts=None, edgecolors=None, hold=None, data=None)
        plt.plot(bp_rp_mag_know,Mgstars,'*',ms=6,color='darkviolet',alpha=0.65)
        plt.colorbar(sc, aspect=25).set_label(label=r'$\varpi $', size=10)
        plt.xlabel('G$_{Bp}$ - G$_{Rp}$',fontsize=11)
        plt.ylabel('M$_{G}$',fontsize=11)
        plt.gca().invert_yaxis()
        plt.tick_params(axis='both', which='major', labelsize=12)
        plt.subplots_adjust(wspace=0.5)

        # Information of the target at this stage
        plt.subplot(2,4,8)
        plt.axis([0, 10, 0, 10])
        plt.text(0.5, 9, "Object= {0}".format(star_name), wrap=True, fontsize=10)
        plt.text(0.5, 8, "Number of objects= {0}".format(len(para_l)), wrap=True, fontsize=10)
        plt.text(0.5, 6, "Known PARALLAX= {:.2f} +/- {:.2f} [mas]".format(para_know,parae_know), wrap=True, fontsize=10)
        plt.text(0.5, 5, "Known PMRA= {:.2f} +/- {:.2f} [mas/yr]".format(pmra_know,pmrae_know), wrap=True, fontsize=10)
        plt.text(0.5, 4, "Known PMDEC= {:.2f} +/- {:.2f} [mas/yr]".format(pmdec_know,pmdece_know), wrap=True, fontsize=10)
        plt.text(0.5, 2, "Mean PARALLAX = {:.2f} [mas]".format(mu_para) , wrap=True, fontsize=10)
        plt.text(0.5, 1, "Mean PMRA= {:.2f} [mas/yr]".format(mu_pmra), wrap=True, fontsize=10)
        plt.text(0.5, 0, "Mean PMDEC= {:.2f} [mas/yr]".format(mu_pmdec), wrap=True, fontsize=10)
        plt.setp(plt.gca(), frame_on=False, xticks=(), yticks=())
        plt.subplots_adjust(wspace=0.6)
    
        plt.savefig(results_dir + '/PARALLAX_' + str(star_name) + '[' +str(para_lower)+ '_' + str(para_upper)+ ';' +str(len(para_l))+ ']' + '.pdf')
        plt.tight_layout(h_pad=1.2, w_pad=1.2)
        plt.show()
        plt.close()
    
        
        again=str(input('Choose a new Parallax range? (y/n):')) 
    
    ## Started the selection of another parameter: Proper motion in RA
    again='y'
    while again=='y':
    
        IDGaia1_l,ID1_l,ra1_l,rae1_l,dec1_l,dece1_l,para1_l,parae1_l,pmra1_l,pmrae1_l,pmdec1_l,pmdece1_l,gmag1_l,BpRp1_l,BpG1_l, GRp1_l=[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]

        print '\033[1;35m'
        print ("Make your selection!\n")
        print '\033[1;30m'
    
        # Define the range to select the data
        pmra_lower = input("PMRA RANGE >=")
        pmra_upper = input("PMRA RANGE<=")
        
        
        for i in range(len(ID_l)):
            if pmra_l[i]>=pmra_lower and pmra_l[i]<=pmra_upper: 
                IDGaia1_l.append(IDGaia_l[i])
                ID1_l.append(ID_l[i])
                ra1_l.append(ra_l[i])
                rae1_l.append(rae_l[i])
                dec1_l.append(dec_l[i])
                dece1_l.append(dece_l[i])
                para1_l.append(para_l[i])
                parae1_l.append(parae_l[i])
                pmra1_l.append(pmra_l[i])            
                pmrae1_l.append(pmrae_l[i])
                pmdec1_l.append(pmdec_l[i]) 
                pmdece1_l.append(pmdece_l[i])
                gmag1_l.append(gmag_l[i])
                BpRp1_l.append(BpRp_l[i])
                BpG1_l.append(BpG_l[i])
                GRp1_l.append(GRp_l[i])
        
        Mg1_l = gmag1_l + 5*np.log10(para1_l) - 10

        
        mu_para1,mu_pmra1,mu_pmdec1 = np.mean(para1_l),np.mean(pmra1_l),np.mean(pmdec1_l)

        print '\033[1;34m'
        print("Number of objects in this selection = {0}".format(len(para1_l)))   
        print '\033[1;30m'
        
        print '\033[1;34m'
        print("Known value of PARALLAX = {0:.3f}[mas]".format(para_know))
        print("Known value of PMRA =  {0:.3f}[mas/yr]".format(pmra_know))
        print("Known value of PMDEC =  {0:.3f}[mas/yr]".format(pmdec_know))
        print '\033[1;30m'

        print '\033[1;34m'
        print("Mean value of PARALLAX = {0:.3f} [mas]".format(mu_para1))
        print("Mean value of PMRA = {0:.3f} [mas/yr]".format(mu_pmra1))
        print("Mean value of PMDEC = {0:.3f} [mas/yr]".format(mu_pmdec1))
        print '\033[1;30m'
    
        
        starRA_list1   = np.linspace(starRA, starRA, len(ra1_l)) 
        starDEC_list1  = np.linspace(starDEC,starDEC, len(dec1_l))  
        DECresult1     = np.asarray(dec1_l)
        DeltaRA1       = np.cos(DECresult1/(180/np.pi)) * (ra1_l - (starRA_list1))
        DeltaDEC1      = dec1_l - (starDEC_list1)
        Distance1      = np.sqrt(np.power(DeltaRA1,2) + np.power(DeltaDEC1,2))
    
        count1, bins1, bars1 = plt.hist(Distance1, 30,color='white',histtype=u'barstacked', rwidth=0.8, lw=2, alpha=0.8) 
        plt.close()
        
        newbin1          = bins1 - np.min(Distance1) 
        dist_lower1      = np.linspace(np.min(newbin1), np.max(newbin1), len(count1))
        diff1            = np.diff(dist_lower1)
        high1            = np.linspace(np.min(diff1), np.max(diff1), len(count1))
        dist_high1       = dist_lower1 + high1
        ringarea1        = np.pi * np.power(dist_high1,2) - np.pi * np.power(dist_lower1,2)
        density1         = count1 / ringarea1
        count_err1       = np.sqrt(count1)
        density_err1     = count_err1 / ringarea1
        
        
        plt.figure(figsize=(12.5,6.5))    

        ## Parallax 
        plt.subplot(2,4,1)
        count_para1, bins_para1, bars_para1 =  plt.hist(para1_l, 35,  color='greenyellow',histtype=u'barstacked', rwidth=0.8, lw=2, alpha=0.8)
        plt.axvline(x=para_know, color='k', linestyle='--', lw=1)
        plt.xlabel(r'$\varpi $',fontsize=14)
        plt.ylabel('Number objects',fontsize=11)
        plt.tick_params(axis='both', which='major', labelsize=11)
        plt.subplots_adjust(hspace=0.5,wspace=0.6)

        ## Proper motion RA
        plt.subplot(2,4,2)
        count_pmra1, bins_pmra1, bars_pmra1 = plt.hist(pmra1_l, 35,  color='aqua',histtype=u'barstacked', rwidth=0.8, lw=2, alpha=0.8)
        plt.axvline(x=pmra_know, color='k', linestyle='--', lw=1)
        plt.xlabel(r'$\mu_{\alpha} $',fontsize=14)
        plt.tick_params(axis='both', which='major', labelsize=11)
        plt.subplots_adjust(hspace=0.5,wspace=0.6)

        ## Proper motion DEC
        plt.subplot(2,4,3)
        count_pmdec1, bins_pmdec1, bars_pmdec1 = plt.hist(pmdec1_l, 35,   color='royalblue',histtype=u'barstacked', rwidth=0.8, lw=2, alpha=0.8)
        plt.axvline(x=pmdec_know, color='k', linestyle='--', lw=1)
        plt.xlabel(r'$\mu_{\delta} $',fontsize=14)
        plt.tick_params(axis='both', which='major', labelsize=11)
        plt.subplots_adjust(hspace=0.5,wspace=0.6)
    
        # Density profile
        plt.subplot(2,4,4)
        plt.plot(dist_high1,density1, '--.', color='midnightblue', alpha=0.8)
        plt.errorbar(dist_high1,density1, yerr=density_err1, fmt='.', color='midnightblue', alpha=0.8, lw=.8)
        plt.xlabel('Distance (deg) ',fontsize=11) 
        plt.ylabel('Density (deg$^{-2}$)',fontsize=11) 
        plt.tick_params(axis='both', which='major', labelsize=11)
        plt.subplots_adjust(hspace=0.5,wspace=0.6)


        ## Spatial distribution
        plt.subplot(2,4,5)
        sc=plt.scatter(ra1_l,dec1_l, c=para1_l, marker='o',s=(np.array(para1_l)*4.5), cmap= 'YlGnBu', norm=None, vmin=np.min(para1_l), vmax=np.max(para1_l), alpha=None, linewidths=None, verts=None, edgecolors=None, hold=None, data=None)
        plt.plot(starRA,starDEC,'*',ms=6,color='darkviolet',alpha=0.65)
        plt.xlabel('RA',fontsize=11)
        plt.ylabel('DEC',fontsize=11)
        plt.gca().invert_xaxis()
        plt.colorbar(sc).remove()
        plt.tick_params(axis='both', which='major', labelsize=11)
        plt.subplots_adjust(wspace=0.5)


        ## Proper motion  distribution
        plt.subplot(2,4,6)
        sc=plt.scatter(pmra1_l,pmdec1_l, c=para1_l,s=(np.array(para1_l)*4.5), marker='o',cmap= 'YlGnBu', zorder=  35, norm=None, vmin=np.min(para1_l), vmax=np.max(para1_l), alpha=None, linewidths=None, verts=None, edgecolors=None, hold=None, data=None)
        plt.plot(pmra_know,pmdec_know,'*',ms=6, color='darkviolet',alpha=0.65)
        plt.colorbar(sc).remove()
        plt.xlabel(r'$\mu_{\alpha*} $',fontsize=14)
        plt.ylabel(r'$\mu_{\delta} $',fontsize=14)
        plt.tick_params(axis='both', which='major', labelsize=11)
        plt.subplots_adjust(wspace=0.5)

        ## CMD
        plt.subplot(2,4,7)
        sc=plt.scatter(BpRp1_l,Mg1_l, c=para1_l, marker='o',s=(np.array(para1_l)*4.5),cmap= 'YlGnBu', zorder=  35,norm=None, vmin=np.min(para1_l), vmax=np.max(para1_l), alpha=None, linewidths=None, verts=None, edgecolors=None, hold=None, data=None)
        plt.plot(bp_rp_mag_know,Mgstars,'*',ms=6,color='darkviolet',alpha=0.65)
        plt.colorbar(sc, aspect=25).set_label(label=r'$\varpi $', size=10)
        plt.xlabel('G$_{Bp}$ - G$_{Rp}$',fontsize=11)
        plt.ylabel('M$_{G}$',fontsize=11)
        plt.gca().invert_yaxis()
        plt.tick_params(axis='both', which='major', labelsize=11)
        plt.subplots_adjust(wspace=0.5)

        # Information of the target at this stage
        plt.subplot(2,4,8)
        plt.axis([0, 10, 0, 10])
        plt.text(0.5, 9, "Object= {0}".format(star_name), wrap=True, fontsize=10)
        plt.text(0.5, 8, "Number of objects= {0}".format(len(para1_l)), wrap=True, fontsize=10)
        plt.text(0.5, 6, "Known PARALLAX= {:.2f} +/- {:.2f} [mas]".format(para_know,parae_know), wrap=True, fontsize=10)
        plt.text(0.5, 5, "Known PMRA= {:.2f} +/- {:.2f} [mas/yr]".format(pmra_know,pmrae_know), wrap=True, fontsize=10)
        plt.text(0.5, 4, "Known PMDEC= {:.2f} +/- {:.2f} [mas/yr]".format(pmdec_know,pmdece_know), wrap=True, fontsize=10)
        plt.text(0.5, 2, "Mean PARALLAX = {:.2f} [mas]".format(mu_para1) , wrap=True, fontsize=10)
        plt.text(0.5, 1, "Mean PMRA= {:.2f} [mas/yr]".format(mu_pmra1), wrap=True, fontsize=10)
        plt.text(0.5, 0, "Mean PMDEC= {:.2f} [mas/yr]".format(mu_pmdec1), wrap=True, fontsize=10)
        plt.setp(plt.gca(), frame_on=False, xticks=(), yticks=())
        plt.subplots_adjust(wspace=0.6)
    
        plt.savefig(results_dir + '/PMRA_' + str(star_name) + '[' +str(pmra_lower)+ '_' + str(pmra_upper)+ ';' +str(len(pmra1_l))+ ']' + '.pdf')
        plt.tight_layout(h_pad=1.2, w_pad=1.2)
        plt.show()
        plt.close()
    
        
        again=str(input('Choose a new Proper Motion in RA range? (y/n):')) 
        
    ## Started the selection of another parameter: Proper motion in DEC
    again='y'
    while again=='y':
        
        IDGaia2_l,ID2_l,ra2_l,rae2_l,dec2_l,dece2_l,para2_l,parae2_l,pmra2_l,pmrae2_l,pmdec2_l,pmdece2_l,gmag2_l,BpRp2_l,BpG2_l, GRp2_l=[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]

        print '\033[1;35m'
        print ("Make your selection!\n")
        print '\033[1;30m'
        
        # Define the range to select the data
        pmdec_lower = input("PMDEC RANGE>=")
        pmdec_upper = input("PMDEC RANGE<=")

        
        for i in range(len(ID1_l)):
            if pmdec1_l[i]>=pmdec_lower and pmdec1_l[i]<=pmdec_upper: 
                IDGaia2_l.append(IDGaia1_l[i])
                ID2_l.append(ID1_l[i])
                ra2_l.append(ra1_l[i])
                rae2_l.append(rae1_l[i])
                dec2_l.append(dec1_l[i])
                dece2_l.append(dece1_l[i])
                para2_l.append(para1_l[i])
                parae2_l.append(parae1_l[i])
                pmra2_l.append(pmra1_l[i])            
                pmrae2_l.append(pmrae1_l[i])
                pmdec2_l.append(pmdec1_l[i]) 
                pmdece2_l.append(pmdece1_l[i])
                gmag2_l.append(gmag1_l[i])
                BpRp2_l.append(BpRp1_l[i])
                BpG2_l.append(BpG1_l[i])
                GRp2_l.append(GRp1_l[i])

        
        Mg2_l = gmag2_l + 5*np.log10(para2_l) - 10

        
        mu_para2,mu_pmra2,mu_pmdec2 = np.mean(para2_l),np.mean(pmra2_l),np.mean(pmdec2_l)

        print '\033[1;34m'
        print("Number of objects in this selection = {0}".format(len(para2_l)))
        print '\033[1;30m'
        
        print '\033[1;34m'
        print("Known value of PARALLAX = {0:.3f}[mas]".format(para_know))
        print("Known value of PMRA = {0:.3f}[mas/yr]".format(pmra_know))
        print("Known value of PMDEC = {0:.3f}[mas/yr]".format(pmdec_know))
        print '\033[1;30m'

        print '\033[1;34m'
        print("Mean value of PARALLAX = {0:.3f} [mas]".format(mu_para2))
        print("Mean value of PMRA = {0:.3f} [mas/yr]".format(mu_pmra2))
        print("Mean value of PMDEC = {0:.3f} [mas/yr]".format(mu_pmdec2))
        print '\033[1;30m'
    
    
        starRA_list2   = np.linspace(starRA, starRA, len(ra2_l)) 
        starDEC_list2  = np.linspace(starDEC,starDEC, len(dec2_l))  
        DECresult2     = np.asarray(dec2_l)
        DeltaRA2       = np.cos(DECresult2/(180/np.pi)) * (ra2_l - (starRA_list2))
        DeltaDEC2      = dec2_l - (starDEC_list2)
        Distance2      = np.sqrt(np.power(DeltaRA2,2) + np.power(DeltaDEC2,2))

        count2, bins2, bars2 = plt.hist(Distance2, 30, color='white',histtype=u'barstacked', rwidth=0.8, lw=2, alpha=0.8) 
        plt.close()
        
        newbin2          = bins2 - np.min(Distance2) 
        dist_lower2      = np.linspace(np.min(newbin2), np.max(newbin2), len(count2))
        diff2            = np.diff(dist_lower2)
        high2            = np.linspace(np.min(diff2), np.max(diff2), len(count2))
        dist_high2       = dist_lower2 + high2
        ringarea2        = np.pi * np.power(dist_high2,2) - np.pi * np.power(dist_lower2,2)
        density2         = count2 / ringarea2
        count_err2       = np.sqrt(count2)
        density_err2     = count_err2 / ringarea2
    
    
        plt.figure(figsize=(12.5,6.5))    

        ## Parallax 
        plt.subplot(2,4,1)
        count_para2, bins_para2, bars_para2 =  plt.hist(para2_l, 35,  color='greenyellow',histtype=u'barstacked', rwidth=0.8, lw=2, alpha=0.8)
        plt.axvline(x=para_know, color='k', linestyle='--', lw=1)
        plt.xlabel(r'$\varpi $',fontsize=14)
        plt.ylabel('Number objects',fontsize=11)
        plt.tick_params(axis='both', which='major', labelsize=11)
        plt.subplots_adjust(hspace=0.5,wspace=0.6)

        ## Proper motion RA
        plt.subplot(2,4,2)
        count_pmra2, bins_pmra2, bars_pmra2 =  plt.hist(pmra2_l, 35,  color='aqua',histtype=u'barstacked', rwidth=0.8, lw=2, alpha=0.8)
        plt.axvline(x=pmra_know, color='k', linestyle='--', lw=1)
        plt.xlabel(r'$\mu_{\alpha*} $',fontsize=14)
        plt.tick_params(axis='both', which='major', labelsize=11)
        plt.subplots_adjust(hspace=0.5,wspace=0.6)

        ## Proper motion DEC
        plt.subplot(2,4,3)
        count_pmdec2, bins_pmdec2, bars_pmdec2 = plt.hist(pmdec2_l, 35,   color='royalblue',histtype=u'barstacked', rwidth=0.8, lw=2, alpha=0.8)
        plt.axvline(x=pmdec_know, color='k', linestyle='--', lw=1)
        plt.xlabel(r'$\mu_{\delta} $',fontsize=14)
        plt.tick_params(axis='both', which='major', labelsize=11)
        plt.subplots_adjust(hspace=0.5,wspace=0.6)
    
        # Density profile
        plt.subplot(2,4,4)
        plt.plot(dist_high2,density2, '--.', color='midnightblue', alpha=0.8)
        plt.errorbar(dist_high2,density2, yerr=density_err2, fmt='.', color='midnightblue', alpha=0.8, lw=.8)
        plt.xlabel('Distance (deg) ',fontsize=11) 
        plt.ylabel('Density (deg$^{-2}$)',fontsize=11) 
        plt.tick_params(axis='both', which='major', labelsize=11)
        plt.subplots_adjust(hspace=0.5, wspace=0.6)


        ## Spatial distribution
        plt.subplot(2,4,5)
        sc=plt.scatter(ra2_l,dec2_l, c=para2_l, marker='o',s=(np.array(para2_l)*4.5), cmap= 'YlGnBu', norm=None, vmin=np.min(para2_l), vmax=np.max(para2_l), alpha=None, linewidths=None, verts=None, edgecolors=None, hold=None, data=None)
        plt.plot(starRA,starDEC,'*',ms=6,color='darkviolet',alpha=0.65)
        plt.xlabel('RA',fontsize=11)
        plt.ylabel('DEC',fontsize=11)
        plt.gca().invert_xaxis()
        plt.colorbar(sc).remove()
        plt.tick_params(axis='both', which='major', labelsize=11)
        plt.subplots_adjust(wspace=0.5)

        ## Proper motion  distribution
        plt.subplot(2,4,6)
        sc=plt.scatter(pmra2_l,pmdec2_l, c=para2_l, marker='o',s=(np.array(para2_l)*4.5),cmap= 'YlGnBu', zorder=  35, norm=None, vmin=np.min(para2_l), vmax=np.max(para2_l), alpha=None, linewidths=None, verts=None, edgecolors=None, hold=None, data=None)
        plt.plot(pmra_know,pmdec_know,'*',ms=6, color='darkviolet',alpha=0.65)
        plt.colorbar(sc).remove()
        plt.xlabel(r'$\mu_{\alpha*} $',fontsize=14)
        plt.ylabel(r'$\mu_{\delta} $',fontsize=14)
        plt.tick_params(axis='both', which='major', labelsize=11)
        plt.subplots_adjust(wspace=0.5)

        ## CMD
        plt.subplot(2,4,7)
        sc=plt.scatter(BpRp2_l,Mg2_l, c=para2_l, marker='o',s=(np.array(para2_l)*4.5),cmap= 'YlGnBu', zorder=  35,norm=None, vmin=np.min(para2_l), vmax=np.max(para2_l), alpha=None, linewidths=None, verts=None, edgecolors=None, hold=None, data=None)
        plt.plot(bp_rp_mag_know,Mgstars,'*',ms=6,color='darkviolet',alpha=0.65)
        plt.colorbar(sc, aspect=25).set_label(label=r'$\varpi $', size=10)
        plt.xlabel('G$_{Bp}$ - G$_{Rp}$',fontsize=11)
        plt.ylabel('M$_{G}$',fontsize=11)
        plt.gca().invert_yaxis()
        plt.tick_params(axis='both', which='major', labelsize=11)
        plt.subplots_adjust(wspace=0.5)

        # Information of the target at this stage
        plt.subplot(2,4,8)
        plt.axis([0, 10, 0, 10])
        plt.text(0.5, 9, "Object= {0}".format(star_name), wrap=True, fontsize=10)
        plt.text(0.5, 8, "Number of objects= {0}".format(len(para2_l)), wrap=True, fontsize=10)
        plt.text(0.5, 6, "Known PARALLAX= {:.2f} +/- {:.2f} [mas]".format(para_know,parae_know), wrap=True, fontsize=10)
        plt.text(0.5, 5, "Known PMRA= {:.2f} +/- {:.2f} [mas/yr]".format(pmra_know,pmrae_know), wrap=True, fontsize=10)
        plt.text(0.5, 4, "Known PMDEC= {:.2f} +/- {:.2f} [mas/yr]".format(pmdec_know,pmdece_know), wrap=True, fontsize=10)
        plt.text(0.5, 2, "Mean PARALLAX = {:.2f} [mas]".format(mu_para2) , wrap=True, fontsize=10)
        plt.text(0.5, 1, "Mean PMRA= {:.2f} [mas/yr]".format(mu_pmra2), wrap=True, fontsize=10)
        plt.text(0.5, 0, "Mean PMDEC= {:.2f} [mas/yr]".format(mu_pmdec2), wrap=True, fontsize=10)
        plt.setp(plt.gca(), frame_on=False, xticks=(), yticks=())
        plt.subplots_adjust(wspace=0.6)

        plt.savefig(results_dir + '/PMDEC_' + str(star_name) + '[' +str(pmdec_lower)+ '_' + str(pmdec_upper)+ ';' + str(len(pmdec2_l))+ ']' + '.pdf')
        plt.tight_layout(h_pad=1.2, w_pad=1.2)
        plt.show()
        plt.close()
    
       
        again=str(input('Choose a new Proper Motion in DEC range? (y/n):')) 
       
    # Up to this point the code would finish, if you select 'y' or continue if 
    # you decide to make more selection over the paramate (Parallax, PMRA or PMDEC)
    # in that order.
    loop=str(input('Would you like to start again the selection? (y/n): '))

    while loop=='n':
        print('Well done! This are your results =D \n')

        print '\033[1;34m'
        print("Number of objects = {0}".format(len(pmra2_l)))
        print '\033[1;30m'
        
        print '\033[1;34m'
        print("Known value of PARALLAX =  {0:.3f}[mas]".format(para_know))
        print("Known value of PMRA =  {0:.3f}[mas/yr]".format(pmra_know))
        print("Known value of PMDEC =  {0:.3f}[mas/yr]".format(pmdec_know))
        print '\033[1;30m'
        
        print '\033[1;34m'
        print("Mean value of PARALLAX = {0:.3f} [mas]".format(mu_para2))
        print("Mean value of PMRA = {0:.3f} [mas/yr]".format(mu_pmra2))
        print("Mean value of PMDEC = {0:.3f} [mas/yr]".format(mu_pmdec2))
        print '\033[1;30m'

        # This would be the output file from CEREAL for the star is running 
        # at the moment.
        file= open(pathdata + '/OutputData_' + str(star_name) +  '.txt','w')
        file.write('#ID\t' 'RA\t' 'RA_Error\t' 'DEC\t' 'DEC_Error\t' 'PARALLAX\t' 'PARALLAX_Error\t' 'PMRA\t' 'PMRA_Error\t' ' PMDEC\t' ' PMDEC_Error\t' ' G_mag\t'  'Bp-Rp_mag\t' 'Bp-G_mag\t' 'G-Rp_mag\n' )
        for i in range(len(ID2_l)):
            file.write("%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n" % (ID2_l[i], ra2_l[i], rae2_l[i], dec2_l[i], dece2_l[i], para2_l[i], parae2_l[i], pmra2_l[i], pmrae2_l[i], pmdec2_l[i], pmdece2_l[i], gmag2_l[i], BpRp2_l[i], BpG2_l[i], GRp2_l[i]))
        file.close()           
        

 
        ## FINAL PLOT WITH THE RESULTS OF THE LAST SELECTION
        ## Estimate de error trought the weighted mean ecuation
        
        #Parallax#
        wpara = 1/np.power(parae2_l,2)
        Sigma_wparasqueare = np.sqrt(1/wpara)
        Xpara = (wpara*para2_l)/(wpara)
        Sigma_wparaT = np.sqrt(1/np.sum(wpara))
        XparaT = np.sum(wpara*para2_l)/np.sum(wpara)

        #PmRA
        wpmra = 1/np.power(pmrae2_l,2)
        Sigma_wpmrasqueare = np.sqrt(1/wpmra)
        Xpmra = (wpmra*pmra2_l)/(wpmra)
        Sigma_wpmraT = np.sqrt(1/np.sum(wpmra))
        XpmraT = np.sum(wpmra*pmra2_l)/np.sum(wpmra)

        #PmDEC
        wpmdec = 1/np.power(pmdece2_l,2)
        Sigma_wpmdecsqueare = np.sqrt(1/wpmdec)
        Xpmdec = (wpmdec*pmdec2_l)/(wpmdec)
        Sigma_wpmdecT = np.sqrt(1/np.sum(wpmdec))
        XpmdecT = np.sum(wpmdec*pmdec2_l)/np.sum(wpmdec)
        
        plt.figure(figsize=(12.5,6.5))
        
        plt.subplot(2,4,1)
        count_para2, bins_para2, bars_para2 =  plt.hist(para2_l, 35,  color='greenyellow',histtype=u'barstacked', rwidth=0.8, lw=2, alpha=0.8)
        plt.axvline(x=para_know, color='k', linestyle='--', lw=1)
        plt.xlabel(r'$\varpi $',fontsize=14)
        plt.ylabel('Number objects',fontsize=11)
        plt.tick_params(axis='both', which='major', labelsize=11)
        plt.subplots_adjust(hspace=0.5,wspace=0.6)

        ## Proper motion RA
        plt.subplot(2,4,2)
        count_pmra2, bins_pmra2, bars_pmra2 =  plt.hist(pmra2_l,   35,  color='aqua',histtype=u'barstacked', rwidth=0.8, lw=2, alpha=0.8)
        plt.axvline(x=pmra_know, color='k', linestyle='--', lw=1)
        plt.xlabel(r'$\mu_{\alpha*} $',fontsize=14)
        plt.tick_params(axis='both', which='major', labelsize=11)
        plt.subplots_adjust(hspace=0.5,wspace=0.6)

        ## Proper motion DEC
        plt.subplot(2,4,3)
        count_pmdec2, bins_pmdec2, bars_pmdec2 = plt.hist(pmdec2_l,  35,   color='royalblue',histtype=u'barstacked', rwidth=0.8, lw=2, alpha=0.8)
        plt.axvline(x=pmdec_know, color='k', linestyle='--', lw=1)
        plt.xlabel(r'$\mu_{\delta} $',fontsize=14)
        plt.tick_params(axis='both', which='major', labelsize=11)
        plt.subplots_adjust(hspace=0.5,wspace=0.6)        
        
        #Density profile
        plt.subplot(2,4,4)
        plt.plot(dist_high2,density2, '--.',color='midnightblue', alpha=0.8)
        plt.errorbar(dist_high2,density2, yerr=density_err2, fmt='.', color='midnightblue', alpha=0.8, lw=.8)
        plt.xlabel('Distance (deg) ',fontsize=11) 
        plt.ylabel('Density (deg$^{-2}$)',fontsize=11) 
        plt.tick_params(axis='both', which='major', labelsize=10)
        plt.subplots_adjust(hspace=0.5,wspace=0.6)
        
         
        ## Spatial distribution
        plt.subplot(2,4,5)
        sc=plt.scatter(ra2_l,dec2_l, c=para2_l, marker='o', s=(np.array(para2_l)*4.5),cmap= 'YlGnBu', norm=None, vmin=np.min(para2_l), vmax=np.max(para2_l), alpha=None, linewidths=None, verts=None, edgecolors=None, hold=None, data=None)
        plt.plot(starRA,starDEC,'*',ms=6,color='darkviolet',alpha=0.65)
        plt.xlabel('RA',fontsize=11)
        plt.ylabel('DEC',fontsize=11)
        plt.gca().invert_xaxis()
        plt.colorbar(sc).remove()
        plt.tick_params(axis='both', which='major', labelsize=11)
        plt.subplots_adjust(hspace=0.5, wspace=0.6)

        ## Proper motion  distribution
        plt.subplot(2,4,6)
        sc=plt.scatter(pmra2_l,pmdec2_l, c=para2_l, marker='o', s=(np.array(para2_l)*4.5),cmap= 'YlGnBu', zorder=  35, norm=None, vmin=np.min(para2_l), vmax=np.max(para2_l), linewidths=None, verts=None, edgecolors=None, hold=None, data=None)
        plt.plot(pmra_know,pmdec_know,'s',ms=5, color='darkviolet',label='Known Value')
        plt.plot(XpmraT,XpmdecT,'o',ms=5,color='darkorange',label='This programme')
        plt.xlabel(r'$\mu_{\alpha*} $',fontsize=14)
        plt.ylabel(r'$\mu_{\delta} $',fontsize=14)
        plt.colorbar(sc).remove()
        plt.tick_params(axis='both', which='major', labelsize=11)
        plt.subplots_adjust(hspace=0.5, wspace=0.6)

        ## CMD
        plt.subplot(2,4,7)
        sc=plt.scatter(BpRp2_l,Mg2_l, c=para2_l, marker='o',s=(np.array(para2_l)*4.5),cmap= 'YlGnBu', zorder=  35,norm=None, vmin=np.min(para2_l), vmax=np.max(para2_l), alpha=None, linewidths=None, verts=None, edgecolors=None, hold=None, data=None)
        plt.plot(bp_rp_mag_know,Mgstars,'*',ms=6,color='darkviolet',alpha=0.65)
        plt.colorbar(sc, aspect=25).set_label(label=r'$\varpi $', size=10)
        plt.xlabel('G$_{Bp}$ - G$_{Rp}$',fontsize=11)
        plt.ylabel('M$_{G}$',fontsize=11)
        plt.gca().invert_yaxis()
        plt.tick_params(axis='both', which='major', labelsize=11)
        plt.subplots_adjust(hspace=0.5, wspace=0.6)

        ## Text of the last result
        plt.subplot(2,4,8)
        plt.axis([0, 10, 0, 10])
        plt.text(0.5, 9, "Object= {0}".format(star_name), wrap=True, fontsize=10)
        plt.text(0.5, 8, "Number of objects= {0}".format(len(para2_l)), wrap=True, fontsize=10)
        plt.text(0.5, 7, "Known PARALLAX= {:.2f} +/- {:.2f} [mas]".format(para_know,parae_know), wrap=True, fontsize=10)
        plt.text(0.5, 6, "Known PMRA= {:.2f} +/- {:.2f} [mas/yr]".format(pmra_know,pmrae_know), wrap=True, fontsize=10)
        plt.text(0.5, 5, "Known PMDEC= {:.2f} +/- {:.2f} [mas/yr]".format(pmdec_know,pmdece_know), wrap=True, fontsize=10)
        plt.text(0.5, 4, "Mean PARALLAX = {:.2f} +/- {:.2f}  [mas]".format(XparaT,Sigma_wparaT) , wrap=True, fontsize=10)
        plt.text(0.5, 3, "Mean PMRA= {:.2f} +/- {:.2f}  [mas/yr]".format(XpmraT,Sigma_wpmraT), wrap=True, fontsize=10)
        plt.text(0.5, 2, "Mean PMDEC= {:.2f} +/- {:.2f}  [mas/yr]".format(XpmdecT,Sigma_wpmdecT), wrap=True, fontsize=10)
        plt.text(0.5, 0, "Square-> Known Value; Circle-> CEREAL", wrap=True, fontsize=8)
        plt.setp(plt.gca(), frame_on=False, xticks=(), yticks=())
        plt.subplots_adjust(wspace=0.6)

        plt.savefig(results_dir + '/OutputData_' + str(star_name) + '.pdf')
        plt.tight_layout(h_pad=1.2, w_pad=1.2)
        plt.show()
        plt.close()

        break
        
    #New loop after the started again    
    while loop=='y':
        
        again=str(input('Introduce a new Parallax range (y): '))
        while again=='y':
            
            IDGaia3_l,ID3_l,ra3_l,rae3_l,dec3_l,dece3_l,para3_l,parae3_l,pmra3_l,pmrae3_l,pmdec3_l,pmdece3_l,gmag3_l,BpRp3_l,BpG3_l,GRp3_l = [],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]
 
            print '\033[1;35m'
            print ("Make your selection!\n")
            print '\033[1;30m'
    
            # Define the range to select the data
            para_lower = input("PARALLAX RANGE>=")
            para_upper = input("PARALLAX RANGE<=")

            
            for i in range(len(ID2_l)):
                if para2_l[i]>=para_lower and para2_l[i]<=para_upper: 
                    IDGaia3_l.append(IDGaia2_l[i])
                    ID3_l.append(ID2_l[i])
                    ra3_l.append(ra2_l[i])
                    rae3_l.append(rae2_l[i])
                    dec3_l.append(dec2_l[i])
                    dece3_l.append(dec2_l[i])
                    para3_l.append(para2_l[i])
                    parae3_l.append(parae2_l[i])
                    pmra3_l.append(pmra2_l[i])            
                    pmrae3_l.append(pmrae2_l[i])
                    pmdec3_l.append(pmdec2_l[i]) 
                    pmdece3_l.append(pmdece2_l[i])
                    gmag3_l.append(gmag2_l[i])
                    BpRp3_l.append(BpRp2_l[i])
                    BpG3_l.append(BpG2_l[i])
                    GRp3_l.append(GRp2_l[i])

            
            Mg3_l = gmag3_l + 5*np.log10(para3_l) - 10
            
            
            mu_para3,mu_pmra3,mu_pmdec3 = np.mean(para3_l),np.mean(pmra3_l),np.mean(pmdec3_l)

            print '\033[1;34m'
            print("Number of objects in this selection = {0}".format(len(para3_l)))
            print '\033[1;30m'
             
            print '\033[1;34m'
            print("Known value of PARALLAX = {0:.3f}[mas]".format(para_know))
            print("Known value of PMRA = {0:.3f}[mas/yr]".format(pmra_know))
            print("Known value of PMDEC =  {0:.3f}[mas/yr]".format(pmdec_know))
            print '\033[1;30m'
            
            print '\033[1;34m'
            print("Mean value of PARALLAX = {0:.3f} [mas]".format(mu_para3))
            print("Mean value of PMRA = {0:.3f} [mas/yr]".format(mu_pmra3))
            print("Mean value of PMDEC = {0:.3f} [mas/yr]".format(mu_pmdec3))
            print '\033[1;30m'


            
            starRA_list3   = np.linspace(starRA, starRA, len(ra3_l)) 
            starDEC_list3  = np.linspace(starDEC,starDEC, len(dec3_l))  
            DECresult3     = np.asarray(dec3_l)
            DeltaRA3       = np.cos(DECresult3/(180/np.pi)) * (ra3_l - (starRA_list3))
            DeltaDEC3      = dec3_l - (starDEC_list3)
            Distance3      = np.sqrt(np.power(DeltaRA3,2) + np.power(DeltaDEC3,2))

            count3, bins3, bars3 = plt.hist(Distance3, 30, color='white',histtype=u'barstacked', rwidth=0.8, lw=2, alpha=0.8) 
            plt.close()
            
            newbin3          = bins3 - np.min(Distance3) 
            dist_lower3      = np.linspace(np.min(newbin3), np.max(newbin3), len(count3))
            diff3            = np.diff(dist_lower3)
            high3            = np.linspace(np.min(diff3), np.max(diff3), len(count3))
            dist_high3       = dist_lower3 + high3
            ringarea3        = np.pi * np.power(dist_high3,2) - np.pi * np.power(dist_lower3,2)
            density3         = count3 / ringarea3
            count_err3       = np.sqrt(count3)
            density_err3     = count_err3 / ringarea3
            

            plt.figure(figsize=(12.5,6.5))

            ## Parallax 
            plt.subplot(2,4,1)
            count_para3, bins_para3, bars_para3 = plt.hist(para3_l,  35,  color='greenyellow',histtype=u'barstacked', rwidth=0.8, lw=2, alpha=0.8)
            plt.axvline(x=para_know, color='k', linestyle='--', lw=1)
            plt.xlabel(r'$\varpi $',fontsize=14)
            plt.ylabel('Number objects',fontsize=11)
            plt.tick_params(axis='both', which='major', labelsize=11)
            plt.subplots_adjust(hspace=0.5,wspace=0.6)

            ## Proper motion RA
            plt.subplot(2,4,2)
            count_pmra3, bins_pmra3, bars_pmra3 = plt.hist(pmra3_l,  35,  color='aqua',histtype=u'barstacked', rwidth=0.8, lw=2, alpha=0.8)
            plt.axvline(x=pmra_know, color='k', linestyle='--', lw=1)
            plt.xlabel(r'$\mu_{\alpha*} $',fontsize=14)
            plt.tick_params(axis='both', which='major', labelsize=11)
            plt.subplots_adjust(hspace=0.5,wspace=0.6)

            ## Proper motion DEC
            plt.subplot(2,4,3)
            count_pmdec3, bins_pmdec3, bars_pmdec3 = plt.hist(pmdec3_l, 35,  color='royalblue',histtype=u'barstacked', rwidth=0.8, lw=2, alpha=0.8)           
            plt.axvline(x=pmdec_know, color='k', linestyle='--', lw=1)
            plt.xlabel(r'$\mu_{\delta} $',fontsize=14)
            plt.tick_params(axis='both', which='major', labelsize=11)
            plt.subplots_adjust(hspace=0.5,wspace=0.6)
            
            ## Density profile
            plt.subplot(2,4,4)
            plt.plot(dist_high3,density3, '--.', color='midnightblue', alpha=0.8)
            plt.errorbar(dist_high3,density3, yerr=density_err3, fmt='.', color='midnightblue', alpha=0.8, lw=.8)
            plt.xlabel('Distance (deg) ',fontsize=11)
            plt.ylabel('Density (deg$^{-2}$)',fontsize=11)
            plt.tick_params(axis='both', which='major', labelsize=11)
            plt.subplots_adjust(hspace=0.5,wspace=0.6)
            
            ## Spatial distribution
            plt.subplot(2,4,5)
            sc=plt.scatter(ra3_l,dec3_l, c=para3_l, marker='o', s=(np.array(para3_l)*4.5),cmap= 'YlGnBu', norm=None, vmin=np.min(para3_l), vmax=np.max(para3_l), alpha=None, linewidths=None, verts=None, edgecolors=None, hold=None, data=None)
            plt.plot(starRA,starDEC,'*',ms=6,color='darkviolet',alpha=0.65)
            plt.xlabel('RA',fontsize=11)
            plt.ylabel('DEC',fontsize=11)
            plt.gca().invert_xaxis()
            plt.colorbar(sc).remove()
            plt.tick_params(axis='both', which='major', labelsize=11)
            plt.subplots_adjust(wspace=0.5)

            ## Proper motion  distribution
            plt.subplot(2,4,6)
            sc=plt.scatter(pmra3_l,pmdec3_l, c=para3_l, marker='o',s=(np.array(para3_l)*4.5),cmap= 'YlGnBu', zorder=  35, norm=None, vmin=np.min(para3_l), vmax=np.max(para3_l), alpha=None, linewidths=None, verts=None, edgecolors=None, hold=None, data=None)
            plt.plot(pmra_know,pmdec_know,'*',ms=6, color='darkviolet',alpha=0.65)
            plt.colorbar(sc).remove()
            plt.xlabel(r'$\mu_{\alpha*} $',fontsize=14)
            plt.ylabel(r'$\mu_{\delta} $',fontsize=14)
            plt.tick_params(axis='both', which='major', labelsize=10)
            plt.subplots_adjust(wspace=0.5)

            ## CMD
            plt.subplot(2,4,7)
            sc=plt.scatter(BpRp3_l,Mg3_l, c=para3_l, marker='o',s=(np.array(para3_l)*4.5),cmap= 'YlGnBu', zorder=  35,norm=None, vmin=np.min(para3_l), vmax=np.max(para3_l), alpha=None, linewidths=None, verts=None, edgecolors=None, hold=None, data=None)
            plt.plot(bp_rp_mag_know,Mgstars,'*',ms=6,color='darkviolet',alpha=0.65)
            plt.colorbar(sc, aspect=25).set_label(label=r'$\varpi $', size=10)
            plt.xlabel('G$_{Bp}$ - G$_{Rp}$',fontsize=11)
            plt.ylabel('M$_{G}$',fontsize=11)
            plt.gca().invert_yaxis()
            plt.tick_params(axis='both', which='major', labelsize=11)
            plt.subplots_adjust(wspace=0.5)

            # Information of the target at this stage
            plt.subplot(2,4,8)
            plt.axis([0, 10, 0, 10])
            plt.text(0.5, 9, "Object= {0}".format(star_name), wrap=True, fontsize=10)
            plt.text(0.5, 8, "Number of objects= {0}".format(len(para3_l)), wrap=True, fontsize=10)
            plt.text(0.5, 6, "Known PARALLAX= {:.2f} +/- {:.2f} [mas]".format(para_know,parae_know), wrap=True, fontsize=10)
            plt.text(0.5, 5, "Known PMRA= {:.2f} +/- {:.2f} [mas/yr]".format(pmra_know,pmrae_know), wrap=True, fontsize=10)
            plt.text(0.5, 4, "Known PMDEC= {:.2f} +/- {:.2f} [mas/yr]".format(pmdec_know,pmdece_know), wrap=True, fontsize=10)
            plt.text(0.5, 2, "Mean PARALLAX = {:.2f} [mas]".format(mu_para3) , wrap=True, fontsize=10)
            plt.text(0.5, 1, "Mean PMRA= {:.2f} [mas/yr]".format(mu_pmra3), wrap=True, fontsize=10)
            plt.text(0.5, 0, "Mean PMDEC= {:.2f} [mas/yr]".format(mu_pmdec3), wrap=True, fontsize=10)
            plt.setp(plt.gca(), frame_on=False, xticks=(), yticks=())
            plt.subplots_adjust(wspace=0.6)
     
            ## The figure save after the "Star again" will have a flag "sa" after the parameter (para, pmra, pmdec)
           
            plt.savefig(results_dir + '/PARALLAX_' + str(star_name) + '[' +str(para_lower)+ '_' + str(para_upper)+ ';' +str(len(para3_l))+ ']' + '_sa' + '.pdf')
            plt.tight_layout(h_pad=1.2, w_pad=1.2)
            plt.show()
            plt.close()

            
            again=str(input('Choose a new Parallax range? (y/n):')) 
           

        again=str(input('Introduce a new Proper Motion in RA range (y): '))
        while again=='y':
            
            IDGaia4_l,ID4_l,ra4_l,rae4_l,dec4_l,dece4_l,para4_l,parae4_l,pmra4_l,pmrae4_l,pmdec4_l,pmdece4_l,gmag4_l,BpRp4_l,BpG4_l ,GRp4_l =[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]
 
            print '\033[1;35m'
            print ("Make your selection!\n")
            print '\033[1;30m'

            # Define the range to select the data
            pmra_lower = input("PMRA RANGE >=")
            pmra_upper = input("PMRA RANGE<=")

            # Beginig of the loop
            for i in range(len(ID3_l)):
                if pmra3_l[i]>=pmra_lower and pmra3_l[i]<=pmra_upper: 
                    IDGaia4_l.append(IDGaia3_l[i])
                    ID4_l.append(ID3_l[i])
                    ra4_l.append(ra3_l[i])
                    rae4_l.append(rae3_l[i])
                    dec4_l.append(dec3_l[i])
                    dece4_l.append(dece3_l[i])
                    para4_l.append(para3_l[i])
                    parae4_l.append(parae3_l[i])
                    pmra4_l.append(pmra3_l[i])            
                    pmrae4_l.append(pmrae3_l[i])
                    pmdec4_l.append(pmdec3_l[i]) 
                    pmdece4_l.append(pmdece3_l[i])
                    gmag4_l.append(gmag3_l[i])
                    BpRp4_l.append(BpRp3_l[i])
                    BpG4_l.append(BpG3_l[i])
                    GRp4_l.append(GRp3_l[i])
            
            
            Mg4_l = gmag4_l + 5*np.log10(para4_l) - 10

            
            mu_para4,mu_pmra4,mu_pmdec4 = np.mean(para4_l),np.mean(pmra4_l),np.mean(pmdec4_l)

            print '\033[1;34m'
            print("Number of objects in this selection = {0}".format(len(para4_l)))
            print '\033[1;30m'
            
            print '\033[1;34m'
            print("Known value of PARALLAX = {0:.3f}[mas]".format(para_know))
            print("Known value of PMRA = {0:.3f}[mas/yr]".format(pmra_know))
            print("Known value of PMDEC ={0:.3f}[mas/yr]".format(pmdec_know))
            print '\033[1;30m'
            
            print '\033[1;34m'
            print("Mean value of PARALLAX = {0:.3f} [mas]".format(mu_para4))
            print("Mean value of PMRA = {0:.3f} [mas/yr]".format(mu_pmra4))
            print("Mean value of PMDEC = {0:.3f} [mas/yr]".format(mu_pmdec4))
            print '\033[1;30m'


            starRA_list4   = np.linspace(starRA, starRA, len(ra4_l)) 
            starDEC_list4  = np.linspace(starDEC,starDEC, len(dec4_l))  
            DECresult4     = np.asarray(dec4_l)
            DeltaRA4       = np.cos(DECresult4/(180/np.pi)) * (ra4_l - (starRA_list4))
            DeltaDEC4      = dec4_l - (starDEC_list4)
            Distance4      = np.sqrt(np.power(DeltaRA4,2) + np.power(DeltaDEC4,2))

            count4, bins4, bars4 = plt.hist(Distance4,30,color='white',histtype=u'barstacked', rwidth=0.8, lw=2, alpha=0.8) 
            plt.close()
            
            newbin4          = bins4 - np.min(Distance4) 
            dist_lower4      = np.linspace(np.min(newbin4), np.max(newbin4), len(count4))
            diff4            = np.diff(dist_lower4)
            high4            = np.linspace(np.min(diff4), np.max(diff4), len(count4))
            dist_high4       = dist_lower4 + high4
            ringarea4        = np.pi * np.power(dist_high4,2) - np.pi * np.power(dist_lower4,2)
            density4         = count4 / ringarea4
            count_err4       = np.sqrt(count4)
            density_err4     = count_err4 / ringarea4


            plt.figure(figsize=(12.5,6.5))

            ## Parallax 
            plt.subplot(2,4,1)
            count_para4, bins_para4, bars_para4 = plt.hist(para4_l, 35,  color='greenyellow',histtype=u'barstacked', rwidth=0.8, lw=2, alpha=0.8)
            plt.axvline(x=para_know, color='k', linestyle='--', lw=1)
            plt.xlabel(r'$\varpi $',fontsize=14)
            plt.ylabel('Number objects',fontsize=11)
            plt.tick_params(axis='both', which='major', labelsize=11)
            plt.subplots_adjust(hspace=0.5,wspace=0.6)

            ## Proper motion RA
            plt.subplot(2,4,2)
            count_pmra4, bins_pmra4, bars_pmra4 = plt.hist(pmra4_l, 35, color='aqua',histtype=u'barstacked', rwidth=0.8, lw=2, alpha=0.8)
            plt.axvline(x=pmra_know, color='k', linestyle='--', lw=1)
            plt.xlabel(r'$\mu_{\alpha*} $',fontsize=14)
            plt.tick_params(axis='both', which='major', labelsize=11)
            plt.subplots_adjust(hspace=0.5,wspace=0.6)

            ## Proper motion DEC
            plt.subplot(2,4,3)
            count_pmdec4, bins_pmdec4, bars_pmdec4 = plt.hist(pmdec4_l, 35,  color='royalblue',histtype=u'barstacked', rwidth=0.8, lw=2, alpha=0.8) 
            plt.axvline(x=pmdec_know, color='k', linestyle='--', lw=1)
            plt.xlabel(r'$\mu_{\delta} $',fontsize=14)
            plt.tick_params(axis='both', which='major', labelsize=11)
            plt.subplots_adjust(hspace=0.5,wspace=0.6)

            #Density profile
            plt.subplot(2,4,4)
            plt.plot(dist_high4,density4, '--.',  color='midnightblue', alpha=0.8)
            plt.errorbar(dist_high4,density4, yerr=density_err4, fmt='.', color='midnightblue', alpha=0.8, lw=.8)
            plt.xlabel('Distance (deg) ',fontsize=11) 
            plt.ylabel('Density (deg$^{-2}$)',fontsize=11) 
            plt.tick_params(axis='both', which='major', labelsize=11)
            plt.subplots_adjust(hspace=0.5,wspace=0.6)

            ## Spatial distribution
            plt.subplot(2,4,5)
            sc=plt.scatter(ra4_l,dec4_l, c=para4_l, marker='o',s=(np.array(para4_l)*4.5), cmap= 'YlGnBu', norm=None, vmin=np.min(para4_l), vmax=np.max(para4_l), alpha=None, linewidths=None, verts=None, edgecolors=None, hold=None, data=None)
            plt.plot(starRA,starDEC,'*',ms=6,color='darkviolet',alpha=0.65)
            plt.xlabel('RA',fontsize=11)
            plt.ylabel('DEC',fontsize=11)
            plt.gca().invert_xaxis()
            plt.colorbar(sc).remove()
            plt.tick_params(axis='both', which='major', labelsize=11)
            plt.subplots_adjust(wspace=0.5)

            ## Proper motion  distribution
            plt.subplot(2,4,6)
            sc=plt.scatter(pmra4_l,pmdec4_l, c=para4_l, marker='o',s=(np.array(para4_l)*4.5),cmap= 'YlGnBu', zorder=  35, norm=None, vmin=np.min(para4_l), vmax=np.max(para4_l), alpha=None, linewidths=None, verts=None, edgecolors=None, hold=None, data=None)
            plt.plot(pmra_know,pmdec_know,'*',ms=6, color='darkviolet',alpha=0.65)
            plt.colorbar(sc).remove()
            plt.xlabel(r'$\mu_{\alpha*} $',fontsize=14)
            plt.ylabel(r'$\mu_{\delta} $',fontsize=14)
            plt.tick_params(axis='both', which='major', labelsize=10)
            plt.subplots_adjust(wspace=0.5)

            ## CMD
            plt.subplot(2,4,7)
            sc=plt.scatter(BpRp4_l,Mg4_l, c=para4_l, marker='o',s=(np.array(para4_l)*4.5),cmap= 'YlGnBu', zorder=  35,norm=None, vmin=np.min(para4_l), vmax=np.max(para4_l), alpha=None, linewidths=None, verts=None, edgecolors=None, hold=None, data=None)
            plt.plot(bp_rp_mag_know,Mgstars,'*',ms=6,color='darkviolet',alpha=0.65)
            plt.colorbar(sc, aspect=25).set_label(label=r'$\varpi $', size=10)
            plt.xlabel('G$_{Bp}$ - G$_{Rp}$',fontsize=11)
            plt.ylabel('M$_{G}$',fontsize=11)
            plt.gca().invert_yaxis()
            plt.tick_params(axis='both', which='major', labelsize=10)
            plt.subplots_adjust(wspace=0.5)

            # Information of the target at this stage
            plt.subplot(2,4,8)
            plt.axis([0, 10, 0, 10])
            plt.text(0.5, 9, "Object= {0}".format(star_name), wrap=True, fontsize=10)
            plt.text(0.5, 8, "Number of objects= {0}".format(len(para4_l)), wrap=True, fontsize=10)
            plt.text(0.5, 6, "Known PARALLAX= {:.2f} +/- {:.2f} [mas]".format(para_know,parae_know), wrap=True, fontsize=10)
            plt.text(0.5, 5, "Known PMRA= {:.2f} +/- {:.2f} [mas/yr]".format(pmra_know,pmrae_know), wrap=True, fontsize=10)
            plt.text(0.5, 4, "Known PMDEC= {:.2f} +/- {:.2f} [mas/yr]".format(pmdec_know,pmdece_know), wrap=True, fontsize=10)
            plt.text(0.5, 2, "Mean PARALLAX = {:.2f} [mas]".format(mu_para4) , wrap=True, fontsize=10)
            plt.text(0.5, 1, "Mean PMRA= {:.2f} [mas/yr]".format(mu_pmra4), wrap=True, fontsize=10)
            plt.text(0.5, 0, "Mean PMDEC= {:.2f} [mas/yr]".format(mu_pmdec4), wrap=True, fontsize=10)
            plt.setp(plt.gca(), frame_on=False, xticks=(), yticks=())
            plt.subplots_adjust(wspace=0.6)
     

            plt.savefig(results_dir + '/PMRA_' + str(star_name) + '[' +str(pmra_lower)+ '_' + str(pmra_upper)+ ';' +str(len(pmra4_l))+ ']' + '_sa' + '.pdf')
            plt.tight_layout(h_pad=1.2, w_pad=1.2)#pad=0.6
            plt.show()
            plt.close()

            
            again=str(input('Choose a new Proper Motion in RA range? (y/n):')) 
            

        again=str(input('Introduce a new Proper Motion in DEC range (y): '))
        while again=='y':
            
            IDGaia5_l,ID5_l,ra5_l,rae5_l,dec5_l,dece5_l,para5_l,parae5_l,pmra5_l,pmrae5_l,pmdec5_l,pmdece5_l,gmag5_l,BpRp5_l,BpG5_l ,GRp5_l =[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]

            print '\033[1;35m'
            print ("Make your selection!\n")
            print '\033[1;30m'
        
        
            pmdec_lower = input("PMDEC RANGE>=")
            pmdec_upper = input("PMDEC RANGE<=")

            
            for i in range(len(ID4_l)):
                if pmdec4_l[i]>=pmdec_lower and pmdec4_l[i]<=pmdec_upper: 
                    IDGaia5_l.append(IDGaia4_l[i])
                    ID5_l.append(ID4_l[i])
                    ra5_l.append(ra4_l[i])
                    rae5_l.append(rae4_l[i])
                    dec5_l.append(dec4_l[i])
                    dece5_l.append(dece4_l[i])
                    para5_l.append(para4_l[i])
                    parae5_l.append(parae4_l[i])
                    pmra5_l.append(pmra4_l[i])            
                    pmrae5_l.append(pmrae4_l[i])
                    pmdec5_l.append(pmdec4_l[i]) 
                    pmdece5_l.append(pmdece4_l[i])
                    gmag5_l.append(gmag4_l[i])
                    BpRp5_l.append(BpRp4_l[i])
                    BpG5_l.append(BpG4_l[i])
                    GRp5_l.append(GRp4_l[i])
            
            
            Mg5_l = gmag5_l + 5*np.log10(para5_l) - 10

            
            mu_para5,mu_pmra5,mu_pmdec5 = np.mean(para5_l),np.mean(pmra5_l),np.mean(pmdec5_l)

            print '\033[1;34m'
            print("Number of objects in this selection = {0}".format(len(para5_l)))
            print '\033[1;30m'
            
            print '\033[1;34m'
            print("Known value of PARALLAX = {0:.3f}[mas]".format(para_know))
            print("Known value of PMRA =  {0:.3f}[mas/yr]".format(pmra_know))
            print("Known value of PMDEC =  {0:.3f}[mas/yr]".format(pmdec_know))
            print '\033[1;30m'
            
            print '\033[1;34m'
            print("Mean value of PARALLAX = {0:.3f} [mas]".format(mu_para5))
            print("Mean value of PMRA = {0:.3f} [mas/yr]".format(mu_pmra5))
            print("Mean value of PMDEC = {0:.3f} [mas/yr]".format(mu_pmdec5))
            print '\033[1;30m'
    
            
            starRA_list5   = np.linspace(starRA, starRA, len(ra5_l)) 
            starDEC_list5  = np.linspace(starDEC,starDEC, len(dec5_l))  
            DECresult5     = np.asarray(dec5_l)
            DeltaRA5       = np.cos(DECresult5/(180/np.pi)) * (ra5_l - (starRA_list5))
            DeltaDEC5      = dec5_l - (starDEC_list5)
            Distance5      = np.sqrt(np.power(DeltaRA5,2) + np.power(DeltaDEC5,2))
   
            count5, bins5, bars5 = plt.hist(Distance5, 30,  color='white',histtype=u'barstacked', rwidth=0.8, lw=2, alpha=0.8) 
            plt.close()
            
            newbin5          = bins5 - np.min(Distance5) 
            dist_lower5      = np.linspace(np.min(newbin5), np.max(newbin5), len(count5))
            diff5            = np.diff(dist_lower5)
            high5            = np.linspace(np.min(diff5), np.max(diff5), len(count5))
            dist_high5       = dist_lower5 + high5
            ringarea5        = np.pi * np.power(dist_high5,2) - np.pi * np.power(dist_lower5,2)
            density5         = count5 / ringarea5
            count_err5       = np.sqrt(count5)
            density_err5     = count_err5 / ringarea5
            
            

            plt.figure(figsize=(12.5,6.5))

            ## Parallax 
            plt.subplot(2,4,1)
            count_para5, bins_para5, bars_para5 = plt.hist(para5_l, 35,  color='greenyellow',histtype=u'barstacked', rwidth=0.8, lw=2, alpha=0.8)
            plt.axvline(x=para_know, color='k', linestyle='--', lw=1)
            plt.xlabel(r'$\varpi $',fontsize=14)
            plt.ylabel('Number objects',fontsize=11)
            plt.tick_params(axis='both', which='major', labelsize=11)
            plt.subplots_adjust(hspace=0.5,wspace=0.6)

            ## Proper motion RA
            plt.subplot(2,4,2)
            count_pmra5, bins_pmra5, bars_pmra5 = plt.hist(pmra5_l,  35,  color='aqua',histtype=u'barstacked', rwidth=0.8, lw=2, alpha=0.8)
            plt.axvline(x=pmra_know, color='k', linestyle='--', lw=1)
            plt.xlabel(r'$\mu_{\alpha*} $',fontsize=14)
            plt.tick_params(axis='both', which='major', labelsize=11)
            plt.subplots_adjust(hspace=0.5,wspace=0.6)

            ## Proper motion DEC
            plt.subplot(2,4,3)
            count_pmdec5, bins_pmdec5, bars_pmdec5 = plt.hist(pmdec5_l, 35,  color='royalblue',histtype=u'barstacked', rwidth=0.8, lw=2, alpha=0.8)
            plt.axvline(x=pmdec_know, color='k', linestyle='--', lw=1)
            plt.xlabel(r'$\mu_{\delta} $',fontsize=14)
            plt.tick_params(axis='both', which='major', labelsize=11)
            plt.subplots_adjust(hspace=0.5, wspace=0.6)

            #Density profile
            plt.subplot(2,4,4)
            plt.plot(dist_high5,density5, '--.', color='midnightblue', alpha=0.8)
            plt.errorbar(dist_high5,density5, yerr=density_err5, fmt='.', color='midnightblue', alpha=0.8, lw=.8)
            plt.xlabel('Distance (deg) ',fontsize=11) 
            plt.ylabel('Density (deg$^{-2}$)',fontsize=11) 
            plt.tick_params(axis='both', which='major', labelsize=11)
            plt.subplots_adjust(hspace=0.5, wspace=0.6)


            ## Spatial distribution
            plt.subplot(2,4,5)
            sc=plt.scatter(ra5_l,dec5_l, c=para5_l, marker='o', s=(np.array(para5_l)*4.5),cmap= 'YlGnBu', norm=None, vmin=np.min(para5_l), vmax=np.max(para5_l), alpha=None, linewidths=None, verts=None, edgecolors=None, hold=None, data=None)
            plt.plot(starRA,starDEC,'*',ms=6,color='darkviolet',alpha=0.65)
            plt.xlabel('RA',fontsize=11)
            plt.ylabel('DEC',fontsize=11)
            plt.gca().invert_xaxis()
            plt.colorbar(sc).remove()
            plt.tick_params(axis='both', which='major', labelsize=11)
            plt.subplots_adjust(wspace=0.5)


            ## Proper motion  distribution
            plt.subplot(2,4,6)
            sc=plt.scatter(pmra5_l,pmdec5_l, c=para5_l, marker='o',s=(np.array(para5_l)*4.5),cmap= 'YlGnBu', zorder=  35, norm=None, vmin=np.min(para5_l), vmax=np.max(para5_l), alpha=None, linewidths=None, verts=None, edgecolors=None, hold=None, data=None)
            plt.plot(pmra_know,pmdec_know,'*',ms=6, color='darkviolet',alpha=0.65)
            plt.colorbar(sc).remove()
            plt.xlabel(r'$\mu_{\alpha*} $',fontsize=14)
            plt.ylabel(r'$\mu_{\delta} $',fontsize=14)
            plt.tick_params(axis='both', which='major', labelsize=11)
            plt.subplots_adjust(wspace=0.5)


            ## CMD
            plt.subplot(2,4,7)
            sc=plt.scatter(BpRp5_l,Mg5_l, c=para5_l, marker='o',s=(np.array(para5_l)*4.5),cmap= 'YlGnBu', zorder=  35,norm=None, vmin=np.min(para5_l), vmax=np.max(para5_l), alpha=None, linewidths=None, verts=None, edgecolors=None, hold=None, data=None)
            plt.plot(bp_rp_mag_know,Mgstars,'*',ms=6,color='darkviolet',alpha=0.65)
            plt.colorbar(sc).set_label(label=r'$\varpi $', size=10)
            plt.xlabel('G$_{Bp}$ - G$_{Rp}$',fontsize=11)
            plt.ylabel('M$_{G}$',fontsize=11)
            plt.gca().invert_yaxis()
            plt.tick_params(axis='both', which='major', labelsize=11)
            plt.subplots_adjust(wspace=0.5)

            # Information of the target at this stage
            plt.subplot(2,4,8)
            plt.axis([0, 10, 0, 10])
            plt.text(0.5, 9, "Object= {0}".format(star_name), wrap=True, fontsize=10)
            plt.text(0.5, 8, "Number of objects= {0}".format(len(para5_l)), wrap=True, fontsize=10)
            plt.text(0.5, 6, "Known PARALLAX= {:.2f} +/- {:.2f} [mas]".format(para_know,parae_know), wrap=True, fontsize=10)
            plt.text(0.5, 5, "Known PMRA= {:.2f} +/- {:.2f} [mas/yr]".format(pmra_know,pmrae_know), wrap=True, fontsize=10)
            plt.text(0.5, 4, "Known PMDEC= {:.2f} +/- {:.2f} [mas/yr]".format(pmdec_know,pmdece_know), wrap=True, fontsize=10)
            plt.text(0.5, 2, "Mean PARALLAX = {:.2f} [mas]".format(mu_para5) , wrap=True, fontsize=10)
            plt.text(0.5, 1, "Mean PMRA= {:.2f} [mas/yr]".format(mu_pmra5), wrap=True, fontsize=10)
            plt.text(0.5, 0, "Mean PMDEC= {:.2f} [mas/yr]".format(mu_pmdec5), wrap=True, fontsize=10)
            plt.setp(plt.gca(), frame_on=False, xticks=(), yticks=())
            plt.subplots_adjust(wspace=0.6)

            plt.savefig(results_dir + '/PMDEC_' + str(star_name) + '[' +str(pmdec_lower)+ '_' + str(pmdec_upper)+ ';' + str(len(pmdec5_l))+ ']' + '_sa' + '.pdf')
            plt.tight_layout(h_pad=1.2, w_pad=1.2)#pad=0.6
            plt.show()
            plt.close()

            
            again=str(input('Choose a new Proper Motion in DEC range? (y/n):')) 
            
        ## Started the selection again, from the last loop. 
        # Selection the parameters in order: Parallax, PMRA and PMDEC.
        loop=str(input('Would you like to start again the selection? (y/n): '))

        while loop=='n':
            print('Well done! This are your results =D \n')

            print '\033[1;34m'
            print("Number of objects for this selection = {0}".format(len(pmra5_l)))
            print '\033[1;30m'
            
            print '\033[1;34m'
            print("Known value of PARALLAX = {0:.3f}[mas]".format(para_know))
            print("Known value of PMRA =  {0:.3f}[mas/yr]".format(pmra_know))
            print("Known value of PMDEC =  {0:.3f}[mas/yr]".format(pmdec_know))
            print '\033[1;30m'
            
            print '\033[1;34m'
            print("Mean value of PARALLAX = {0:.3f} [mas]".format(mu_para5))
            print("Mean value of PMRA = {0:.3f} [mas/yr]".format(mu_pmra5))
            print("Mean value of PMDEC = {0:.3f} [mas/yr]".format(mu_pmdec5))
            print '\033[1;30m'

            file= open(pathdata + '/OutputData_' + str(star_name) + '_sa' +  '.txt','w')
            file.write('#ID\t' 'RA\t' 'RA_Error\t' 'DEC\t' 'DEC_Error\t' 'PARALLAX\t' 'PARALLAX_Error\t' 'PMRA\t' 'PMRA_Error\t' ' PMDEC\t' ' PMDEC_Error\t' ' G_mag\t'  'Bp-Rp_mag\t' 'Bp-G_mag\t' 'G-Rp_mag\n' )
            for i in range(len(ID5_l)):
                file.write("%f %f %f %f %f %f %f %f %f %f %f %f  %f %f %f\n" % (ID5_l[i], ra5_l[i], rae5_l[i], dec5_l[i], dece5_l[i], para5_l[i], parae5_l[i], pmra5_l[i], pmrae5_l[i], pmdec5_l[i], pmdece5_l[i], gmag5_l[i], BpRp5_l[i], BpG5_l[i], GRp5_l[i]))
            file.close()           

            ## FINAL PLOT WITH THE RESULTS OF THE LAST SELECTION
            ## Estimate de error trought the weighted mean ecuation
            
            #Parallax#
            wpara = 1/np.power(parae5_l,2)
            Sigma_wparasqueare = np.sqrt(1/wpara)
            Xpara = (wpara*para5_l)/(wpara)
            Sigma_wparaT = np.sqrt(1/np.sum(wpara))
            XparaT = np.sum(wpara*para5_l)/np.sum(wpara)

            #PmRA
            wpmra = 1/np.power(pmrae5_l,2)
            Xpmra = (wpmra*pmra5_l)/(wpmra)
            Sigma_wpmraT = np.sqrt(1/np.sum(wpmra))
            XpmraT = np.sum(wpmra*pmra5_l)/np.sum(wpmra)

            #PmDEC
            wpmdec = 1/np.power(pmdece5_l,2)
            Sigma_wpmdecsqueare = np.sqrt(1/wpmdec)
            Xpmdec = (wpmdec*pmdec5_l)/(wpmdec)
            Sigma_wpmdecT = np.sqrt(1/np.sum(wpmdec))
            XpmdecT = np.sum(wpmdec*pmdec5_l)/np.sum(wpmdec)

 
            plt.figure(figsize=(12.5,6.5))
            
            ## Parallax 
            plt.subplot(2,4,1)
            count_para5, bins_para5, bars_para5 = plt.hist(para5_l,  35,  color='greenyellow',histtype=u'barstacked', rwidth=0.8, lw=2, alpha=0.8)
            plt.axvline(x=para_know, color='k', linestyle='--', lw=1)
            plt.xlabel(r'$\varpi $',fontsize=14)
            plt.ylabel('Number objects',fontsize=11)
            plt.tick_params(axis='both', which='major', labelsize=11)
            plt.subplots_adjust(hspace=0.5,wspace=0.6)

            ## Proper motion RA
            plt.subplot(2,4,2)
            count_pmra5, bins_pmra5, bars_pmra5 = plt.hist(pmra5_l, 35,  color='aqua',histtype=u'barstacked', rwidth=0.8, lw=2, alpha=0.8)
            plt.axvline(x=pmra_know, color='k', linestyle='--', lw=1)
            plt.xlabel(r'$\mu_{\alpha*} $',fontsize=14)
            plt.tick_params(axis='both', which='major', labelsize=11)
            plt.subplots_adjust(hspace=0.5,wspace=0.6)

            ## Proper motion DEC
            plt.subplot(2,4,3)
            count_pmdec5, bins_pmdec5, bars_pmdec5 = plt.hist(pmdec5_l,35,  color='royalblue',histtype=u'barstacked', rwidth=0.8, lw=2, alpha=0.8)
            plt.axvline(x=pmdec_know, color='k', linestyle='--', lw=1)
            plt.xlabel(r'$\mu_{\delta} $',fontsize=14)
            plt.tick_params(axis='both', which='major', labelsize=11)
            plt.subplots_adjust(hspace=0.5,wspace=0.6)
            
            #Density profile
            plt.subplot(2,4,4)
            plt.plot(dist_high5,density5, '--.',ms=4, color='midnightblue', alpha=0.8)
            plt.errorbar(dist_high5,density5, yerr=density_err5, fmt='.', color='midnightblue', alpha=0.8, lw=.8)
            plt.xlabel('Distance (deg) ',fontsize=11) 
            plt.ylabel('Density (deg$^{-2}$)',fontsize=11) 
            plt.tick_params(axis='both', which='major', labelsize=11)
            plt.subplots_adjust(wspace=0.8)
            
            
            ## Spatial distribution
            plt.subplot(2,4,5)
            sc=plt.scatter(ra5_l,dec5_l, c=para5_l, marker='o', s=(np.array(para5_l)*4.5),cmap= 'YlGnBu', norm=None, vmin=np.min(para5_l), vmax=np.max(para5_l), alpha=None, linewidths=None, verts=None, edgecolors=None, hold=None, data=None)
            plt.plot(starRA,starDEC,'*',ms=6,color='darkviolet',alpha=0.65)
            plt.xlabel('RA',fontsize=11)
            plt.ylabel('DEC',fontsize=11)
            plt.gca().invert_xaxis()
            plt.colorbar(sc).remove()
            plt.tick_params(axis='both', which='major', labelsize=11)
            plt.subplots_adjust(hspace=0.5, wspace=0.6)

            ## Proper motion  distribution
            plt.subplot(2,4,6)
            sc=plt.scatter(pmra5_l,pmdec5_l, c=para5_l, marker='o', s=(np.array(para5_l)*4.5),cmap= 'YlGnBu', zorder=  35, norm=None, vmin=np.min(para5_l), vmax=np.max(para5_l), linewidths=None, verts=None, edgecolors=None, hold=None, data=None)#, **kwargs cmap =cool, RdPu, magma,inferno, viridis, YlGnB
            plt.plot(pmra_know,pmdec_know,'s',ms=5, color='darkviolet',label='Known Value')
            plt.plot(XpmraT,XpmdecT,'o',ms=5,color='darkorange',label='This programme')
            plt.xlabel(r'$\mu_{\alpha*} $',fontsize=14)
            plt.ylabel(r'$\mu_{\delta} $',fontsize=14)
            plt.colorbar(sc).remove()
            plt.tick_params(axis='both', which='major', labelsize=11)
            plt.subplots_adjust(hspace=0.5, wspace=0.6)


            ## CMD
            plt.subplot(2,4,7)
            sc=plt.scatter(BpRp5_l,Mg5_l, c=para5_l, marker='o',s=(np.array(para5_l)*4.5),cmap= 'YlGnBu', zorder=  35,norm=None, vmin=np.min(para5_l), vmax=np.max(para5_l), alpha=None, linewidths=None, verts=None, edgecolors=None, hold=None, data=None)#, **kwargs cmap =cool, RdPu, magma,inferno, viridis, YlGnBu
            plt.plot(bp_rp_mag_know,Mgstars,'*',ms=6,color='darkviolet',alpha=0.65)
            plt.colorbar(sc, aspect=25).set_label(label=r'$\varpi $', size=10)#3,labelpad=10) 
            plt.xlabel('G$_{Bp}$ - G$_{Rp}$',fontsize=11)
            plt.ylabel('M$_{G}$',fontsize=11)
            plt.gca().invert_yaxis()
            plt.tick_params(axis='both', which='major', labelsize=11)
            plt.subplots_adjust(hspace=0.5, wspace=0.6)


            ## Text of the last results
            plt.subplot(2,4,8)
            plt.axis([0, 10, 0, 10])
            plt.text(0.5, 9, "Object= {0}".format(star_name), wrap=True, fontsize=10)
            plt.text(0.5, 8, "Number of objects= {0}".format(len(para5_l)), wrap=True, fontsize=10)
            plt.text(0.5, 7, "Known PARALLAX= {:.2f} +/- {:.2f} [mas]".format(para_know,parae_know), wrap=True, fontsize=10)
            plt.text(0.5, 6, "Known PMRA= {:.2f} +/- {:.2f} [mas/yr]".format(pmra_know,pmrae_know), wrap=True, fontsize=10)
            plt.text(0.5, 5, "Known PMDEC= {:.2f} +/- {:.2f} [mas/yr]".format(pmdec_know,pmdece_know), wrap=True, fontsize=10)
            plt.text(0.5, 4, "Mean PARALLAX = {:.2f} +/- {:.2f}  [mas]".format(XparaT,Sigma_wparaT) , wrap=True, fontsize=10)
            plt.text(0.5, 3, "Mean PMRA= {:.2f} +/- {:.2f}  [mas/yr]".format(XpmraT,Sigma_wpmraT), wrap=True, fontsize=10)
            plt.text(0.5, 2, "Mean PMDEC= {:.2f} +/- {:.2f}  [mas/yr]".format(XpmdecT,Sigma_wpmdecT), wrap=True, fontsize=10)
            plt.text(0.5, 0, "Square-> Known Value; Circle-> CEREAL", wrap=True, fontsize=8)
            plt.setp(plt.gca(), frame_on=False, xticks=(), yticks=())
            plt.subplots_adjust(wspace=0.6)

            plt.savefig(results_dir + '/OutputData_' + str(star_name) + '_sa' +'.pdf')
            plt.tight_layout(h_pad=1.2, w_pad=1.2)
            plt.show()
            plt.close()


            break

            ## Started the selection again, from the last loop. 
            # Selection the parameters in order: Parallax, PMRA and PMDEC.
            loop=str(input('Would you like to start again the selection? (y/n): '))
            
            while loop=='n':
                break

    
    ### Cluster clasification file ###
    # Here would need to insert a number to classified the presence of a cluster
    # around the stars that have been studied: Yes=1; No=0; Maybe=2. This would be
    # save in a file with the rest of the known information of each star. 
    
    print '\033[1;31m'
    inaclusternumb = input("This star is in a cluster? (Yes=1; No=0; Maybe=2) = " ) 
    print '\033[1;30m'

    outputnum = len(para5_l) # Final number of object after the selection process

    #Saving the date of each stars from the input file
    Star.append(star_name)
    A1.append(starRA)
    A2.append(starDEC)
    A3.append(para_know)
    A4.append(parae_know)
    A5.append(pmra_know)
    A6.append(pmrae_know)
    A7.append(pmdec_know)
    A8.append(pmdece_know)
    ClusterNum.append(inaclusternumb)
    B1.append(outputnum)
    
    file = open(pathdata + '/ClusterClasification' +  '.txt','w')
    file.write('#Name\t'  'RA\t'  'DEC\t'  'PARALLAX\t' 'PARALLAX_Error\t' 'PMRA\t' 'PMRA_Error\t' ' PMDEC\t' ' PMDEC_Error\t' 'ClusterClassNum\t' 'NumObjSelect\n' )
    for i in range(len(ClusterNum)):
        file.write("%s %f %f %0.4f %0.4f %0.4f %0.4f %0.4f %0.4f %0.1f %0.1f\n" % (Star[i],  A1[i], A2[i],  A3[i], A4[i], A5[i],  A6[i], A7[i], A8[i], ClusterNum[i],  B1[i]))
    file.close()


    NumStar += 1

now = time.time() #Time after it finished
print('It took: ', now-then, ' seconds')

