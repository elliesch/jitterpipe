def jitterpipe(dirpath, psrname, NANOdir, MJDint, clearoutput=True, mkfiles=True, 
               cal=True, zap=True, scrunch=True, timing=True, resids=True):
    
    '''
    Jitterpipe
    6/27/2016

    Modified version of Michael Lam's pipeline by Ellianna Schwab with help from Michael Lam and Scott Ransom.
    Michael Lam's original pipeline is here: http://astro.cornell.edu/~mlam/files/pipeline.py

    Reduction pipeline for any NANOGrav input data. Time scrunches to two intervals for jitter analysis.
    This assumes all the incoming parfiles are .gls.par!

    Takes arguments as 
    jitterpipe(dirpath, psrname, NANOdir, MJDint, clearoutput=True, mkfiles=True, 
               cal=True, zap=True, scrunch=True, timing=True, resids=True)

        dirpath is the directory that contains the pulsar's files. to run, store quasar calib files in a folder /folded.

        psrname is the name of the psr as displayed in par and sum.sm files, such as J1713+0747
        
        NANOdir is the location of the original fits files on bowser. Leave out the final slash

        MJDint takes in the brightest MJD date and creates a plot interval on that date

        clearoutput takes True or False, clears all prior reduction folders and files
        
        mkfiles takes True or False, creates cf and rf files from fits files

        cal takes True or False, makes calib files and stores them in calib/

        zap takes True or False, removes the RFI and makes zap files from the calibrated files

        scrunch takes True or False, scrunches the files to 10s and 80s subints, with 8 subchannels

        timing takes True or False, creates TOAs and tim files for  both sets of scrunched files

        resids takes True or False, creates plots to show the jitter on the MJD day with brightest flux for that object

    '''
    
    ## ==============================
    ## Imports and Definitions
    ## ==============================

    import os #provides uniform interface to a number of OS functions
    import sys #contains useful functions and variables
    import subprocess #processes that run as independent entities
    import glob #finds all pathnames matching a specified pattern

    import residuals as r #brings in Scott's residuals
    import numpy as np
    import matplotlib
    #matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    #%matplotlib inline
    import math

    OUTPUT_FRONT = "jitterpipe: "
    DIR = dirpath #ex. /nimrod1/eschwab/B1937_data/
    PARFILE = DIR + '%s_NANOGrav_11yv0.par' %psrname 
    TEMPLATEFILE = DIR + '%s.L-wide.PUPPI.11y.x.sum.sm' %psrname 

    def call(x): 
        subprocess.call(x,shell=True)
    def printer(x):
        print(OUTPUT_FRONT+x)
    def add_header (filename, line):
        with open(filename, 'r+') as f:
            content=f.read()
            f.seek(0,0)
            f.write(line.rstrip('\r\n') + '\n' + content)
        

    ## ==============================
    ## Choose Sections to Run
    ## ==============================

    printer("Running %s jitterpipe reduction pipeline" %psrname)

    CLEARFLAG = clearoutput
    FILEFLAG = mkfiles
    CALFLAG = cal
    ZAPFLAG = zap
    SPLITFLAG = scrunch
    TIMFLAG = timing
    RESIDFLAG = resids

    

    ## ==============================
    ## Clear old version
    ## ==============================

    ## Runs a bash script to rm old files and prints out deletion complete

    if CLEARFLAG:
        printer("Deleting old version")
        call("rm -rf %scal" %DIR)
        call("rm -rf %scalib" %DIR)
        call("rm -rf %szap" %DIR)
        call("rm -rf %stiming" %DIR)
        call("rm -rf %sproducts" %DIR)
        printer("Deletion complete")

        
    ## ==============================
    ## File Creation
    ## ==============================
    
    ##Generates rf and cf files from the raw fits files 
    ##First time scrunches to 10s subintervals, and then psradds combines all instances of each scan 
    ##and renames to .cf and .rf
        
    if FILEFLAG:
        printer("Tscrunching fits files")
        call("mkdir %sfits" %DIR)
        MJDarg='*' + str(MJDint) + '*'
        call("pam -u %sfits -e tfits --settsub 10 %s/%s.fits" %(DIR, NANOdir, MJDarg))
	letter=psrname[0]
        coords=psrname[1:]
        call("rename _%s _%s%s %sfits/*_????+*.tfits" %(coords, letter, coords, DIR))

        printer("Adding separated pulses into cf and rf files")
        
        #First separating out the MJDarray to loop over
        files = sorted(glob.glob("%sfits/*.tfits" %DIR))
        rfarray=[]
        cfarray=[]
        for f in files:
            col=str(f.split('_')[4])
            if col == 'cal':
                MJD = str(f.split('_')[1])
                Scannum = str(f.split('_')[3])
                cfarray.append([MJD, Scannum])
            else:
                MJD = str(f.split('_')[1])
                Scannum = str(f.split('_')[3])
		if [MJD, Scannum] not in rfarray:
                    rfarray.append([MJD, Scannum])

        #Now generating the rf files
        for scan in rfarray:
            call("psradd -o %sfolded/puppi_%s_%s_%s.11y.rf %sfits/puppi_%s_%s_%s_????.tfits" 
                 %(DIR, scan[0], psrname, scan[1], DIR, scan[0], psrname, scan[1]))
        
        #Now generating the cf files
        for scan in cfarray:
            call("psradd -o %sfolded/puppi_%s_%s_%s.11y.cf %sfits/puppi_%s_%s_%s_cal_????.tfits" 
                 %(DIR, scan[0], psrname, scan[1], DIR, scan[0], psrname, scan[1]))
                
        printer("rf and cf files generated")
        
    
    ## ==============================
    ## Calibration
    ## ==============================

    ## Generates cal files using pam on the cf files and pac on the rf files
    ## Setup so the cal and rf files should be in a folder called "folded"

    if CALFLAG:
#QUESTION FOR SCOTT
        #Should I be fully time scrunching the calibrators and do I need to save a copy of the database?
        #printer("Generating pulsar calibrators")
        call("mkdir %scal" %DIR)
        #call("cp %sfolded/*cf %scal/." %(DIR, DIR))
        #call("pam -m -T %scal/*cf" %DIR)

        #printer("Generating calibrator database")
        #call("pac -w -k %scal/caldatabase.txt -p %scal/" %(DIR, DIR, DIR))
        
        printer("Calibrating fits files")
        call("pac -Tx %sfolded/*rf -k %scal/caldatabase.txt" %(DIR, DIR)) #pac creates cal files along with database file
        printer("Calibration complete")

        printer("Moving calibrated files")
        call("mkdir %scalib" %DIR)
        call("mv %sfolded/*calib* %scalib/" %(DIR, DIR)) # in bash.cp *word* means anything before and anything after
        
        #printer("Fixing negative bandwidth")
        #call("pam -m --reverse_freqs calib/*calib")


    ## ==============================
    ## Zapping
    ## ==============================

    ## Generates cal and rf files using pam on the cf files and pac on the rf files
    
    if ZAPFLAG:
        printer("Removing RFI (this may take awhile)") #paz uses manual and automatic modes for interference excision
                                                       #for lots of files, this crashes due to memory problems
        files = sorted(glob.glob("%scalib/*.calib*" %DIR))
        for f in files: 
            call("paz -v -e zap -j 'zap median exp={$off:max-$off:min},zap median' %s" %f) 
            #this should generate zap files 
        call("mkdir %szap" %DIR)
        call("mv %scalib/*zap %szap/." %(DIR, DIR))
        printer("Zapping complete")


    ## ==============================
    ## Split into different frequency/time subintervals
    ## ==============================

    ## Given original files with 10 second intervals, gives files that are 80s subintervals + 8 subchannels and files
    ## that are 8 subchannels without being timescrunched.

    if SPLITFLAG:
        printer("Making 10s-->80s, 256-->8 channel files")

#         call("pam -e zap80F8 -t 8 -f 32 %szap/*zap" %DIR)
        call("pam -e zap80F8 --settsub 80 --setnchn 8 %szap/*zap" %DIR)

#         call("rename .zap80F8 _80F8.zap %szap/*zap80F8" %DIR)
        
        printer("Making 256-->8 channel, no tscrunch files")
        call("pam -e zapNTF8 --setnchn 8 %szap/*zap" %DIR)
#         call("rename .zapNTF8 _NTF8.zap %szap/*zapNTF8" %DIR)  
        
        printer("Splitting complete")


    ## ==============================
    ## Timing
    ## ==============================

    ## Calculates TOAs for both time subintervals and saves into tim files

    if TIMFLAG:
        printer("Producing timing solution")
        
        call("mkdir %stiming" %DIR)
        

        #Creating a for loop to pull out the unique MJD dates to an array
        files = sorted(glob.glob("%szap/*.zapNTF8" %DIR))
        MJDarray=[]
        for f in files:
            MJD = int(f.split('_')[1])
            if MJD not in MJDarray:
                MJDarray.append(MJD)
        
        #This will create tim files for the 80S SUBINT + 8 channel files
 
        #Creating forloop to loop over the MJD titles
        for date in MJDarray:
            call("pat -A FDM -e mcmc=0 -C chan -C subint -C snr -C wt -f 'tempo2 IPTA' -s %s %szap/puppi_%s_%s_????.11y.zap80F8 > %stiming/%s_%s_NANOGrav_11y_80F8.tim" %(TEMPLATEFILE, DIR, date, psrname, DIR, date, psrname))
            add_header('%stiming/%s_%s_NANOGrav_11y_80F8.tim' %(DIR, date, psrname), "MODE 1")
        
#         printer("Writing master tim file for 80s subint + 8 channel files")      
#         call("pat -A FDM -e mcmc=0 -C chan -C subint -C snr -C wt -f 'tempo2 IPTA' -s %s %szap/*.zap80F8 > %stiming/master_%s_NANOGrav_11y_80F8.tim" %(TEMPLATEFILE, DIR, DIR, psrname))
#         add_header('%smaster_%s_NANOGrav_11y_80F8.tim' %(DIR, psrname), "MODE 1")
        
        ## ==============================
        
        #This will create tim files for the 10S SUBINT (no tscrunch) + 8 channel files
        
        printer("Writing daily tim file for 10s subint + 8 channel files")  
                
        #Creating forloop to loop over the MJD titles
        for date in MJDarray:
            call("pat -A FDM -e mcmc=0 -C chan -C subint -C snr -C wt -f 'tempo2 IPTA' -s %s %szap/puppi_%s_%s_????.11y.zapNTF8 > %stiming/%s_%s_NANOGrav_11y_NTF8.tim" %(TEMPLATEFILE, DIR, date, psrname, DIR, date, psrname))
            add_header('%stiming/%s_%s_NANOGrav_11y_NTF8.tim' %(DIR, date, psrname), "MODE 1")
        
#         printer("Writing master tim file for 10s subint + 8 channel files")     
#         call("pat -A FDM -e mcmc=0 -C chan -C subint -C snr -C wt -f 'tempo2 IPTA' -s %s %szap/*.zapNTF8 > %stiming/master_%s_NANOGrav_11y_NTF8.tim" %(TEMPLATEFILE, DIR, DIR, psrname))
#         add_header('%smaster_%s_NANOGrav_11y_NTF8.tim' %(DIR, psrname), "MODE 1")  


        printer("Timing solution complete")


    ## ==============================
    ## Residuals
    ## ==============================


    if RESIDFLAG:
        printer("Running tempo")

        call("mkdir %sproducts" %DIR)
        
        #This turns off all the calibrators in the NANOGrav parfiles
        par=open(DIR + "%s_NANOGrav_11yv0.gls.par" %psrname, 'r')
        pardata=par.read()
        par.close()
        
        par=open(PARFILE, 'w')
        par.write(pardata.replace("  1  ", "  0  "))
        par.close
        

        #This runs tempo on 80s + 8 channel files, using the master tim files        
#         call("tempo -G -f %s %stiming/master_%s_NANOGrav_11y_80F8.tim " %(PARFILE, DIR, psrname) )
#         call("mv resid2.tmp %sproducts/resid_80F8.tmp" %DIR)
        
	date = str(MJDint)
#         print "tempo -G -f %s %stiming/%s_%s_NANOGrav_11y_80F8.tim" %(PARFILE, DIR, date, psrname) 
        call("tempo -G -f %s %stiming/%s_%s_NANOGrav_11y_80F8.tim" %(PARFILE, DIR, date, psrname) )
        call("mv resid2.tmp %sproducts/%s_resid_%s_80F8.tmp" %(DIR, psrname, date))

        #This runs tempo on 10s + 8 channel files, recursively using the daily tim files
#         print "tempo -G -f %s %stiming/%s_%s_NANOGrav_11y_NTF8.tim" %(PARFILE, DIR, date, psrname) 
        call("tempo -G -f %s %stiming/%s_%s_NANOGrav_11y_NTF8.tim" %(PARFILE, DIR, date, psrname) )
        call("mv resid2.tmp %sproducts/%s_resid_%s_NTF8.tmp" %(DIR, psrname, date))
        
        printer("Residuals generated")
    
        ## ==========================
    
#        printer("Generating Plots for specified MJD")
#        
#        #First Plotting 80s + 8 channel plot
#        
#        #Defining the residuals
#        x=r.read_residuals(filename= DIR + "products/%s_resid_%s_80F8.tmp" %(psrname, date))
#        
#        #Defining the colormap
#        from matplotlib.colors import LinearSegmentedColormap
#        interval=np.hstack([np.linspace(0.05, 0.25), np.linspace(0.35,0.9)])
#        colors=plt.cm.gist_rainbow(interval)
#        cmap=LinearSegmentedColormap.from_list('name', colors, 8)
#
#        #Defining Tick numbers for the colorbar
#        from matplotlib import ticker
#        tick_locator = ticker.MaxNLocator(nbins=9)
#        
#        #Creating the plot
#        fig,ax = plt.subplots(figsize=(17,8))
#        cax = ax.scatter(x.bary_TOA, x.prefit_sec, c=x.bary_freq, s=20, edgecolor='#262626', linewidth='0.35', cmap=cmap)
#        ax.set_title('MJD %s, All Frequency Bands, 80s subintervals, 8subchannels' %date, fontsize='16')
#        ax.set_xlim(MJDint, (MJDint + 1))
#        ax.set_ylim(-0.00001, 0.00001)
#        cb=fig.colorbar(cax)
#        cb.locator = tick_locator
#        cb.update_ticks()
#
#        plt.savefig(DIR + 'products/%s_%s_80s_8chan.png' %(psrname,date))
#        plt.show()
#        
#        #Next Plotting 10s + 8 channel plot
#        
#        #Defining the residuals
#        y=r.read_residuals(filename= DIR + "products/%s_resid_%s_NTF8.tmp" %(psrname, date))
#        
#        #Creating the plot
#        fig,ax = plt.subplots(figsize=(17,8))
#        cax = ax.scatter(y.bary_TOA, y.prefit_sec, c=y.bary_freq, s=20, edgecolor='#262626', linewidth='0.35', cmap=cmap)
#        ax.set_title('MJD %s, All Frequency Bands, 10s subintervals, 8subchannels' %date, fontsize='16')
#        ax.set_xlim(MJDint, (MJDint + 1))
#        ax.set_ylim(-0.00001, 0.00001)
#        cb=fig.colorbar(cax)
#        cb.locator = tick_locator
#        cb.update_ticks()
#
#        plt.savefig(DIR + 'products/%s_%s_10s_8chan.png' (%psrname, date))
#        plt.show()
#        
#        printer("Plots generated")


    printer("Jitterpipe complete")
