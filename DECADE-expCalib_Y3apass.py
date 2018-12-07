#/usr/bin/env python
#
"""
    expCalib.py THIS VERSION ONLY works for Carina2-3
    Express Calibration, 
    This code will estimate the zero-points   
    v.3 Apr20, 2016
    NOW using apass_2massInDES.sorted.csv via APASS/2MASS.
    v.2 Feb25, 2016:
    Now use APASS Dr7 and tested with ALex.. MagsLite
    Date Thu Sep 24, 2015
    NOTE that this catalog is only for the officical DES-foot print, some ccd will have no ZP
    Example:   
    expCalib_Y3apass.py --help

    GW-expCalib_Y3apass.py -s db-desoper --expnum 475956 --reqnum 04 --attnum 11 
    
    """
import os
import numpy as np
import argparse
import sys
import healpy as hp
import pandas as pd
##################################

def main():
    print " Start with DECADE-expCalib.py \n"
    """Create command line arguments"""
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--caldir', help='caldir is the calibration directory',
                        default='/des002/devel/emorgan2/APASS_TWOMASS/', type=str)
    parser.add_argument('--dir', help='dir is the production directory',
                        default='/archive_data/desarchive/DEC/finalcut/Y5A1/HITS/', type=str)
    parser.add_argument('--outdir', help='dir is the production directory', default='.', type=str)
    parser.add_argument('--expnum', help='expnum is queried', default=288940, type=int)
    parser.add_argument('--reqnum', help='reqnum is queried', default=3505, type=str)
    parser.add_argument('--attnum', help='attnum is queried', default=1, type=int)
    parser.add_argument('--magType', help='mag type to use (mag_psf, mag_auto, mag_aper_8, ...)', default='mag_psf')
    parser.add_argument('--sex_mag_zeropoint',
                        help='default sextractor zeropoint to use to convert fluxes to sextractor mags \
                             (mag_sex = -2.5log10(flux) + sex_mag_zeropoint)',
                        type=float, default=25.0)
    parser.add_argument('--verbose', help='verbosity level of output to screen (0,1,2,...)', default=0, type=int)
    parser.add_argument('--debug', help='debugging option', dest='debug', action='store_true', default=False)
                        
    args = parser.parse_args()
                    
    if args.verbose > 0: print args

    # Create all *std files
 #   getallccdfromAPASS92MASS(args)
    #exit()
    # -- ADDED NEW
  #  doset(args)
    #exit()
    # WHEN NEEDED
    # plot ra,dec of Sex-vs Y2Q1 for each CCD
    if args.verbose >0 :
        plotradec_sexvsY2Q1(args)

    # Estimate 3sigma Clipped Zeropoint for each CCD
#    sigmaClipZP(args)
#    exit()

#    sigmaClipZPallCCDs(args)
#    exit()

    ZP_OUTLIERS(args)
    exit()

    # --
    Onefile(args)

    # --
    # plot ra,dec of matched stars for ALL CCDs
    #" Comment this line for grid production "
#    plotradec_ZP(args)


def are_you_here(fname):
    # Check if file exists, if not - exit
    if not os.path.exists(fname):
        print '%s does not seem to exist... exiting now...' % fname
        sys.exit(1)
    else:
        pass

def read_catlist(args):
    # Open catalog created with make_red_catlist.py 
    catlistFile="""D%08d_r%sp%02d_red_catlist.csv""" % (args.expnum, str(args.reqnum), args.attnum)

    # Check if file exists, exit if it does not                                                                                                                                                                
    are_you_here(catlistFile)

    # Read and return the catalog  
    return catlistFile, pd.read_csv(catlistFile)

def get_corners(dat, off=0.):
    minra = dat[['RA_CENT', 'RAC1', 'RAC2', 'RAC3', 'RAC4']].min() - off
    maxra = dat[['RA_CENT', 'RAC1', 'RAC2', 'RAC3', 'RAC4']].max() + off
    mindec = dat[['DEC_CENT', 'DECC1', 'DECC2', 'DECC3', 'DECC4']].min() - off
    maxdec = dat[['DEC_CENT', 'DECC1', 'DECC2', 'DECC3', 'DECC4']].max() + off

    return minra, maxra, mindec, maxdec

##################################
# Get data from  PROD tables EXPOSURE IMAGE, WDF, and CATALOG,
# then Convert the Fits table to csv and save it
def doset(args):
    import csv
    
    if args.verbose >0 : print args

    catlistFile, data = read_catlist(args)

    for i, row in data.iterrows():
        if os.path.isfile(row['FILENAME']):
            Read_Sexcatalogfitstocsv(args,row['FILENAME'],row['BAND'])
            
            minra, maxra, mindec, maxdec = get_corners(row, off=0.)
            
            desipixc = getipix(128, row['RA_CENT'], row['DEC_CENT'])
            desipix1 = getipix(128, row['RAC1'], row['DECC1'])
            desipix2 = getipix(128, row['RAC2'], row['DECC2'])
            desipix3 = getipix(128, row['RAC3'], row['DECC3'])
            desipix4 = getipix(128, row['RAC4'], row['DECC4'])

            desipix12 = getipix(128, row['RAC1'], row['DEC_CENT'])
            desipix23 = getipix(128, row['RA_CENT'], row['DECC2'])
            desipix34 = getipix(128, row['RAC3'], row['DEC_CENT'])
            desipix14 = getipix(128, row['RA_CENT'], row['DECC4']) 

            desipixlist = desipixc,desipix1,desipix2,desipix3,desipix4,desipix12,desipix23,desipix34,desipix14
            desipixlist = uniqlist(desipixlist)
            
            matchlistout="""%s_match.csv""" % (row['FILENAME'])
            matchlistout = args.outdir+'/'+matchlistout.split('/')[-1]
            objlistFile ="""%s_Obj.csv"""   % (row['FILENAME'])
            objlistFile = args.outdir+'/'+objlistFile.split('/')[-1]
            stdlistFile ="""%s_std.csv"""   % (row['FILENAME'])        
            stdlistFile = args.outdir+'/'+stdlistFile.split('/')[-1]

            are_you_here(objlistFile)
            are_you_here(stdlistFile)

            matchSortedStdwithObsCats(stdlistFile, objlistFile, matchlistout,
                                      stdracol=1, stddeccol=2,
                                      obsracol=1, obsdeccol=2,
                                      matchTolArcsec=1.0, verbose=2)
 
      
##################################
#get_data_home for NOW it is for all CCDs:
#
def Wget_data_home(args):
    import csv
    import glob
    import sys

    catname="""D%08d_%s_%s_r%sp%02d_fullcat.fits""" % (args.expnum,"%","%",args.reqnum,args.attnum)
    myfile="""D%08d_*_r%sp%02d_fullcat.fits""" % (args.expnum,args.reqnum,args.attnum)

    #Check first if file exists...
    if  glob.glob(myfile):
        #Print '%s does seem to exist... exiting now...' % catname
        print "relevant cat files already exist in the current directory... no need to wget..."
        #sys.exit(1)
        return 1
    else:
        print "relevant cat files are not in directory... wgetting them from archive..."
        sys.exit(1)
    
##################################
# Quick Read SEX_table filemane_fullcat.fits then select subsame
# and write it as filemane_fullcat.fits_Obj.csv
def Read_Sexcatalogfitstocsv(args,fitsname,band):

    import fitsio
    import string
    import math
    import csv
 
    catFilename=fitsname
    outFile="""%s_Obj.csv""" % (catFilename)
    outFile = args.outdir+'/'+outFile.split('/')[-1]

    extension=2
    hdr = ["OBJECT_NUMBER","RA","DEC","MAG","MAGERR","ZEROPOINT","MAGTYPE","BAND"]

    magType = args.magType.upper()
    magType = magType.strip()
    fluxType = magType.replace('MAG','FLUX')
    fluxerrType = magType.replace('MAG','FLUXERR') 

    SEXdata=[]
    columns=['NUMBER','ALPHAWIN_J2000','DELTAWIN_J2000',fluxType,fluxerrType,'SPREAD_MODEL','SPREADERR_MODEL','FWHM_WORLD', 'CLASS_STAR', 'FLAGS']

    Fdata = fitsio.read(catFilename,  columns=columns, ext=extension)[:]
    #w0=( Fdata['FLUX_PSF'] > 2000) _OLD  
    w0=( Fdata['FLUX_PSF'] > 1000)  #NEW
    w1=( Fdata['FLAGS'] <= 3)
    #w2=( (Fdata['CLASS_STAR'] > 0.8 ) | (np.abs(Fdata['SPREAD_MODEL'] + 3.*Fdata['SPREADERR_MODEL'] <0.003 )))
    # NEW May16,16
    w2=( (Fdata['CLASS_STAR'] > 0.8 ) & (Fdata['SPREAD_MODEL']  <0.01 ) )

    SEXdata = Fdata[w0 & w1 & w2]
    SEXdata = SEXdata[np.argsort(SEXdata['ALPHAWIN_J2000'])]
    fwhm_arcsec=3600.*SEXdata['FWHM_WORLD']
    mag= -2.5*np.log10(SEXdata[fluxType]) + args.sex_mag_zeropoint
    magerr = (2.5/math.log(10.))*(SEXdata[fluxerrType]/SEXdata[fluxType])
    zeropoint=args.sex_mag_zeropoint*(SEXdata[fluxType]/SEXdata[fluxType])
  
    with open(outFile,'w') as csvFile:
            writer = csv.writer(csvFile,delimiter=',',  quotechar='|',
                                lineterminator='\n', quoting=csv.QUOTE_MINIMAL)

            writer.writerow(hdr)
            line=[]
            for i in range(SEXdata.size):
                line=SEXdata['NUMBER'][i],SEXdata['ALPHAWIN_J2000'][i], SEXdata['DELTAWIN_J2000'][i], mag[i], magerr[i], zeropoint[i], magType,band 
                writer.writerow(line)

    SEXdata=[]
    Fdata=[]

##################################
def radec2thetaphi(ra, dec):
    import numpy as np
    return (90-dec)*np.pi/180., ra*np.pi/180.

##################################
#for DES--nside=128
def getipix(nside,ra,dec):
    import healpy as hp
    #nside=128
    theta, phi = radec2thetaphi(ra, dec)
    ipix = hp.pixelfunc.ang2pix(nside, theta, phi, nest=True)
    return ipix

##################################
#Uniq List with order preserving
def uniqlist(seq): 
   noDupes = []
   [noDupes.append(i) for i in seq if not noDupes.count(i)]
   return noDupes



###################################
#NEW  July 14,2016 
#This is a FULL SKY 
#/data/des20.b/data/sallam/pyPSM_Year2/TWOMASS/ALL-2MASS/2arcsec 
#catalog now in stash cache in dcache
# /pnfs/des//persistent/stash/ALLSKY_STARCAT/....
#apass_TWO_MASS_*.csv  in totall 768 files  
#filters are u_des,g_des,r_des,i_des,z_des,Y_des
#################################
def getallccdfromAPASS92MASS(args):
    import pandas as pd
    import string, glob

    # Read the catalog  
    catlistFile, data = read_catlist(args)
    
    # Get band 
    BAND = data['BAND'][0]
    
    # Transform sky coordinates to pixel cooridnates 
    data['desipixc'] = getipix(8, data['RA_CENT'], data['DEC_CENT'])

    data['desipix1'] = getipix(8, data['RAC1'], data['DECC1']) 
    data['desipix2'] = getipix(8, data['RAC2'], data['DECC2']) 
    data['desipix3'] = getipix(8, data['RAC3'], data['DECC3']) 
    data['desipix4'] = getipix(8, data['RAC4'], data['DECC4']) 

    data['desipix12'] = getipix(8, data['RAC1'], data['DEC_CENT']) 
    data['desipix23'] = getipix(8, data['RA_CENT'], data['DECC2']) 
    data['desipix34'] = getipix(8, data['RAC3'], data['DEC_CENT']) 
    data['desipix14'] = getipix(8, data['RA_CENT'], data['DECC4']) 
 
    # Get unique values for desipix__
    desipixlist = pd.unique(data[['desipixc','desipix1','desipix2','desipix3','desipix4','desipix12','desipix23','desipix34','desipix14']].values.ravel())
    
    # RA for standards
    stdRA = np.std(data['RA_CENT'])
    
    # Round ra ??? why?
    if ( stdRA >20 ) :
        data['RA_CENT'] = [roundra(x) for x in data['RA_CENT']]
        data['RAC1']    = [roundra(x) for x in data['RAC1']]
        data['RAC2']    = [roundra(x) for x in data['RAC2']]
        data['RAC3']    = [roundra(x) for x in data['RAC3']]
        data['RAC4']    = [roundra(x) for x in data['RAC4']]

    # Define corners of frame??
    minra =  data[['RA_CENT', 'RAC1', 'RAC2', 'RAC3', 'RAC4']].min().min() - 0.1
    mindec = data[['DEC_CENT','DECC1', 'DECC2', 'DECC3', 'DECC4']].min().min() - 0.1
    maxra =  data[['RA_CENT', 'RAC1', 'RAC2', 'RAC3', 'RAC4']].max().max() + 0.1
    maxdec = data[['DEC_CENT','DECC1', 'DECC2', 'DECC3', 'DECC4']].max().max() + 0.1
   
    # Create string with output file name
    outfile = """STD%s""" % catlistFile
    
    # Create empty list/tables that will be used to create output file
    good_data  = []

    BANDname = BAND+"_des"

    for i in  desipixlist:
        myfile="""/des002/devel/emorgan2/APASS_TWOMASS/apass_TWO_MASS_%d.csv""" %i

        are_you_here(myfile)

        df= pd.read_csv(myfile)
        good_data.append(df)

    # Put some limits on data, 
    # cut the frame definded with min/max Ra/De from catalog of standards

    # reorganize - list of pd frames to one pd frame
    chunk = pd.concat(good_data, ignore_index=True).sort(['RAJ2000_APASS'], ascending=True)
    datastd = chunk.loc[(chunk['RAJ2000_2MASS'] > minra) & (chunk['RAJ2000_2MASS'] < maxra) &
                        (chunk['DEJ2000_2MASS'] > mindec) & (chunk['DEJ2000_2MASS'] < maxdec) &
                        (chunk[BANDname] > 0)]
    
    datastd1= pd.DataFrame({'MATCHID':datastd['MATCHID'],
                            'RA':datastd['RAJ2000_2MASS'],
                            'DEC':datastd['DEJ2000_2MASS'],
                            'WAVG_MAG_PSF':datastd[BANDname]})

    col = ["MATCHID", "RA","DEC", "WAVG_MAG_PSF"]  # could be replaced by datastd1 keywords
   
    # Save to csv file
    datastd1.to_csv(outfile, columns=col, sep=',', index=False)

    # Create std files - compare standard to observations (?)
    for i, row in data.iterrows():
        stdlistFile = """%s_std.csv"""   % row['FILENAME']
        stdlistFile = args.outdir+'/'+stdlistFile.split('/')[-1] 
        minra, maxra, mindec, maxdec = get_corners(row, off=0.1)
        df = datastd1.loc[(datastd1['RA'] > minra) & (datastd1['RA'] < maxra) &
                          (datastd1['DEC'] > mindec) & (datastd1['DEC'] < maxdec)].sort(['RA'], ascending=True)
        df.to_csv(stdlistFile, columns=col, sep=',', index=False)

##################################
#Matching FILES MAKE SURE the Cols are still in the same order
def matchSortedStdwithObsCats(f1, f2, outfile,
                              stdracol=1, stddeccol=2, obsracol=1, obsdeccol=2, matchTolArcsec=1,verbose=1):
    
    # Open the standard star CSV file and read the first line as list...
    fs = pd.read_csv(f1, delimiter=',')
    
    # Open CSV file of observed data and read it line by line...
    fd2 = open(f2)
    h2 = fd2.readline().strip().split(',')  # read header

    # Create an empty data frame that will be saved later as CSV file...
    #  Note that the column names from the standard star file
    #  now have a suffix of "_1", and that column names from
    #  the observed star file now have a suffix of "_2".
    outputHeader = ['MATCHID']
    for colhead in list(fs):
        outputHeader.append(colhead.upper() + '_1')
    for colhead in h2:
        outputHeader.append(colhead.upper() + '_2')
    out = pd.DataFrame(index=[], columns=outputHeader)
    
    # Initialize some variables
    done_obs = False
    tol = matchTolArcsec / 3600.0  # sky angular separation tolerance (in degrees)
    tol2 = tol*tol  # square for tol

    linecnt = 0  # line count for Obj file
    m_id = 0  # match id

    # Loop through file of observed data...
    while not done_obs:

        # Increment line count from observed data file...
        linecnt += 1
        if linecnt/1000.0 == int(linecnt/1000.0) and verbose > 1 :
            print '\r'+'Progress (lines read from observed catalog):  ',linecnt,
            sys.stdout.flush()

        # Read line from observed data file...
        obsline = fd2.readline().strip().split(',')

        if obsline == ['']:
            done_obs = True
        else:
            obsra = float(obsline[obsracol])
            obsdec = float(obsline[obsdeccol])
            cosd = np.cos(np.radians(obsdec))
            # Check if there is a match for object coordinates in standard stars file
            # calculate distance and check if it is within given limit
            delta2 = np.array((obsra - fs['RA'])*(obsra - fs['RA'])*cosd*cosd+(obsdec-fs['DEC'])*(obsdec-fs['DEC']))
            mtch = fs.iloc[np.where(delta2 < tol2)]
            
            if not mtch.empty:
                for m in mtch.iterrows():
                    m_id += 1
                    out_row = pd.DataFrame([[m_id] + mtch.values.tolist()[0] + obsline], columns=outputHeader)
                    out = out.append(out_row, ignore_index=True)
 
    # Change dtype so output file looks nice
    out.MATCHID = out.MATCHID.astype(int)
    out.MATCHID_1 = out.MATCHID_1.astype(int)
    out.OBJECT_NUMBER_2 = out.OBJECT_NUMBER_2.astype(int)
    
    # Drop matches to output file
    out.to_csv(outfile, index=False)
    
    # close the input object file
    fd2.close()
  

def get_new_cuts(in_file):
    from astropy.stats import sigma_clip

    # add new cuts for Apass9-2mass data set
    mdata = pd.read_csv(in_file, delimiter=',')

    # New CUTS
    data = mdata.loc[(mdata['MAG_2'] - mdata['WAVG_MAG_PSF_1'] - 25.0 < -10) &
                     (mdata['MAG_2'] - mdata['WAVG_MAG_PSF_1'] - 25.0 > -40)]

    if data.shape[0] != 0:
        delt_mag_data = data['MAG_2'] - data['WAVG_MAG_PSF_1'] - 25.0
        filtered_data = sigma_clip(delt_mag_data, sigma=3, iters=3, cenfunc=np.mean, copy=True)
        NumStarsClipped = filtered_data.count()
        NumStarsAll = len(delt_mag_data)

        if NumStarsClipped > 2:
            sigclipZP = np.mean(filtered_data)
            stdsigclipzp = np.std(filtered_data) / np.sqrt(NumStarsClipped)
        else:
            sigclipZP = -999
            stdsigclipzp = -999
    else:
        sigclipZP = -999
        stdsigclipzp = -999
        NumStarsClipped = 0
        NumStarsAll = 0

    return sigclipZP, stdsigclipzp, NumStarsClipped, NumStarsAll

##################################
# Get 3sigma clipped Zero point and iterater
##################################
def sigmaClipZP(args):
    if args.verbose > 0:
        print args
    
    catlistFile, data1 = read_catlist(args)
    
    ZeroListFile = """Zero_D%08d_r%sp%02d.csv""" % (args.expnum, args.reqnum, args.attnum)
    fout = open(ZeroListFile,'w')
    hdr = "FILENAME,Nall,Nclipped,ZP,ZPrms,magType\n"
    fout.write(hdr)

    for row in data1['FILENAME']:
        catFilename = os.path.basename(row)
        matchListFile = "%s_match.csv" % (catFilename)        

        are_you_here(matchListFile)

        sigclipZP, stdsigclipzp, NumStarsClipped, NumStarsAll = get_new_cuts(matchListFile)

        line = """%s,%d,%d,%f,%f,%s""" % (row, NumStarsAll, NumStarsClipped, sigclipZP, stdsigclipzp, args.magType)
        fout.write(line + '\n')

    fout.close()        

    ZeroListFile="""Zero_D%08d_r%sp%02d.csv""" % (args.expnum,args.reqnum,args.attnum)
    are_you_here(ZeroListFile)

    MergedFile="""Merged_D%08d_r%sp%02d.csv""" % (args.expnum,args.reqnum,args.attnum)
    jointwocsv(catlistFile,ZeroListFile,MergedFile) 
    
##################################
#Given two csv join both output to MergedFile
def jointwocsv(file1, file2, MergedFile):    
    import csv
    from collections import OrderedDict

    filenames = file1, file2
    data = OrderedDict()
    fieldnames = []
    for filename in filenames:
        with open(filename, "rb") as fp: 
            reader = csv.DictReader(fp)
            fieldnames.extend(reader.fieldnames)
            for row in reader:
                data.setdefault(row["FILENAME"], {}).update(row)

    fieldnames = list(OrderedDict.fromkeys(fieldnames))
    with open(MergedFile, "wb") as fp:
        writer = csv.writer(fp)
        writer.writerow(fieldnames)
        for row in data.itervalues():
            writer.writerow([row.get(field, '') for field in fieldnames])
#
##################################
# Modified at 10/09/2017 
# changed from 357 to 350 
#
# round ra 360.0--0.0
def roundra(ra):
#    if ra < 357:
    if ra < 356:
        return ra
    else:
        return ra-360.

##################################
# Get 3sigma clipped Zero point and iterater
##################################
def sigmaClipZPallCCDs(args):
    import os,glob


    if args.verbose >0 : print args

    allZPout = """allZP_D%08d_r%sp%02d.csv""" % (args.expnum, args.reqnum, args.attnum)
    stdfile = """STDD%08d_r%sp%02d_red_catlist.csv""" % (args.expnum, args.reqnum,args.attnum)
    objfile = """ObjD%08d_r%sp%02d_red_catlist.csv""" % (args.expnum, args.reqnum,args.attnum)
    outfile = """OUTD%08d_r%sp%02d_red_catlist.csv""" % (args.expnum, args.reqnum,args.attnum)
    
    # Read STD file, sort it and rewrite
    stddf = pd.read_csv(stdfile).sort(['RA'], ascending=True)
    stddf.to_csv(stdfile, sep=',', index=False)
    path = './'
    all_files = glob.glob(os.path.join(path, "*Obj.csv"))     
    df = pd.concat((pd.read_csv(f) for f in all_files)).sort(['RA'], ascending=True)

    #read all file and sort and save 
    df.to_csv(objfile, sep=',', index=False)

    matchSortedStdwithObsCats(stdfile, objfile, outfile, 
                              stdracol=1, stddeccol=2,
                              obsracol=1, obsdeccol=2,
                              matchTolArcsec=1, verbose=2)

    are_you_here(outfile)
    sigclipZP, stdsigclipzp, NumStarsClipped, NumStarsAll = get_new_cuts(outfile)

    hdr = "EXPNUM,REQNUM,ATTNUM,NumStarsAll,NumStarsClipped,sigclipZP,stdsigclipzp\n"    
    line = """%d,%s,%d,%d,%d,%f,%f""" % (args.expnum, args.reqnum, args.attnum, NumStarsAll, NumStarsClipped,
                                         sigclipZP, stdsigclipzp)
    print line
    fout = open(allZPout, 'w')
    fout.write(hdr)
    fout.write(line + '\n')

###################################
# Detect OUTLIERS- and update files
# NEEDS WORK--
# NEED to ADD LOF-see JoinZP3
###################################
def ZP_OUTLIERS(args):
    import os,sys,glob
    import math
    import sklearn
    from sklearn.neighbors import NearestNeighbors
    import matplotlib.pyplot as plt
    import matplotlib as mpl

    if args.verbose >0 : print args
    
    allZeroFile = """allZP_D%08d_r%sp%02d.csv""" % (args.expnum, args.reqnum, args.attnum)
    
    MergedFile="""Merged_D%08d_r%sp%02d.csv""" % (args.expnum, args.reqnum, args.attnum)
    are_you_here(MergedFile)

    fout = """Merg_allZP_D%08d_r%sp%02d.csv""" % (args.expnum, args.reqnum, args.attnum)

    df1 = pd.read_csv(MergedFile)
    df2 = np.genfromtxt(allZeroFile,dtype=None,delimiter=',',names=True)

    w0= ( df1['Nclipped'] < 4 ) | ( df1['ZP'] < -100 ) | ( df1['ZPrms'] > 0.3 )
    df1['NewZP']=np.where(w0 , df2['sigclipZP'],df1['ZP'])
    df1['NewZPrms']=np.where(w0, df2['stdsigclipzp'],df1['ZPrms'])
    df1['NewZPFlag1']=np.where(w0, np.int16(1),np.int16(0))
                    
    df1['DiffZP']=df1['NewZP']- df1.NewZP.median()

    #Currently Diff ZP=0.3mag That is TOO MUCH    
    w2= ( abs(df1['DiffZP']) < 0.15 )
    df1['NewZPFlag2']=np.where(w2 ,np.int16(0),np.int16(-1000))        
    df1['Percent1']=100.0*np.count_nonzero(df1['NewZPFlag2'])/len(df1['NewZP'])

    #if 20% CCD (i.e 12 CCDs out 60)
    w3= ( df1['Percent1'] >= 20 )  
    df1['NewZPFlag3']=np.where(w3 ,np.int16(-9999),np.int16(0))

    df1['NewZPFlag']=df1['NewZPFlag1']+df1['NewZPFlag2']+df1['NewZPFlag3']
    
    df1['DiffZP1']=1000.0*df1['DiffZP']                                            

    df1.to_csv(MergedFile,sep=',',index=False)

    # That will not work around 0-360
    # New
    w= ( df1['RA_CENT'] >150 ) 
    df1['roundRA_CENT']=np.where(w, 360.-df1['RA_CENT'], df1['RA_CENT'])
 
    #As still listed in *_immask.fits Fits heareder
    #PIXSCAL1=0.27 / [arcsec/pixel] Pixel scale, axis1
    #PIXSCAL1=0.27 / [arcsec/pixel] Pixel scale, axis2
    df1['x'] = (1./0.27)*3600*(df1['roundRA_CENT']-df1['roundRA_CENT'].median())*math.cos(math.radians(df1['DEC_CENT'].median()))                                                                                
    df1['y'] = (1./0.27)*3600*(df1['DEC_CENT'] - df1['DEC_CENT'].median())

    cols=['FILENAME','EXPNUM','CCDNUM','NewZP','NewZPrms','NewZPFlag']
    df1.to_csv(fout,columns=cols,sep=',',index=False)
#
##################################
# Read SEX_table filemane_fullcat.fits then 
# apply_ZP with FLAGS
# and write ONE file for all CCDs as
#    filemane_fullcat.fits_Obj.csv
#

def apply_ZP_Sexcatalogfitstocsv(catFilename,EXPNUM,CCDNUM,zeropoint,zeropoint_rms,ZPFLAG,outdir,dir):
    import fitsio
    import string,math,csv
    
    outFile="""%s_Obj.csv""" % (catFilename)   
    outFile = outdir+'/'+outFile.split('/')[-1] 
    extension=2
    
    col=['EXPNUM','CCDNUM','NUMBER','ALPHAWIN_J2000','DELTAWIN_J2000','FLUX_AUTO','FLUXERR_AUTO','FLUX_PSF','FLUXERR_PSF','MAG_AUTO','MAGERR_AUTO','MAG_PSF','MAGERR_PSF','SPREAD_MODEL','SPREADERR_MODEL','FWHM_WORLD','FWHMPSF_IMAGE','FWHMPSF_WORLD','CLASS_STAR','FLAGS','IMAFLAGS_ISO','ZeroPoint','ZeroPoint_rms','ZeroPoint_FLAGS']

    hdr=['NUMBER','ALPHAWIN_J2000','DELTAWIN_J2000','FLUX_AUTO','FLUXERR_AUTO','FLUX_PSF','FLUXERR_PSF','MAG_AUTO','MAGERR_AUTO','MAG_PSF','MAGERR_PSF','SPREAD_MODEL','SPREADERR_MODEL','FWHM_WORLD','FWHMPSF_IMAGE','FWHMPSF_WORLD','CLASS_STAR','FLAGS','IMAFLAGS_ISO']
    data = fitsio.read(catFilename,  columns=hdr, ext=extension)[:]
    data = data[np.argsort(data['ALPHAWIN_J2000'])]

    w1=( data['FLUX_AUTO'] >0. )
    data['MAG_AUTO'] = np.where(w1 , (-2.5*np.log10(data['FLUX_AUTO']) - zeropoint) ,np.int16(-9999))
    data['MAGERR_AUTO'] = np.where(w1 ,(2.5/math.log(10.))*(data['FLUXERR_AUTO']/data['FLUX_AUTO']) ,np.int16(-9999))
    w1=( data['FLUX_PSF'] >0. )
    data['MAG_PSF']= np.where(w1 , (-2.5*np.log10(data['FLUX_PSF']) - zeropoint )   ,np.int16(-9999)) 
    data['MAGERR_PSF'] =  np.where(w1 ,(2.5/math.log(10.))*(data['FLUXERR_PSF']/data['FLUX_PSF']) ,np.int16(-9999))

    with open(outFile,'w') as csvFile:
            writer = csv.writer(csvFile,delimiter=',',  quotechar='|',
                                lineterminator='\n', quoting=csv.QUOTE_MINIMAL)

            writer.writerow(col)

            line=[]
            for i in range(data.size):
                line=EXPNUM,CCDNUM,data['NUMBER'][i],data['ALPHAWIN_J2000'][i],data['DELTAWIN_J2000'][i],data['FLUX_AUTO'][i],data['FLUXERR_AUTO'][i],data['FLUX_PSF'][i],data['FLUXERR_PSF'][i],data['MAG_AUTO'][i],data['MAGERR_AUTO'][i],data['MAG_PSF'][i],data['MAGERR_PSF'][i],data['SPREAD_MODEL'][i],data['SPREADERR_MODEL'][i],data['FWHM_WORLD'][i],data['FWHMPSF_IMAGE'][i],data['FWHMPSF_WORLD'][i],data['CLASS_STAR'][i],data['FLAGS'][i],data['IMAFLAGS_ISO'][i],zeropoint,zeropoint_rms,ZPFLAG

                writer.writerow(line)

    data=[]

#############################################
#To Do
#Still needs new args for type of output csv/fits
def Onefile(args):
    import os,glob
    import pandas as pd
    from astropy.io import fits

    if args.verbose >0 : print args

    catlistFile="""Merg_allZP_D%08d_r%sp%02d.csv""" % (args.expnum,args.reqnum,args.attnum)
    are_you_here(catlistFile)

    fout = """D%08d_r%sp%02d_ZP.csv""" % (args.expnum, args.reqnum, args.attnum)
    fitsout = """D%08d_r%sp%02d_ZP.fits""" % (args.expnum, args.reqnum, args.attnum)

    data = np.genfromtxt(catlistFile, dtype=None, delimiter=',', names=True)
        
    for i in range(data['FILENAME'].size):
        apply_ZP_Sexcatalogfitstocsv(data['FILENAME'][i],
                                     data['EXPNUM'][i], data['CCDNUM'][i],
                                     data['NewZP'][i], data['NewZPrms'][i], data['NewZPFlag'][i],
                                     args.outdir,args.dir)
        
    path = './'
    all_files = glob.glob(os.path.join(path, "*Obj.csv"))     

    big_frame = pd.concat((pd.read_csv(f) for f in all_files)).sort(['ALPHAWIN_J2000'], ascending=True)

    big_frame['ID'] = list(range(len(big_frame['ALPHAWIN_J2000'].index)))
    big_frame['ID'] = 1+big_frame['ID']
    
    Cols = ["ID","EXPNUM","CCDNUM","NUMBER","ALPHAWIN_J2000","DELTAWIN_J2000","FLUX_AUTO","FLUXERR_AUTO","FLUX_PSF","FLUXERR_PSF","MAG_AUTO","MAGERR_AUTO","MAG_PSF","MAGERR_PSF","SPREAD_MODEL","SPREADERR_MODEL","FWHM_WORLD","FWHMPSF_IMAGE","FWHMPSF_WORLD","CLASS_STAR","FLAGS","IMAFLAGS_ISO","ZeroPoint","ZeroPoint_rms","ZeroPoint_FLAGS"]

    big_frame.to_csv(fout, sep=',', columns=Cols,index=False)

    # Later Please ADD new args for args.fits/args.csv  with if one/or and
    # Currently BOTH csv and fits are written to disk with  NO ARGS!
    
    col1   = fits.Column(name='ID', format='J', array=np.array(big_frame['ID']))
    col2   = fits.Column(name='EXPNUM', format='I',array=np.array(big_frame['EXPNUM']))
    col3   = fits.Column(name='CCDNUM', format='I', array=np.array(big_frame['CCDNUM']))
    col4   = fits.Column(name='NUMBER', format='I', array=np.array(big_frame['NUMBER']))
    col5   = fits.Column(name='ALPHAWIN_J2000', format='D', array=np.array(big_frame['ALPHAWIN_J2000']))
    col6   = fits.Column(name='DELTAWIN_J2000', format='D', array=np.array(big_frame['DELTAWIN_J2000']))
    col7   = fits.Column(name='FLUX_AUTO', format='D', array=np.array(big_frame['FLUX_AUTO']))
    col8   = fits.Column(name='FLUXERR_AUTO', format='D', array=np.array(big_frame['FLUXERR_AUTO']))
    col9   = fits.Column(name='FLUX_PSF', format='D', array=np.array(big_frame['FLUX_PSF']))
    col10  = fits.Column(name='FLUXERR_PSF', format='D', array=np.array(big_frame['FLUXERR_PSF']))
    col11  = fits.Column(name='MAG_AUTO', format='D', array=np.array(big_frame['MAG_AUTO']))
    col12  = fits.Column(name='MAGERR_AUTO', format='D', array=np.array(big_frame['MAGERR_AUTO'])) 
    col13  = fits.Column(name='MAG_PSF', format='D', array=np.array(big_frame['MAG_PSF']))
    col14  = fits.Column(name='MAGERR_PSF', format='D', array=np.array(big_frame['MAGERR_PSF']))
    col15  = fits.Column(name='SPREAD_MODEL', format='D', array=np.array(big_frame['SPREAD_MODEL']))
    col16  = fits.Column(name='SPREADERR_MODEL', format='D', array=np.array(big_frame['SPREADERR_MODEL']))
    col17  = fits.Column(name='FWHM_WORLD', format='D', array=np.array(big_frame['FWHM_WORLD']))
    col18  = fits.Column(name='FWHMPSF_IMAGE', format='D', array=np.array(big_frame['FWHMPSF_IMAGE'])) 
    col19  = fits.Column(name='FWHMPSF_WORLD', format='D', array=np.array(big_frame['FWHMPSF_WORLD'])) 
    col20  = fits.Column(name='CLASS_STAR', format='D', array=np.array(big_frame['CLASS_STAR']))
    col21  = fits.Column(name='FLAGS', format='I', array=np.array(big_frame['FLAGS'])) 
    col22  = fits.Column(name='IMAFLAGS_ISO', format='I', array=np.array(big_frame['IMAFLAGS_ISO']))
    col23  = fits.Column(name='ZeroPoint', format='D', array=np.array(big_frame['ZeroPoint']))
    col24  = fits.Column(name='ZeroPoint_rms', format='D', array=np.array(big_frame['ZeroPoint_rms']))
    col25  = fits.Column(name='ZeroPoint_FLAGS', format='I', array=np.array(big_frame['ZeroPoint_FLAGS']))

    cols = fits.ColDefs([col1, col2, col3, col4, col5, col6, col7, col8, col9, col10, col11, col12, col13, col14, col15, col16, col17, col18, col19, col20, col21, col22, col23, col24,col25])

    tbhdu = fits.BinTableHDU.from_columns(cols)
    tbhdu.writeto(fitsout)

##################################
#
#
####################################
#knn function gets the dataset and calculates K-Nearest neighbors and distances
def knn(df,k):
    from sklearn.neighbors import NearestNeighbors
    
    nbrs = NearestNeighbors(n_neighbors=3)
    nbrs.fit(df)
    distances, indices = nbrs.kneighbors(df)
    return distances, indices

##################
#reachDist calculates the reach distance of each point to MinPts around it
def reachDist(df, MinPts, knnDist):
    from sklearn.neighbors import NearestNeighbors
    
    nbrs = NearestNeighbors(n_neighbors=MinPts)
    nbrs.fit(df)
    distancesMinPts, indicesMinPts = nbrs.kneighbors(df)
    distancesMinPts[:, 0] = np.amax(distancesMinPts, axis=1)
    distancesMinPts[:, 1] = np.amax(distancesMinPts, axis=1)
    distancesMinPts[:, 2] = np.amax(distancesMinPts, axis=1)
    return distancesMinPts, indicesMinPts

##################
#lrd calculates the Local Reachability Density
def lrd(MinPts,knnDistMinPts):
    return (MinPts / np.sum(knnDistMinPts, axis=1))

##################
#Finally lof calculates lot outlier scores
def lof(Ird, MinPts, dsts):
    lof = []
    for item in dsts:
       tempIrd = np.divide(Ird[item[1:]], Ird[item[0]])
       lof.append(tempIrd.sum()/MinPts)
    return lof

##################
#We flag anything with outlier score greater than 1.2 as outlier
#This is just for charting purposes
def returnFlag(x):
    if x['Score'] > 1.2:
       return 1
    else:
       return 0
    
##################################

if __name__ == "__main__":
    main()

##################################
