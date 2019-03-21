#!/usr/bin/env python
# $Id: tag_firstpass.py 39826 2015-09-03 22:22:22Z mgower $
# $Rev:: 39826                            $:  # Revision of last commit.
# $LastChangedBy:: mgower                 $:  # Author of last commit.
# $LastChangedDate:: 2015-09-03 17:22:22 #$:  # Date of last commit.

""" Examine a set of exposures to make a first pass decision about 
    whether it should be included in a processing tag.
"""

import argparse
import despydb.desdbi as desdbi
from datetime import datetime, timedelta
import csv
import sys
import os
import glob
import re
import time
import math

verbose=0
updateDB=False

#################################
def CSVtoDict(filename,verbose=0):
    """ Read CSV file and instantiate as a dictionary """

    Dict={}
    t0=time.time()

    if (os.path.isfile(filename)):
        if (verbose > 0):
            print("Reading CSV file: {:s} and rendering as a dict of dicts".format(filename))
        csvfile=open(filename,'r')
        csvreader=csv.DictReader(csvfile)
 
        for row in csvreader:
            fname=row['FILENAME'].split("/")[-1]
            row['FILENAME']=fname
            for d in row:
                if (row[d] in ["","-9999.0"]):
                    row[d]=None
            Dict[fname]=row  
        csvfile.close()
    else:
        print("File {:s} does not exist".format(filename))
        print("Aborting!")
        exit(1)
    if (verbose > 0):
        print("Read {:d} records.".format(len(Dict)))


    return Dict


#################################
def GetImgData(Dict,dbh,dbSchema,verbose=0):
    """ Obtain image data and update Dict based on catalog files (keys) in Dict"""

#
#   Prepare GTT_FILENAME table with list of possible inputs 
#
    CatList=[]
    for key in Dict:
        CatList.append([key.split("/")[-1]])

    # Make sure the GTT_FILENAME table is empty
    curDB=dbh.cursor()
    curDB.execute('delete from GTT_FILENAME')
    # load img ids into opm_filename_gtt table
    print("# Loading GTT_FILENAME table for secondary queries with entries for {:d} images".format(len(CatList)))
    dbh.insert_many('GTT_FILENAME',['FILENAME'],CatList)
#
#   Obtain associated images (red_immask).
#
    query="""SELECT 
        i.filename as imagefile,
        c.filename as catalogfile,
        i.pfw_attempt_id as pfw_attempt_id,
        i.band as band,
        i.expnum as expnum,
        i.ccdnum as ccdnum
    FROM image i, catalog c, GTT_FILENAME gtt
    WHERE c.filename=gtt.filename
        and i.pfw_attempt_id=c.pfw_attempt_id
        and i.filetype='red_immask'
        and i.ccdnum=c.ccdnum
    """.format(schema=dbSchema)

    if (verbose > 0):
        print("# Executing query to obtain red_immask images corresponding to the cat_finalcut catalogs")
        if (verbose == 1):
            print("# sql = {:s} ".format(" ".join([d.strip() for d in query.split('\n')])))
        if (verbose > 1):
            print("# sql = {:s}".format(query))
    t1=time.time()
    curDB.execute(query)
    desc = [d[0].lower() for d in curDB.description]

    for row in curDB:
        rowd = dict(zip(desc, row))
        CatName=rowd['catalogfile']
        ImgName=rowd['imagefile']

        if (CatName in Dict):
            CatRow=Dict[CatName]
            if ((rowd['ccdnum']==int(CatRow['CCDNUM']))and
                (rowd['expnum']==int(CatRow['EXPNUM']))):

                Dict[CatName]['IMAGEFILE']=rowd['imagefile']
                Dict[CatName]['BAND']=rowd['band']
            else:
                print("Catalog mismatch to Image: expnum=({cexp:7s} vs {iexp:7d}), ccdnum=({cccd:2s} vs {iccd:02d}))".format(
                    cexp=Dict[CatName]['EXPNUM'],iexp=rowd['expnum'],
                    cccd=Dict[CatName]['CCDNUM'],iccd=rowd['ccdnum']))
        else:
            print('Warning: Catalog name not in Dict (something is really wrong!!!).  Catname={:s}'.format(CatName))
    t2=time.time()
    print("Timing for IMG query: {:.2f}".format(t2-t1))

    return Dict
    

#################################
def ReplaceImgCatData(Dict,ProcTag,dbh,dbSchema,verbose=0):
    """Use catalog files (keys in Dict) to obtain an alternate set of
       of image and catalog based on a proctag.  A new dictionary is 
       formed with replacement of the original catalog filename.
    """

#
#   Prepare GTT_FILENAME table with list of possible inputs 
#
    CatList=[]
    for key in Dict:
        CatList.append([key[10:]])
    
    # Make sure the GTT_FILENAME table is empty
    curDB=dbh.cursor()
    curDB.execute('delete from GTT_FILENAME')
    # load img ids into opm_filename_gtt table
    print("# Loading GTT_FILENAME table for secondary queries with entries for {:d} images".format(len(CatList)))
    dbh.insert_many('GTT_FILENAME',['FILENAME'],CatList)
#
#   Obtain associated images (red_immask).
#
    query="""SELECT 
        i.filename as imagefile,
        d.filename as catalogfile,
        c.filename as oldcatalogfile,
        i.pfw_attempt_id as pfw_attempt_id,
        i.band as band,
        i.expnum as expnum,
        i.ccdnum as ccdnum
    FROM {schema:s}image i, {schema:s}catalog c, {schema:s}catalog d, {schema:s}proctag t, GTT_FILENAME gtt
    WHERE c.filename=gtt.filename
        and c.expnum=d.expnum
        and c.ccdnum=d.ccdnum
        and d.filetype='cat_finalcut'
        and d.pfw_attempt_id=t.pfw_attempt_id
        and t.tag='{ptag:s}'
        and i.pfw_attempt_id=d.pfw_attempt_id
        and i.filetype='red_immask'
        and i.ccdnum=c.ccdnum
    """.format(schema=dbSchema,ptag=ProcTag)

    if (verbose > 0):
        print("# Executing query to replace catalogfiles with those from PROCTAG={:s} (and obtain associated red_immask images)".format(ProcTag))
        if (verbose == 1):
            print("# sql = {:s} ".format(" ".join([d.strip() for d in query.split('\n')])))
        if (verbose > 1):
            print("# sql = {:s}".format(query))
    t1=time.time()
    curDB.execute(query)
    desc = [d[0].lower() for d in curDB.description]

    NewDict={}
    for row in curDB:
        rowd = dict(zip(desc, row))
        OldCatName=rowd['oldcatalogfile']
        CatName=rowd['catalogfile']
        ImgName=rowd['imagefile']

        if (OldCatName in Dict):
            CatRow=Dict[OldCatName]
            if ((rowd['ccdnum']==int(CatRow['CCDNUM']))and
                (rowd['expnum']==int(CatRow['EXPNUM']))):

                NewDict[CatName]=Dict[OldCatName]
                NewDict[CatName]['FILENAME']=rowd['catalogfile']
                NewDict[CatName]['IMAGEFILE']=rowd['imagefile']
                NewDict[CatName]['BAND']=rowd['band']
                if (verbose > 2):
                    print("Replacing {:s} with {:s} (and linking image to {:s})".format(OldCatName,CatName,ImgName))
            else:
                print("Warning: Catalog mismatch to Image: expnum=({cexp:7s} vs {iexp:7d}), ccdnum=({cccd:2s} vs {iccd:02d}))".format(
                    cexp=Dict[CatName]['EXPNUM'],iexp=rowd['expnum'],
                    cccd=Dict[CatName]['CCDNUM'],iccd=rowd['ccdnum']))
                print("Warning: This really should not be possible")
        else:
            print('Warning: Old Catalog name not in Dict (something is really wrong!!!).  Catname={:s}'.format(OldCatName))

    t2=time.time()
    print("Timing for Replacement CAT/IMG query: {:.2f}".format(t2-t1))

    return NewDict
    

#################################
def ingest_zeropoint(ZPT_Dict,DBtable,DBorder,DtoC,dbh,updateDB,verbose=0):
    """Ingest a set of zeropoints"""
    t0=time.time()
    InsertCnt=0
    print("Preparing lists of list to ingest many zeropoints")
    
#   Preliminary sanity check (Do I know which data comes from where)
    CheckCheckIt=True
    for col in DBorder:
        if (col not in DtoC):
            print('No entry for column ({:s}) in DtoC.  Unable to Proceed!'.format(col))
            CheckCheckIt=False

    if (not(CheckCheckIt)):
        print('Aborting!')
        exit(1)

    new_data=[]
    for entry in ZPT_Dict:
        new_row=[]
        for col in DBorder:
            if (DtoC[col]['src'] == 'CSV'):
                new_row.append(ZPT_Dict[entry][DtoC[col]['col']])
            elif (DtoC[col]['src'] == 'arg'):
                new_row.append(DtoC[col]['val'])
            elif (DtoC[col]['src'] == 'cal'):
                if (col == "MAG_ZERO"):
                    zptval="{:.6f}".format(-1.0*float(ZPT_Dict[entry]['NewZP']))
                    new_row.append(zptval)
                elif (col == "FLAG"):
                    flag_val=0
                    if (int(ZPT_Dict[entry]['NewZPFlag'])>=0):
                        flag_val=int(ZPT_Dict[entry]['NewZPFlag'])
                    else:
#                       The following case has been deprecated (in that a correction
#                       occurs in the main code.                        
                        if (int(ZPT_Dict[entry]['NewZPFlag']) < 0):
                            flag_val=abs(int(ZPT_Dict[entry]['NewZPFlag']))
#                            print("Updated flag for {:s}".format(entry))
                    new_row.append(flag_val)
                else:
                    new_row.append(0)
            else:
                new_row.append(None)
        new_data.append(new_row)
    t1=time.time()

    print("Successfully Formed list for Ingestion (Nrows={:d}, Timing: {:.2f})".format(len(new_data),(t1-t0)))

    if (updateDB):
        print("# Loading {:s} with {:d} entries".format(DBtable,len(new_data)))
        t1=time.time()
        dbh.insert_many(DBtable,DBorder,new_data)
        t2=time.time()
        print("Commit of {:d} rows, timing={:.2f} of {:.2f})".format(len(new_data),(t2-t1),(t2-t0)))
        dbh.commit()
    ningest=len(new_data)        

    return ningest


#################################
def parse_argv(argv):
    """ Return command line arguments as dict """
    parser = argparse.ArgumentParser(description='Examine a set of exposures to determine fitness for a processing campaign tag')
#    parser.add_argument('-e','--exptag', action='store', required=True)
    parser.add_argument('--des_services', action='store')
    parser.add_argument('-s','--des_db_section', action='store', required=True)
    parser.add_argument('-S','--Schema', action='store', required=True)
    parser.add_argument('-u','--updateDB', action='store_true', default=False)
    parser.add_argument('-D','--DBtable', action='store',default='ZEROPOINT')
    parser.add_argument('-c','--common_time', action='store_true', default=False, help='Flag to cause all inserts to share a common timestamp')

    parser.add_argument('--file', action='store', type=str, default=None, required=True, help='File containing ZPT data from PGCM')

    parser.add_argument('--source',   action='store', type=str, required=True, help='Source of zeropoints being ingested (default=GCM)')
    parser.add_argument('--version',  action='store', type=str, required=True, help='Version responsible for these Zeropoints')
    parser.add_argument('--tag',      action='store', type=str, required=True, help='TAG approriatiate for these Zeropoints')

    parser.add_argument('--proctag_replace', action='store', type=str, default=None, help='Replace image and catalog files from a different proctag (using expnum, ccdnum info')
    parser.add_argument('--partial_replace', action='store_true', default=False, help='Flag to allow partial replacement to occur')

    parser.add_argument('-v', '--verbose', action='store', type=str, default=0, help='Verbosity (defualt:0; currently values up to 2)')
    args = vars(parser.parse_args(argv))   # convert dict

    return args


#################################
def main():
    """ Program entry point """
    args = parse_argv(sys.argv[1:])

    if args['verbose'] is not None:
        verbose = int(args['verbose'])
    else:
        verbose=0

    des_services=args['des_services']
    des_db_section=args['des_db_section']
    dbSchema="%s." % (args['Schema'])
    updateDB=args['updateDB']

    DBtable=args['DBtable']
    
    if (args['common_time']):
        common_time=datetime.now()
    else:
        common_time=None

    dbh = desdbi.DesDbi(des_services, des_db_section)
    cur = dbh.cursor()

    ZPT_Dict=CSVtoDict(args['file'],verbose=verbose)

    if (args['proctag_replace'] is None):

#       Find associated Image to the original catalog files

        ZPT_Dict=GetImgData(ZPT_Dict,dbh,dbSchema,verbose=verbose)    
    else:
#
#       This version allows for the replacement of catalogs from a different processing set
#       (useful as a hack when Zpt Info comes from a different processing attempt
#
        orig_len=len(ZPT_Dict)
        print("Attempting to replace catalog(and image files) with versions from PROCTAG={:s}".format(args['proctag_replace']))

        ZPT_Dict=ReplaceImgCatData(ZPT_Dict,args['proctag_replace'],dbh,dbSchema,verbose=verbose)

        new_len=len(ZPT_Dict)
        if (new_len == orig_len):
            print("Size of new dictionary indicates a full replacement has occurred")
        else:
            print("Warning: Original Dictionary had {:d} entries.  New Dictionary has {:d} entries.".format(orig_len,new_len))
            print("Warning: It is almost certain that a full, one-for-one, replacement is not possible")
            if (not(args.partial_replace)):
                print("Warning: In order to make a partial replacement use --partial_replace")
                print("Aborting!")
                exit(1)
            else:
                print("Warning: Proceeding with partial replacement")

#
#   Prepare to translate data into the DB.
#
    DBorder=['IMAGENAME','SOURCE','VERSION','CATALOGNAME','EXPNUM','CCDNUM','BAND',
             'INSERT_DATE','MAG_ZERO','SIGMA_MAG_ZERO','MAG_ONE','SIGMA_MAG_ONE','TAG','FLAG']

    DBtoCSV_ColDict={
        'IMAGENAME':{'src':'CSV','col':'IMAGEFILE'},
        'SOURCE':   {'src':'arg','val':args['source']},
        'VERSION':  {'src':'arg','val':args['version']},
        'CATALOGNAME':{'src':'CSV','col':'FILENAME'},
        'EXPNUM':   {'src':'CSV','col':'EXPNUM'},
        'CCDNUM':   {'src':'CSV','col':'CCDNUM'},
        'BAND':     {'src':'CSV','col':'BAND'},
        'INSERT_DATE':{'src':'arg','val':common_time},
        'MAG_ZERO': {'src':'cal'},
        'SIGMA_MAG_ZERO':{'src':'CSV','col':'NewZPrms'},
        'MAG_ONE':  {'src':'arg','val':None},
        'SIGMA_MAG_ONE':{'src':'arg','val':None},
        'TAG':      {'src':'arg','val':args['tag']},
        'FLAG':     {'src':'cal'}
    }

#
#   Fix any Null/None values in MAG_ZERO_MEAN_ERR'
#   Also, current version of expCalib can produce Negative Flags 
#      which are really MAG_ZEROs that should not be used so need to be positive for COADD pipeline.
#   Therefore, replace negative values for NewZPFlag with their absolute value
#
    FixCnt=0
    FixCntFlag=0
    for key in ZPT_Dict:
        if (ZPT_Dict[key]['NewZPrms'] is None):
            FixCnt=FixCnt+1
            ZPT_Dict[key]['NewZPrms']="-1.0"
            ZPT_Dict[key]['NewZPFlag']='-2000'
        if (int(ZPT_Dict[key]['NewZPFlag']) < 0):
            FixCntFlag=FixCntFlag+1
            ZPT_Dict[key]['NewZPFlag']='{:d}'.format(abs(int(ZPT_Dict[key]['NewZPFlag'])))


    print("Fixed {:d} records with MAG_ZERO_MEAN_ERR are NoneType".format(FixCnt))
    print("Fixed {:d} records with NewZPFlag values are negative".format(FixCntFlag))

#
#   Go on and ingest
#
    ningest=ingest_zeropoint(ZPT_Dict,DBtable,DBorder,DBtoCSV_ColDict,dbh,updateDB,verbose=verbose)
    
    
    print "Finished! Closing up shop."
    dbh.close()


if __name__ == "__main__":
    main()
