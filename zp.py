import numpy as np
import pandas as pd
from despydb import DesDbi
import argparse
from make_red_catlist import MakeRCat
from expCalib_Y3apass import main as expcalib
from ingest_expCalib_zpt import main as ingest
import glob
import os

dbh = DesDbi(None, 'db-decade', retry=True)
cur = dbh.cursor()

##############
"""Create command line arguments"""

parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)

#make_red_catlist arguments:
parser.add_argument('--filein', help='create catalog for this file', default='.', type=str)

# expCalib arguments:
parser.add_argument('--caldir', help='caldir is the calibration directory',
                    default='/des002/devel/emorgan2/APASS_TWOMASS/', type=str)
parser.add_argument('--outdir', help='dir is the production directory', default='./', type=str)
parser.add_argument('--expnum', help='expnum is queried',  type=int)
parser.add_argument('--reqnum', help='reqnum is queried',  type=str)
parser.add_argument('--attnum', help='attnum is queried',  type=int)
parser.add_argument('--magType', help='mag type to use (mag_psf, mag_auto, mag_aper_8, ...)', default='mag_psf')
parser.add_argument('--sex_mag_zeropoint',
                    help='default sextractor zeropoint to use to convert fluxes to sextractor mags \
                         (mag_sex = -2.5log10(flux) + sex_mag_zeropoint)',
                    type=float, default=25.0)
parser.add_argument('--verbose', help='verbosity level of output to screen (0,1,2,...)', default=0, type=int)
parser.add_argument('--debug', help='debugging option', dest='debug', action='store_true', default=False)

# ingest arguments
parser.add_argument('--des_services', action='store')
parser.add_argument('-s','--des_db_section', action='store', required=True)
parser.add_argument('-S','--Schema', action='store', required=True)
parser.add_argument('-u','--updateDB', action='store_true', default=False)
parser.add_argument('-D','--DBtable', action='store',default='ZEROPOINT')
parser.add_argument('-c','--common_time', action='store_true', default=False, help='Flag to cause all inserts to share a common timestamp')
parser.add_argument('--file', action='store', type=str, default=None, help='File containing ZPT data from expCalib: Merged_*')
parser.add_argument('--source',   action='store', type=str, required=True, help='Source of zeropoints being ingested)')
parser.add_argument('--version',  action='store', type=str, required=True, help='Version responsible for these Zeropoints')
parser.add_argument('--tag',      action='store', type=str, default=None, help='TAG approriatiate for these Zeropoints')
parser.add_argument('--proctag_replace', action='store', type=str, default=None, help='Replace image and catalog files from a different proctag (using expnum, ccdnum info')
parser.add_argument('--partial_replace', action='store_true', default=False, help='Flag to allow partial replacement to occur')

# for querying database in particular band
parser.add_argument('--band', help='band', default='all', type=str)


args = parser.parse_args()
args.magType = args.magType.strip().upper()

# Check what is processed
i = 0
query_db = True
todo = pd.DataFrame(columns=['expnum', 'band','path'])

while query_db:
    if args.band == 'all':
        query = "select o.expnum, e.band from OPS_AUTO_QUEUE o, exposure e where processed=1 and o.expnum=e.expnum and not exists(select 1 from zeropoint zp where zp.expnum=o.expnum)"
        query+= " and e.band not in ('VR', 'N964') and e.band is not NULL  order by o.expnum" 
        query+= " offset %i rows  fetch next 1000 rows only" % (i*1000)
    else:
        query = "select o.expnum, e.band from OPS_AUTO_QUEUE o, exposure e where processed=1 and o.expnum=e.expnum and not exists(select 1 from zeropoint zp where zp.expnum=o.expnum)"
        query+= " and  e.band='%s'  order by o.expnum" % (args.band)                                                                                                   
        query+= " offset %i rows  fetch next 1000 rows only" % (i*1000)
    cur.execute(query)
    try:
        temp = pd.DataFrame(cur.fetchall(), columns=['expnum', 'band'])
    except ValueError:
        temp = pd.DataFrame(columns=['expnum', 'band'])
    if temp.shape[0]<1000:
        query_db = False

    # Add paths to files 
    unitnames = ['D00'+str(e) for e  in temp.expnum.values]
    path_query = "select unitname, archive_path from pfw_attempt a, task t,pfw_request r where r.reqnum=a.reqnum "
    path_query += "and t.id=a.task_id and r.project='DEC' and unitname in ('%s')" % ("','".join(unitnames))

    cur.execute(path_query)
    path_df = pd.DataFrame(cur.fetchall(), columns=['unitname', 'path'])
    path_df['expnum']=''
    for u in path_df['unitname']:
        e_no = int(u.split('D00')[1])
        path_df.loc[path_df['unitname']==u, 'expnum']=e_no

    temp_with_path = pd.merge(temp, path_df, on=['expnum'], how='inner')
    
    todo = todo.append(temp_with_path)
    i += 1

dbh.close()

total = todo.shape[0]

if total == 0:
    print "Nothing to do"
    exit()

def cleanup(exposure):
    for f in glob.glob("*%s*.csv" % (exposure)):
        os.remove(f)

for no, row in todo.iterrows():
    print "\n##########\n %i / %i %s band=%s \n##########\n" % (no+1, total, row['unitname'], row['band'])
    try:
        args.filein =  "/deca_archive/"+row['path']
    except TypeError:
        print "Something wrong with path"
        continue
    print args.filein
    args.expnum = int(row['expnum'])
    args.reqnum = row['path'].split("/")[3].split("r")[1][:4]
    args.attnum = int(row['path'].split("/")[-1].split("p")[1])
    
    # Create catalog
    try:
        mrc = MakeRCat(args)
        mrc.runIt()
    except IndexError:
        print "make_red_catalog: list is empty"
        continue
    
    # Run expCalib
    try:
        expcalib(args)
    except IOError:
        print "Some missing file!"
        cleanup(args.expnum)
        continue    
    # Ingest to database
    merged = "Merged_%s_r%sp%02d.csv" % (row['unitname'], args.reqnum, args.attnum)
    args.file = merged
    ingest(args)
        
    # Move Merged_* file to ZPs catalog
    os.rename(merged, "./ZPs/"+merged)
    # Clean up
    cleanup(args.expnum)
 


