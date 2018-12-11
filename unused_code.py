# New plots ONLY if NEEDED!!
def plotradec_sexvsY2Q1(args):
    import matplotlib.pyplot as plt

    if args.verbose > 0: print args

    catlistFile = """D%08d_r%sp%02d_red_catlist.csv""" % (args.expnum, args.reqnum, args.attnum)

    are_you_here(catlistFile):
    print '%s does not seem to exist... exiting now...' % catlistFile
    sys.exit(1)


data = np.genfromtxt(catlistFile, dtype=None, delimiter=',', names=True)

for i in range(data['FILENAME'].size):
    ra1 = [];
    ra2 = [];
    dec1 = [];
    dec2 = []
    catFilename = os.path.basename(data['FILENAME'][i])
    rac = data['RA_CENT'][i];
    decc = data['DEC_CENT'][i]
    rac1 = data['RAC1'][i];
    decc1 = data['DECC1'][i]
    rac2 = data['RAC2'][i];
    decc2 = data['DECC2'][i]
    rac3 = data['RAC3'][i];
    decc3 = data['DECC3'][i]
    rac4 = data['RAC4'][i];
    decc4 = data['DECC4'][i]
    CCDpoints = [[rac2, decc2], [rac, decc2], [rac3, decc3], [rac3, decc], [rac4, decc4], [rac, decc4], [rac1, decc1],
                 [rac1, decc]]
    ccdline = plt.Polygon(CCDpoints, fill=None, edgecolor='g')

    pnglistout = """%s.png""" % (catFilename)
    objlistFile = """%s_Obj.csv""" % (catFilename)
    stdlistFile = """%s_std.csv""" % (catFilename)

    are_you_here(objlistFile)
    are_you_here(stdlistFile)

    # Read in the file...
    ra1 = np.genfromtxt(stdlistFile, dtype=float, delimiter=',', skiprows=1, usecols=(1))
    dec1 = np.genfromtxt(stdlistFile, dtype=float, delimiter=',', skiprows=1, usecols=(2))
    ra2 = np.genfromtxt(objlistFile, dtype=float, delimiter=',', skiprows=1, usecols=(1))
    dec2 = np.genfromtxt(objlistFile, dtype=float, delimiter=',', skiprows=1, usecols=(2))
    plt.axes()
    plt.gca().add_patch(ccdline)
    plt.scatter(ra1, dec1, marker='.')
    plt.scatter(ra2, dec2, c='r', marker='+')
    line = plt.Polygon(CCDpoints, fill=None, edgecolor='r')

    plt.title(catFilename, color='#afeeee')
    plt.savefig(pnglistout, format="png")

    ra1 = [];
    ra2 = [];
    dec1 = [];
    dec2 = []

    plt.clf()

#     That was originally in ZP_outliers()
#     #########################################################
#     ######### Extra ONLY if the Full Exposure have problem!
#     #########################################################
#     #FIND OUTLIER via Local Outlier Factor(LOF) see
#     #Ref:
#     #  http://www.dbs.ifi.lmu.de/Publikationen/Papers/LOF.pdf
#     #  THIS is a TEST
#     #########################################################
#     #try
#     #mask = ( df1['NewZPFlag3'] < 0 )
#     #if ( df1[mask]['x'].size > 0 ):
#         #data = df1[['x','y','DiffZP1']]
#         ##You can change below value for different MinPts
#         ##m=5,10,15,30,35,40,45,50,55 testing
#         #m=10
#
#         #knndist, knnindices = knn(data,3)
#         #reachdist, reachindices = reachDist(data,m,knndist)
#
#         #irdMatrix = lrd(m,reachdist)
#         #lofScores = lof(irdMatrix,m,reachindices)
#         #scores= pd.DataFrame(lofScores,columns=['Score'])
#         #mergedData=pd.merge(data,scores,left_index=True,right_index=True)
#         #mergedData['flag'] = mergedData.apply(returnFlag,axis=1)
#         #Outliers = mergedData[(mergedData['flag']==1)]
#         #Normals = mergedData[(mergedData['flag']==0)]
#
#         #mergedData1=pd.concat([df1, mergedData], axis=1)
#         #mergedData1.to_csv(fout,sep=',',index=False)
#
#         #print Outliers
#
#         #pnglistout0="""%s_ZP_Warning.png""" % (catlistFile)
#         #pnglistout1="""%s_ZP_Score.png""" % (catlistFile)
#         ################
#         # New plot in X,Y
#         # NEED to convert the RA, DEC vs the expCal Zero-point mag
#         # New plot SCORE - histogram
#         ################
#
#
#         #l1=plt.scatter(Normals['x'],Normals['y'],Normals['DiffZP1'],c='b',marker='o')
#         #l2=plt.scatter(Outliers['x'],Outliers['y'],Outliers['DiffZP1'],c='r',marker='*')
#         #plt.legend((l1,l2),('Regular','Outlier'),scatterpoints=1,loc='upper left',ncol=1, fontsize=9)
#         #plt.title('Warning D%08d_r%sp%02d %s-Band' %(args.expnum,args.reqnum,args.attnum,BAND))
#         #plt.xlabel(r"$X$", size=18)
#         #plt.ylabel(r"$Y$", size=18)
#         #plt.xlim([min(data['x']),max(data['x'])])
#         #plt.ylim([min(data['y']),max(data['y'])])
#         #plt.savefig(pnglistout0, format="png" )
#         #plt.clf()
#
#         #SCORE - histogram
#         #plt.hist(mergedData['Score'],bins=100,facecolor='red')
#         #plt.xlabel('LOF Score')
#         #plt.ylabel('Frequency')
#         #plt.title('Outlier Scores')
#         #plt.savefig(pnglistout1, format="png" )
#         #plt.clf()

##################################
# New plots for Zero-Points
def plotradec_ZP(args):
    import numpy as np
    import matplotlib.pyplot as plt
    import matplotlib as mpl
    import sys

    if args.verbose > 0: print args

    catlistFile = """D%08d_r%sp%02d_red_catlist.csv""" % (args.expnum, args.reqnum, args.attnum)
    are_you_here(catlistFile)

    # ZeroListFile="""Zero_D%08d_r%sp%02d.csv""" % (args.expnum,args.reqnum,args.attnum)
    # if not os.path.isfile(catlistFile):
    #    print '%s does not seem to exist... exiting now...' % ZeroListFile
    #    sys.exit(1)
    # Mergedout="""Merged_D%08d_r%sp%02d.csv""" % (args.expnum,args.reqnum,args.attnum)
    # print catlistFile,ZeroListFile,Mergedout
    # jointwocsv(catlistFile,ZeroListFile,Mergedout)
    # MergedFile="""Merg_allZP_D%08d_r%sp%02d.csv""" % (args.expnum,args.reqnum,args.attnum)

    MergedFile = """Merged_D%08d_r%sp%02d.csv""" % (args.expnum, args.reqnum, args.attnum)
    data = np.genfromtxt(MergedFile, dtype=None, delimiter=',', names=True)
    stdRA = np.std(data['RA_CENT'])
    if stdRA > 20:
        data['RA_CENT'] = [roundra(x) for x in data['RA_CENT']]
        data['RAC1'] = [roundra(x) for x in data['RAC1']]
        data['RAC2'] = [roundra(x) for x in data['RAC2']]
        data['RAC3'] = [roundra(x) for x in data['RAC3']]
        data['RAC4'] = [roundra(x) for x in data['RAC4']]

    #    print " Unexpected value of the RA spred stdRA=%f \n" % stdRA
    #    sys.exit(1)
    w0 = (data['ZP'] == -999)
    w1 = (data['ZP'] > -999)

    data0 = data[w0]
    data1 = data[w1]

    if (len(data1) == 0):
        sys.exit(1)

    BAND = data1['BAND'][0]
    zpmedian = np.median(data1['ZP'])

    pnglistout0 = """%s_ZP.png""" % (catlistFile)
    pnglistout1 = """%s_deltaZP.png""" % (catlistFile)
    pnglistout2 = """%s_NumClipstar.png""" % (catlistFile)
    pnglistout3 = """%s_CCDsvsZPs.png""" % (catlistFile)

    w2 = (data['NewZPFlag'] == 0)
    w3 = (data['NewZPFlag'] == 1)

    # w3=(data('NewZPFlag') >0 )
    data2 = data[w2]
    data3 = data[w3]
    pnglistout4 = """%s_NewZP.png""" % (catlistFile)
    pnglistout5 = """%s_NewdeltaZP.png""" % (catlistFile)

    minra = min(min(data['RA_CENT']), min(data['RAC1']), min(data['RAC2']), min(data['RAC3']), min(data['RAC4'])) - .075
    mindec = min(min(data['DEC_CENT']), min(data['DECC1']), min(data['DECC2']), min(data['DECC3']),
                 min(data['DECC4'])) - .075
    maxra = max(max(data['RA_CENT']), max(data['RAC1']), max(data['RAC2']), max(data['RAC3']), max(data['RAC4'])) + .075
    maxdec = max(max(data['DEC_CENT']), max(data['DECC1']), max(data['DECC2']), max(data['DECC3']),
                 max(data['DECC4'])) + .075

    ################
    # New plot the RA, DEC vs the expCal Zero-point mag
    ################
    l1 = plt.scatter(data0['RA_CENT'], data0['DEC_CENT'], c=data0['ZP'], s=15, marker=(25, 0), cmap=mpl.cm.spectral,
                     vmin=np.min(data1['ZP']), vmax=np.max(data1['ZP']))
    l2 = plt.scatter(data1['RA_CENT'], data1['DEC_CENT'], c=data1['ZP'], s=500, marker=(5, 0), cmap=mpl.cm.spectral,
                     vmin=np.min(data1['ZP']), vmax=np.max(data1['ZP']))
    cbar = plt.colorbar(ticks=np.linspace(np.min(data1['ZP']), np.max(data1['ZP']), 4))
    cbar.set_label('Zero-Point Mag')
    # plt.legend((l1,l2),('No Y2Q1 data','ExpCal'),scatterpoints=1,loc='upper left',ncol=1, fontsize=9)
    plt.legend((l1, l2), ('No APASS data', 'ExpCal'), scatterpoints=1, loc='upper left', ncol=1, fontsize=9)

    for i in range(data['RA_CENT'].size):
        CCDpoints = [[data['RAC2'][i], data['DECC2'][i]], [data['RAC3'][i], data['DECC3'][i]],
                     [data['RAC4'][i], data['DECC4'][i]], [data['RAC1'][i], data['DECC1'][i]]]
        ccdline = plt.Polygon(CCDpoints, fill=None, edgecolor='k')
        plt.gca().add_patch(ccdline)

    plt.title('D%08d_r%sp%02d %s-Band ZP_Median=%.3f ' % (args.expnum, args.reqnum, args.attnum, BAND, zpmedian))
    plt.xlabel(r"$RA$", size=18)
    plt.ylabel(r"$DEC$", size=18)
    plt.xlim([minra, maxra])
    plt.ylim([mindec, maxdec])
    plt.savefig(pnglistout0, format="png")
    plt.clf()

    ################
    # New plot the RA, DEC vs the expCal Delta Zero-point mag from median
    ################

    l1 = plt.scatter(data0['RA_CENT'], data0['DEC_CENT'], c=data0['ZP'], s=15, marker=(25, 0), cmap=mpl.cm.spectral,
                     vmin=np.min(data1['ZP']), vmax=np.max(data1['ZP']))
    l2 = plt.scatter(data1['RA_CENT'], data1['DEC_CENT'], c=1000 * (data1['ZP'] - zpmedian), s=500, marker=(5, 0),
                     cmap=mpl.cm.spectral, vmin=min(1000 * (data1['ZP'] - zpmedian)),
                     vmax=max(1000 * (data1['ZP'] - zpmedian)))
    cbar = plt.colorbar(
        ticks=np.linspace(min(1000 * (data1['ZP'] - zpmedian)), max(1000 * (data1['ZP'] - zpmedian)), 4))
    cbar.set_label('Delta Zero-Point mili-Mag')
    plt.legend((l1, l2), ('No APASS data', 'ExpCal'), scatterpoints=1, loc='upper left', ncol=1, fontsize=9)

    for i in range(data['RA_CENT'].size):
        CCDpoints = [[data['RAC2'][i], data['DECC2'][i]], [data['RAC3'][i], data['DECC3'][i]],
                     [data['RAC4'][i], data['DECC4'][i]], [data['RAC1'][i], data['DECC1'][i]]]
        ccdline = plt.Polygon(CCDpoints, fill=None, edgecolor='k')
        plt.gca().add_patch(ccdline)

    plt.title('D%08d_r%sp%02d %s-Band ZP_Median=%.3f ' % (args.expnum, args.reqnum, args.attnum, BAND, zpmedian))

    plt.xlabel(r"$RA$", size=18)
    plt.ylabel(r"$DEC$", size=18)
    plt.xlim([minra, maxra])
    plt.ylim([mindec, maxdec])
    plt.savefig(pnglistout1, format="png")
    plt.clf()

    ################
    # New plot RA DEC vs Number of stars clipped stars from expCal
    ################
    l1 = plt.scatter(data0['RA_CENT'], data0['DEC_CENT'], c=data0['Nclipped'], s=15, marker=(25, 0),
                     cmap=mpl.cm.spectral)
    l2 = plt.scatter(data1['RA_CENT'], data1['DEC_CENT'], c=data1['Nclipped'], s=500, marker=(5, 0),
                     cmap=mpl.cm.spectral)

    cbar = plt.colorbar()
    cbar.set_label('No. 3 $\sigma$ clipped Stars')
    plt.legend((l1, l2), ('No APASS data', 'expCal'), scatterpoints=1, loc='upper left', ncol=1, fontsize=9)

    for i in range(data['RA_CENT'].size):
        CCDpoints = [[data['RAC2'][i], data['DECC2'][i]], [data['RAC3'][i], data['DECC3'][i]],
                     [data['RAC4'][i], data['DECC4'][i]], [data['RAC1'][i], data['DECC1'][i]]]
        ccdline = plt.Polygon(CCDpoints, fill=None, edgecolor='k')
        plt.gca().add_patch(ccdline)

    plt.title('D%08d_r%sp%02d %s-Band ZP_Median=%.3f ' % (args.expnum, args.reqnum, args.attnum, BAND, zpmedian))

    plt.xlabel(r"$RA$", size=18)
    plt.ylabel(r"$DEC$", size=18)
    plt.xlim([minra, maxra])
    plt.ylim([mindec, maxdec])
    plt.savefig(pnglistout2, format="png")
    plt.clf()

    ################
    # New plot CCDs vs ZP from expCal
    plt.errorbar(data0['CCDNUM'], data0['ZP'], data0['ZPrms'], fmt='o', label='No APASS data')
    plt.errorbar(data1['CCDNUM'], data1['ZP'], data1['ZPrms'], fmt='o', label='expCal')
    legend = plt.legend(loc='upper center', shadow=None, fontsize=12)
    legend.get_frame().set_facecolor('#00FFCC')
    plt.title('D%08d_r%sp%02d %s-Band ZP_Median=%.3f ' % (args.expnum, args.reqnum, args.attnum, BAND, zpmedian))
    plt.xlabel(r"$CCDs$", size=18)
    plt.ylabel(r"$Zero Points$", size=18)
    plt.ylim(min(data1['ZP']) - .01, max(data1['ZP']) + .02)
    plt.xlim(min(data1['CCDNUM']) - 1.5, max(data1['CCDNUM']) + 1.5)
    plt.savefig(pnglistout3, format="png")
    plt.clf()

    ################
    # New plot the RA, DEC vs the NEW Zero-point mag
    ################
    l1 = plt.scatter(data2['RA_CENT'], data2['DEC_CENT'], c=data2['NewZP'], s=500, marker=(5, 0), cmap=mpl.cm.spectral,
                     vmin=np.min(data2['NewZP']), vmax=np.max(data2['NewZP']))
    l2 = plt.scatter(data3['RA_CENT'], data3['DEC_CENT'], c=data3['NewZP'], s=25, marker=(25, 0), cmap=mpl.cm.spectral,
                     vmin=np.min(data2['NewZP']), vmax=np.max(data2['NewZP']))
    cbar = plt.colorbar(ticks=np.linspace(np.min(data2['NewZP']), np.max(data2['NewZP']), 4))
    cbar.set_label('Zero-Point Mag')
    # CHANGE
    plt.legend((l1, l2), ('CCD', 'allEXP'), scatterpoints=1, loc='upper left', ncol=1, fontsize=9)
    for i in range(data['RA_CENT'].size):
        CCDpoints = [[data['RAC2'][i], data['DECC2'][i]], [data['RAC3'][i], data['DECC3'][i]],
                     [data['RAC4'][i], data['DECC4'][i]], [data['RAC1'][i], data['DECC1'][i]]]
        ccdline = plt.Polygon(CCDpoints, fill=None, edgecolor='k')
        plt.gca().add_patch(ccdline)

    plt.title('D%08d_r%sp%02d %s-Band ' % (args.expnum, args.reqnum, args.attnum, BAND))
    plt.xlabel(r"$RA$", size=18)
    plt.ylabel(r"$DEC$", size=18)
    plt.xlim([minra, maxra])
    plt.ylim([mindec, maxdec])
    plt.savefig(pnglistout4, format="png")
    plt.clf()

    ################
    # New plot the RA, DEC vs the NEW Delta Zero-point mag from median
    ################
    l1 = plt.scatter(data2['RA_CENT'], data2['DEC_CENT'], c=data2['DiffZP1'], s=500, marker=(5, 0),
                     cmap=mpl.cm.spectral, vmin=np.min(data2['DiffZP1']), vmax=np.max(data2['DiffZP1']))
    l2 = plt.scatter(data3['RA_CENT'], data3['DEC_CENT'], c=data3['DiffZP1'], s=25, marker=(25, 0),
                     cmap=mpl.cm.spectral, vmin=min(data2['DiffZP1']), vmax=max(data2['DiffZP1']))
    cbar = plt.colorbar(ticks=np.linspace(min(data2['DiffZP1']), max(data2['DiffZP1']), 4))
    cbar.set_label('Delta Zero-Point mili-Mag')
    plt.legend((l1, l2), ('CDD', 'allExP'), scatterpoints=1, loc='upper left', ncol=1, fontsize=9)

    for i in range(data['RA_CENT'].size):
        CCDpoints = [[data['RAC2'][i], data['DECC2'][i]], [data['RAC3'][i], data['DECC3'][i]],
                     [data['RAC4'][i], data['DECC4'][i]], [data['RAC1'][i], data['DECC1'][i]]]
        ccdline = plt.Polygon(CCDpoints, fill=None, edgecolor='k')
        plt.gca().add_patch(ccdline)

    plt.title('D%08d_r%sp%02d %s-Band ' % (args.expnum, args.reqnum, args.attnum, BAND))

    plt.xlabel(r"$RA$", size=18)
    plt.ylabel(r"$DEC$", size=18)
    plt.xlim([minra, maxra])
    plt.ylim([mindec, maxdec])
    plt.savefig(pnglistout5, format="png")
    plt.clf()
    #