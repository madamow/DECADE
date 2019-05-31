import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import scipy.stats as ss
import os, re
import fnmatch
from scipy.spatial import cKDTree
import matplotlib
import astropy.io.fits as pyfits
plt.switch_backend('qt5Agg')

# punkty_zielone = array(n,2)
# punkty_niebieskie = array(m, 2)

A = np.sqrt(np.pi/2.)

def fits_to_csv_cp(cat='./', pattern=''):
    flist = find_files(cat, pattern=pattern)
    dict = {}


    for file in flist:
        f_to_list = file.split("/")[-1].split("_")
        out_file = "_".join(f_to_list[:5])
        f_id = f_to_list[1]+"_"+f_to_list[2]+"_"+f_to_list[3]
        if f_id in dict:
            dict[f_id].append(file)
        else:
            dict[f_id ] =[file]

    for fid in dict:
        df = pd.DataFrame()
        outname = "cp_"+fid+".csv"
        print outname
        for fits_file in dict[fid]:
            ccdnum = fits_file.split("/")[-1].split("_")[6][1:]

            keys = ['ALPHAWIN_J2000', 'DELTAWIN_J2000', 'MAG_PSF', 'MAGERR_PSF', 'SPREAD_MODEL']
            f = pyfits.open(fits_file)
            temp_df = pd.DataFrame()

            for key in keys:
                temp_df[key] = f[2].data.field(key)
            f.close()
            temp_df['ccdnum']=ccdnum
            df = df.append(temp_df)
        df = df.sort_values('ccdnum')
        df = df.reset_index(drop=True)
        df.to_csv(outname)


def find_neighbour_with_tree(t1, t2):
    tree_green = cKDTree(t1)
    tree_blue = cKDTree(t2)
    result = tree_green.query_ball_tree(tree_blue, 0.00001)

    final = [
        (i, o[0]) for i, o in enumerate(result) if o
    ]

    return final

def find_files(cat, pattern=''):
    files = []

    for root, dirnames, filenames in os.walk(cat):
        for filename in fnmatch.filter(filenames, pattern):
            files.append(os.path.join(root, filename))
    return files


def read_cat(cat_file, add_fields=False):
    # read cat and remove bad values for magnitudes
    c = pd.read_csv(cat_file)

    c = c.loc[(c['MAG_PSF'] > -10.) & (np.abs(c['SPREAD_MODEL']) < 0.003) & (c['MAG_PSF'] < 99.)]
    c = c.reset_index()

    if add_fields:
        c['MAG_PSF_F2'] = -999.
        c['MAGERR_PSF_F2'] = -999.
        c['ZeroPoint_F2'] = -999.
        c['ZeroPoint_rms_F2'] = -999.

    return c


def get_comp(c_ref, f2, outname='out.csv'):
    c2 = read_cat(f2)

    #Get two numpy tables with coordinates
    a1 = np.column_stack((c_ref['ALPHAWIN_J2000'],c_ref['DELTAWIN_J2000']))
    a2 = np.column_stack((c2['ALPHAWIN_J2000'], c2['DELTAWIN_J2000']))

    pairs_ind = find_neighbour_with_tree(a1, a2)

    for i,j in pairs_ind:
        c_ref.at[i, 'MAG_PSF_F2'] = c2.loc[j, 'MAG_PSF']
        c_ref.at[i, 'MAGERR_PSF_F2'] = c2.loc[j, 'MAGERR_PSF']
        try:
            c_ref.at[i, 'ZeroPoint_F2'] = c2.loc[j, 'ZeroPoint']
            c_ref.at[i, 'ZeroPoint_rms_F2'] = c2.loc[j, 'ZeroPoint_rms']
        except KeyError:
            pass

    out = c_ref.loc[(c_ref['MAG_PSF_F2'] > 0.) & (c_ref['MAGERR_PSF_F2'] > 0.)]
    out = out.reset_index()
    out.to_csv(outname)


def get_bins(file_in):
    f = pd.read_csv(file_in)
    bin_no = 20
    print file_in
    dict = {}

    f['err'] = np.sqrt(f['MAGERR_PSF'] ** 2 + f['MAGERR_PSF_F2'] ** 2)
    f['dmag'] = A * np.abs(f['MAG_PSF'] - f['MAG_PSF_F2'])

    try:
        f['total_err'] = np.sqrt(f['MAGERR_PSF'] ** 2 + f['MAGERR_PSF_F2'] ** 2 +
                             f['ZeroPoint_rms'] ** 2 + f['ZeroPoint_rms_F2'] ** 2)
    except KeyError:
        f['total_err'] = np.sqrt(f['MAGERR_PSF'] ** 2 + f['MAGERR_PSF_F2'] ** 2)

    bin_err = ss.binned_statistic(f['MAG_PSF'], f['err'], statistic='mean', bins=bin_no)
    bin_dmag = ss.binned_statistic(f['MAG_PSF'], f['dmag'], statistic='mean', bins=bin_no)
    bin_terr = ss.binned_statistic(f['MAG_PSF'], f['total_err'], statistic='mean', bins=bin_no)
    #bin_zp = ss.binned_statistic(f['MAG_PSF'], zp_rms, statistic='mean', bins=bin_no)

    dict['mag_bin'] = bin_err[1][:-1] + (bin_err[1][2] - bin_err[1][1]) * 0.5  # mag value in the center of bin
    dict['err_bin'] = bin_err[0]
    dict['dmag_bin'] = bin_dmag[0]
    dict['t_err'] = bin_terr[0]
    return dict




def make_plot(df):
    colors = ['r','g','b']
    symbols = ['o', '^', 's']
    source = (df['source'].unique()).tolist()

    for froot, grp in df.groupby(by='comp'):

        fig, ax = plt.subplots()
        ax.set_title(froot)
        ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
        ax.get_yaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
        ax.set_yticks([0.01, 0.02, 0.03, 0.05, 0.1, 1])
        ax.set_xlabel(r'$\sqrt{\sigma_1^2 +\sigma_2^2}$')
        ax.set_ylabel(r'$\langle\mid m_1-m_2 \mid\rangle$')
        out_name = froot + "_p1.pdf"
        txt = ''

        fig2, ax2 = plt.subplots()
        ax2.set_title(froot)
        ax2.set_xlabel(r'mag')
       # ax2.set_ylim(0,5)
       # ax2.set_xlim(14, 25)
        ax2.set_ylabel(r'$\langle\mid m_1 - m_2 \mid\rangle / err $')
        ax2.text(0.05, 0.9, r'err$ = \sqrt{\sigma_1^2 +\sigma_2^2 + ZP_{rms,1}^2 + ZP_{rms,2}^2}$',
                 transform=ax2.transAxes)
        out_name2 = froot + "_p2.pdf"


        for i, row in grp.iterrows():
            idx = source.index(row['source'])
            bins = get_bins(row['filename'])
            ax.set_title(froot)
            ab = np.polyfit(bins['err_bin'], bins['dmag_bin'], 1)

    #       FIG 1 linear fit to binned errors and delta mags
            ax.loglog(bins['err_bin'], bins['dmag_bin'], '.', color=colors[idx])
            ax.loglog(bins['err_bin'], np.polyval(ab, bins['err_bin']), color=colors[idx], label=row['source'])
            txt += '%s slope= %6.4f, intercept=%6.4f\n' % (row['source'], ab[0], ab[1])

            # FIG 2
            ax2.errorbar(bins['mag_bin'], bins['dmag_bin'] / bins['t_err'], yerr=bins['t_err'],
                         fmt='o', color=colors[idx], label=row['source'])
            ax2.axhline(1, c='m')

        ax.text(0.05, 0.8, txt, transform=ax.transAxes)
        ax.legend(loc=4)
        ax2.legend()

        fig.savefig(out_name, facecolor="white", bbox_inches='tight')


        fig2.savefig(out_name2, facecolor="white", bbox_inches='tight')
        plt.clf()
        plt.close()


    return ab[0], ab[1], 0# np.mean(zp_rms)


def compare_catalogs(cat='./', pattern='', fits=False, refexp=''):
    flist = find_files(cat, pattern=pattern)

    ref = [s for s in flist if refexp in s][0]

    flist.remove(ref)
    c1 = read_cat(ref, add_fields=True)

    for f2 in flist:
        print f2
        no = int(re.findall(r'\d+', f2.split("/")[-1].split("_")[0])[0])
        out_name = refexp + "_" + str(no) + "_"+ f2.split("/")[1]+".csv"
        print out_name
        try:
            get_comp(c1, f2, outname=out_name)
        except ValueError:
            print "Fail"
            pass


def summary_plot(tab, tname, labels):
    # FIG 3 slopes and ZP for all epochs
    fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True)
    ax1.set_title("REF: " + tname)
    ax1.plot(tab['slope'], '.', label='slopes', color='g')
    ax1.legend(loc=4)
    ax2.plot(tab['intr'], 'x', label='intercepts')
    ax2.plot(tab['zp'], '+', label='$\sqrt{ZP_{rms,1}^2 + ZP_{rms,2}^2}$')
   # ax2.set_ylim(0.001, 0.033)
    ax2.legend(loc=1)
    plt.subplots_adjust(wspace=0, hspace=0)
    xs = range(tab.shape[0])
    plt.xticks(xs, labels, rotation='vertical')
    fig.savefig('summary.pdf', facecolor="white", bbox_inches='tight')
    plt.close()


def get_plots():
    flist = find_files('./catalog_comparison', pattern='*HITS*.csv')
    # Build dataframe with all objects
    df = pd.DataFrame(columns=['filename', 'ref', 'comp', 'source'])
    for i, item in enumerate(flist):
        f = item.split("/")[-1]
        ref = f.split("_")[0]
        comp = f.split("_")[1]
        source = f.split("_")[2].split(".csv")[0]
        df.loc[i] = [item, ref, comp, source]

    df = make_plot(df)
#    columns = ['file', 'slope', 'intr', 'zp']
#    res = pd.DataFrame(columns=columns)
#    labels = []

   # for item in flist:
   #     a, b, zp = make_plot(item)
   #     new = pd.DataFrame([(item, a, b, zp)], columns=columns)
   #     res = res.append(new, ignore_index=True)
   #     labels.append(item.split('/')[-1][10:19])

    #ref_star = flist[0].split('/')[-1][0:9]
    #summary_plot(res, ref_star, labels)
    #plt.show()



def main():
 #   fits_to_csv_cp(cat='./HITS_CP', pattern='*g*.fits')
  #  compare_catalogs(cat='./HITScp', pattern='ccp*_g.csv', refexp='410936')
    get_plots()


if __name__ == "__main__":
    main()
