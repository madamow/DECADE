import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import scipy.stats as ss
import os
import fnmatch
from scipy.spatial import cKDTree

# punkty_zielone = array(n,2)
# punkty_niebieskie = array(m, 2)

A = np.sqrt(np.pi/2.)


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
    c = c.loc[(c['MAG_PSF'] > -10.) & (np.abs(c['SPREAD_MODEL']) < 0.003)]

    if add_fields:
        c['MAG_PSF_F2'] = -999.
        c['MAGERR_PSF_F2'] = -999.
        c['ZeroPoint_F2'] = -999.
        c['ZeroPoint_rms_F2'] = -999.

    return new


def get_comp(c_ref, f2, outname='out.csv'):
    c2 = read_cat(f2)

    for i, r1 in c_ref.iterrows():
        close = c2.loc[(np.abs(r1['ALPHAWIN_J2000'] - c2['ALPHAWIN_J2000']) < 0.00001) &
                       (np.abs(r1['DELTAWIN_J2000'] - c2['DELTAWIN_J2000']) < 0.00001)]
        if close.shape[0] == 1:
            c_ref.at[i, 'MAG_PSF_F2'] = close['MAG_PSF']
            c_ref.at[i, 'MAGERR_PSF_F2'] = close['MAGERR_PSF']
            c_ref.at[i, 'ZeroPoint_F2'] = close['ZeroPoint']
            c_ref.at[i, 'ZeroPoint_rms_F2'] = close['ZeroPoint_rms']

    out = c_ref.loc[(c_ref['MAG_PSF_F2'] > 0.) & (c_ref['MAGERR_PSF_F2'] > 0.)]
    out.to_csv(outname)


def make_plot(file_in):
    f = pd.read_csv(file_in)
    bin_no = 20

    f['err'] = np.sqrt(f['MAGERR_PSF']**2 + f['MAGERR_PSF_F2']**2)
    f['dmag'] = A*np.abs(f['MAG_PSF']-f['MAG_PSF_F2'])
    f['total_err'] = np.sqrt(f['MAGERR_PSF']**2 + f['MAGERR_PSF_F2']**2 +
                             f['ZeroPoint_rms']**2 + f['ZeroPoint_rms_F2']**2)

    zp_rms = np.sqrt(f['ZeroPoint_rms']**2+f['ZeroPoint_rms_F2']**2)

    bin_err = ss.binned_statistic(f['MAG_PSF'], f['err'], statistic='mean', bins=bin_no)
    bin_dmag = ss.binned_statistic(f['MAG_PSF'], f['dmag'], statistic='mean', bins=bin_no)
    bin_terr = ss.binned_statistic(f['MAG_PSF'], f['total_err'], statistic='mean', bins=bin_no)

    mag_bin = bin_err[1][:-1] + (bin_err[1][2]-bin_err[1][1])*0.5  # mag value in the center of bin
    err_bin = bin_err[0]
    dmag_bin = bin_dmag[0]
    t_err = bin_terr[0]

    ab = np.polyfit(err_bin, dmag_bin, 1)
    froot = file_in.split('/')[-1].split('.')[0]

    # FIG 1 linear fit to binned errors and delta mags
    fig, ax = plt.subplots()
    ax.plot(np.log10(err_bin), np.log10(dmag_bin), '.')
    ax.plot(np.log10(err_bin), np.log10(np.polyval(ab, err_bin)))
    ax.set_title(froot)
    ax.text(0.05, 0.9, 'slope= %6.4f, intercept=%6.4f' % (ab[0], ab[1]), transform=ax.transAxes)
    ax.set_xlabel(r'$\log_{10}\sqrt{\sigma_1^2 +\sigma_2^2}$')
    ax.set_ylabel(r'$\log_{10}\langle\mid m_1-m_2 \mid\rangle$')

    out_name = froot + "_1.pdf"

    fig.savefig(out_name, facecolor="white", bbox_inches='tight')
    plt.clf()
    plt.close()

    # FIG 2
    fig, ax = plt.subplots()
    ax.errorbar(mag_bin, dmag_bin/t_err, yerr=t_err, fmt='o')
    ax.set_title(froot)
    ax.axhline(1, c='r')
    ax.set_xlabel(r'mag')
    ax.set_ylabel(r'$\langle\mid m_1 - m_2 \mid\rangle / err $')
    plt.text(0.05, 0.9, r'err$ = \sqrt{\sigma_1^2 +\sigma_2^2 + ZP_{rms,1}^2 + ZP_{rms,2}^2}$', transform=ax.transAxes)
    out_name = froot + "_2.pdf"
    fig.savefig(out_name, facecolor="white", bbox_inches='tight')
    plt.clf()
    plt.close()

    return ab[0], ab[1], np.mean(zp_rms)


def compare_catalogs():
    flist = find_files('./res', pattern='D0*.csv')

    ref = flist[0]
    c1 = read_cat(ref, add_fields=True)
    ref_name = ref.split('/')[-1].split('_')[0]

    for f2 in flist[1:]:
        out_name = ref_name + "_" + f2.split('/')[-1].split('_')[0]+".csv"
        get_comp(c1, f2, outname=out_name)


def summary_plot(tab, tname, labels):
    # FIG 3 slopes and ZP for all epochs
    fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True)
    ax1.set_title("REF: " + tname)
    ax1.plot(tab['slope'], '.', label='slopes', color='g')
    ax1.legend(loc=4)
    ax2.plot(tab['intr'], 'x', label='intercepts')
    ax2.plot(tab['zp'], '+', label='$\sqrt{ZP_{rms,1}^2 + ZP_{rms,2}^2}$')
    ax2.set_ylim(0.001, 0.033)
    ax2.legend(loc=1)
    plt.subplots_adjust(wspace=0, hspace=0)
    xs = range(tab.shape[0])
    plt.xticks(xs, labels, rotation='vertical')
    fig.savefig('summary.pdf', facecolor="white", bbox_inches='tight')
    plt.close()


def get_plots():
    flist = find_files('./', pattern='D*D*.csv')
    columns = ['file', 'slope', 'intr', 'zp']
    res = pd.DataFrame(columns=columns)
    labels = []
    for item in flist:
        a, b, zp = make_plot(item)
        new = pd.DataFrame([(item, a, b, zp)], columns=columns)
        res = res.append(new, ignore_index=True)
        labels.append(item.split('/')[-1][10:19])

    ref_star = flist[0].split('/')[-1][0:9]
    summary_plot(res, ref_star, labels)


def main():
    # compare_catalogs()
    get_plots()


if __name__ == "__main__":
    main()
