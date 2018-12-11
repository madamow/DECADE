import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import scipy.stats as ss
import os
import fnmatch


def find_files(cat, ptr=''):
    files = []

    for root, dirnames, filenames in os.walk(cat):
        for filename in fnmatch.filter(filenames, ptr):
            files.append(os.path.join(root, filename))
    return files


def read_cat(cat_file):
    # read cat and remove bad values for magnitudes
    c = pd.read_csv(cat_file)
    new = c.loc[(c['MAG_PSF'] > -10.) & (np.abs(c['SPREAD_MODEL']) < 0.003)]
    new['diff'] = 100.
    new['err2'] = -9.
    return new


def get_comp(c_ref, f2, outname=''):

    c2 = read_cat(f2)

    for i, r1 in c_ref.iterrows():
        close = c2.loc[(np.abs(r1['ALPHAWIN_J2000'] - c2['ALPHAWIN_J2000']) < 0.00001) &
                       (np.abs(r1['DELTAWIN_J2000'] - c2['DELTAWIN_J2000']) < 0.00001)]
        if close.shape[0] == 1:
            diff = r1['MAG_PSF']-close['MAG_PSF']
            err2 = np.sqrt(close['ZeroPoint_rms']**2 + close['MAGERR_PSF']**2 +
                           r1['ZeroPoint_rms']**2 + r1['MAGERR_PSF']**2)
            c_ref.at[i, 'diff'] = diff
            c_ref.at[i, 'err2'] = err2

    out = c_ref.loc[(c_ref['diff'] < 10.) & (np.abs(c_ref['err2']) > 0.)]
    out.to_csv(outname)


def make_plot(file_in):
    print file_in
    f = pd.read_csv(file_in)
    out_name = file_in.split('/')[-1].split('.')[0] + ".pdf"

    fig = plt.figure()
    ax = fig.add_subplot(111)

    bin_err = ss.binned_statistic(f['MAG_PSF'], f['err2'], statistic='median', bins=100)
    ax.plot(f['MAG_PSF'], np.abs(f['diff']), '.')
    ax.plot(bin_err[1][:-1], bin_err[0], '.')
    ax.set_xlabel('MAG_PDF')

    fig.savefig(out_name, facecolor="white", bbox_inches='tight')


def compare_catalogs():
    flist = find_files('./res', ptr='D0*')

    ref = flist[0]
    c1 = read_cat(ref)
    ref_name = ref.split('/')[-1].split('_')[0]

    for f2 in flist[1:]:
        print f2
        out_name = ref_name + "_" + f2.split('/')[-1].split('_')[0]+".csv"
        get_comp(c1, f2, outname=out_name)


def get_plots():
    flist = find_files('./', ptr='D*D*.csv')
    for item in flist:
        make_plot(item)


def main():
   # compare_catalogs()
    get_plots()


if __name__ == "__main__":
    main()
