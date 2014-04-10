import numpy as np
import matplotlib.pylab as pl
from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica'], 'size':13})
pl.ion()
datadir='./data/'


def make_baseline():
    main(savename='bb.png')


def main(lpower=1.5, include_polarbear=True, force_crop=True, 
         filetype='png', savepath='./', savename=None):

    # load theory
    l, cl_bb_lens, cl_bb_r = load_theory(r=0.2)

    # load data
    spt = load_data('sptpol', lpower)
    bi = load_data('bicep2', lpower)
    pb = load_data('polarbear', lpower)
    
    # Set up a few things for plotting.
    lscal = lscaling(l, lpower)
    fs = 18
    lw_theory = 2
    lw_data = 2
    fplot = pl.semilogx
    pl.clf()
    pl.xlim(3,6e3)
    pl.xlabel(r'$\ell$', fontsize=fs)
    pl.ylabel(get_ylabel(lpower), fontsize=fs)

    # Plot the theory spectra.
    fplot(l, cl_bb_lens*lscal, 'k-.', linewidth=lw_theory)
    fplot(l, cl_bb_r*lscal, 'k--', linewidth=lw_theory)
    fplot(l, (cl_bb_lens+cl_bb_r)*lscal, 'k', linewidth=lw_theory)

    # Plot the data.
    exp2plot = [bi, spt]
    if include_polarbear: exp2plot.append(pb)
    ax = pl.gca()
    for e in exp2plot:
        if (force_crop) & (e['name']=='polarbear'):
            pl.ylim(ax.get_ylim())
        pl.errorbar(e['l'], e['plot'], yerr=e['dplot'], fmt=' ', 
                    linewidth=lw_data, color=e['color'], capsize=3)
        pl.plot(e['l'], e['plot'], e['symbol'], color=e['color'], ms=9)
        if e['name']=='sptpol': 
            pl.plot(e['l'], e['plot']+e['dplot'], e['symbol'], color=e['color'], ms=7)
            pl.plot(e['l'], e['plot']-e['dplot'], e['symbol'], color=e['color'], ms=7)

    # Add a legend.
    yl = ax.get_ylim()
    xl = ax.get_xlim()
    fs_legend = 19
    xt = 2.*xl[0]
    yt = 0.85*(yl[1]-yl[0]) + yl[0]
    dyt = yt/8.
    for i,e in enumerate(exp2plot):
        pl.text(xt,yt-i*dyt,e['legend_name'], color=e['color'], fontsize=fs_legend)


    # Save the figure.
    if savename==None:
        savename = savepath+'bb_l%0.1f_bicep2_sptpol'%lpower
        if include_polarbear: savename += '_pbear'
        savename += '.'+filetype
    print 'making '+savename
    pl.savefig(savename, dpi=300)



def get_ylabel(lpower):
    if lpower==0: o=r'$C_\ell^{BB}$'
    if lpower==0.5: o=r'$\ell^{0.5} C_\ell^{BB}$'
    if lpower==1: o=r'$\ell C_\ell^{BB}$'
    if lpower==1.5: o=r'$\ell^{0.5} ( \ell + 1)  C_\ell^{BB}$'
    if lpower==2: o=r'$\ell ( \ell + 1) C_\ell^{BB} /(2 \pi)$'
    o += '  '
    o += '$[\mu K^2]$'
    return o


def load_theory(r=0.2):
    # Read in the theory spectra.
    cambfile = datadir+'planck2013_TableIICol4_lensedCls.dat'
    tmp = np.loadtxt(cambfile)
    l = tmp[:,0]
    nl = len(l)
    dl_bb_lens = tmp[:,3]
    cl_bb_lens = dl_bb_lens/l/(l+1.)*2.*np.pi

    tmp = np.loadtxt(datadir+'delta_clbb_delta_r_uk2.txt')
    r=0.2
    cl_bb_r = r*tmp[0:nl, 1]
    return l, cl_bb_lens, cl_bb_r


def load_data(exp, lpower):
    if exp=='sptpol':
        # Load SPTpol Hanson et al bandpowers.
        tmp=np.loadtxt(datadir+'sptpol_hanson13_CL_BB.txt',delimiter=',')
        l = tmp[:,0]
        cl_bb = tmp[:,1]/l/1e4
        sigma_cl_bb = tmp[:,2]/l/1e4
        legend_name = 'SPTpol (lensing)'
        name = 'sptpol'
        color = 'blue'
        symbol = 'v'
        
    if exp=='bicep2':
        # Load BICEP2 bandpowers.
        tmp=np.loadtxt(datadir+'B2_3yr_bandpowers_20140314.txt')
        l = tmp[:,1]
        cl_bb = tmp[:,6]/l/(l+1.)*2.*np.pi
        sigma_cl_bb = tmp[:,12]/l/(l+1.)*2.*np.pi
        legend_name = 'BICEP2'
        name = 'bicep2'
        color = 'darkred'
        symbol = 'o'


    if exp=='polarbear':
        # Load PolarBear bandpowers.
        l = np.array([700., 1100., 1500., 1900.])
        dl_bb = np.array([0.093, 0.149, -0.317, 0.487])
        sigma_dl_bb = np.array([0.056, 0.117, 0.236, 0.482])
        cl_bb = dl_bb/l/(l+1.)*2.*np.pi
        sigma_cl_bb = sigma_dl_bb/l/(l+1.)*2.*np.pi
        legend_name = 'POLARBEAR'
        name = 'polarbear'
        color = 'green'
        symbol = 'o'

    return {'l':l, 'cl':cl_bb, 'dcl':sigma_cl_bb, 
            'name':name, 'legend_name':legend_name, 
            'color':color, 
            'plot':cl_bb*lscaling(l, lpower),
            'dplot':sigma_cl_bb*lscaling(l, lpower), 
            'symbol':symbol}

def lscaling(l, lpower=1.):
    if (lpower==0) | (lpower==0.5) | (lpower==1): output = l**(lpower)
    if (lpower==1.5): output = l**(0.5)*(l+1.)
    if (lpower==2): output = l*(l+1.)/2./np.pi
    return output


