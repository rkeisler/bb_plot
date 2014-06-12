import numpy as np
import matplotlib.pylab as pl
from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica'], 'size':13})
pl.ion()
datadir='./data/'


def make_many():
    make_one(force_crop=True, include_actpol=False, include_polarbear=True)
    make_one(force_crop=True, include_actpol=False, include_polarbear=False)
    make_one(force_crop=False, include_actpol=False, include_polarbear=True)
    make_one(force_crop=False, include_actpol=True, include_polarbear=True)



def make_one(lpower=2, include_polarbear=True, 
             include_actpol=False,
             force_crop=True, 
             logy=False, filetype='pdf', 
             savepath='./', savename=None):
    '''
    The main program for making this figure.
    LPOWER [1.5] controls the power of the multipole scaling.
    INCLUDE_POLARBEAR [True] determines whether to include the POLARBEAR data.
    INCLUDE_ACTPOL [True] determines whether to include the ACTPOL data.
    FORCE_CROP [True] determines whether to use the "BICEP2+SPTpol" vertical range.
    FILETYPE ['pdf'] is the type of image produced.
    SAVEPATH ['./'] is the save path.
    SAVENAME [None] is available if you want to force a particular filename.
    '''

    # load theory
    l, cl_bb_lens, cl_bb_r = load_theory(r=0.2)

    # load data
    spt = load_data('sptpol', lpower)
    bi = load_data('bicep2', lpower)
    pb = load_data('polarbear', lpower)
    actpol = load_data('actpol',lpower)

    # Set up a few things for plotting.
    lscal = lscaling(l, lpower)
    fs = 18
    lw_theory = 2
    lw_data = 2
    fplot = pl.semilogx
    if logy: fplot = pl.loglog
    pl.clf()
    #pl.xlim(3,6e3)
    pl.xlim(20,5e3)
    pl.xlabel(r'$\ell$', fontsize=fs)
    pl.ylabel(get_ylabel(lpower), fontsize=fs)

    # Plot the theory spectra.
    fplot(l, cl_bb_lens*lscal, 'k-.', linewidth=lw_theory)
    fplot(l, cl_bb_r*lscal, 'k--', linewidth=lw_theory)
    fplot(l, (cl_bb_lens+cl_bb_r)*lscal, 'k', linewidth=lw_theory)

    # Plot the data.
    exp2plot = [bi, spt]
    if include_polarbear: exp2plot.append(pb)
    if include_actpol: exp2plot.append(actpol)
    ax = pl.gca()
    ms_scale = 1.
    if not(force_crop): ms_scale *= 0.5
    if (lpower==2): ms_scale *= 0.5
    for e in exp2plot:
        #if (force_crop) & (e['name']=='polarbear') & (not(logy)):
        #    ylim_current = np.array(ax.get_ylim())
        #    pl.ylim(1.*ylim_current)
        pl.errorbar(e['l'], e['plot'], yerr=e['dplot'], fmt=' ', 
                    linewidth=lw_data, color=e['color'], capsize=e['capsize'])
        pl.plot(e['l'], e['plot'], e['symbol'], color=e['color'], ms=9*ms_scale)
        if e['name']=='sptpol': 
            pl.plot(e['l'], e['plot']+e['dplot'], e['symbol'], color=e['color'], ms=7*ms_scale)
            pl.plot(e['l'], e['plot']-e['dplot'], e['symbol'], color=e['color'], ms=7*ms_scale)

    if not(logy):
        if force_crop: pl.ylim(-0.02,0.14)

    if logy:
        theory_max = np.max((cl_bb_lens+cl_bb_r)*lscal)
        if force_crop: pl.ylim(theory_max/150., theory_max*6.)
        else: pl.ylim(theory_max/150., theory_max*15.)

    # Add a legend.
    yl = ax.get_ylim()
    xl = ax.get_xlim()
    fs_legend = 19
    if logy:
        xt = 2.*xl[0]
        yt = 0.3*yl[1]
        dyt = 1.8
    else:
        xt = 2.*xl[0]
        yt = 0.85*(yl[1]-yl[0]) + yl[0]
        dyt = (yl[1]-yl[0])/14.
    this_yt = yt
    for i,e in enumerate(exp2plot):
        if logy: this_yt /= dyt
        else: this_yt = yt-i*dyt
        pl.text(xt,this_yt,e['legend_name'], color=e['color'], fontsize=fs_legend)


    # Save the figure.
    if savename==None:
        savename = savepath+'bb_l%0.1f_bicep2_sptpol'%lpower
        if include_polarbear: savename += '_pbear'
        if include_actpol: savename += '_actpol'
        if not(force_crop): savename += '_nocrop'
        if logy: savename += '_logy'
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
        legend_name = 'SPTpol (lensing only)'
        name = 'sptpol'
        color = 'Navy'
        symbol = '^'
        capsize=0
        
    if exp=='bicep2':
        # Load BICEP2 bandpowers.
        tmp=np.loadtxt(datadir+'B2_3yr_bandpowers_20140314.txt')
        l = tmp[:,1]
        cl_bb = tmp[:,6]/l/(l+1.)*2.*np.pi
        sigma_cl_bb = tmp[:,12]/l/(l+1.)*2.*np.pi
        legend_name = 'BICEP2'
        name = 'bicep2'
        color = 'DarkRed'
        symbol = 'o'
        capsize=3

    if exp=='polarbear':
        # Load PolarBear bandpowers.
        l = np.array([700., 1100., 1500., 1900.])
        dl_bb = np.array([0.093, 0.149, -0.317, 0.487])
        sigma_dl_bb = np.array([0.056, 0.117, 0.236, 0.482])
        cl_bb = dl_bb/l/(l+1.)*2.*np.pi
        sigma_cl_bb = sigma_dl_bb/l/(l+1.)*2.*np.pi
        legend_name = 'POLARBEAR'
        name = 'polarbear'
        color = 'DarkOrange'
        symbol = 'o'
        capsize=3


    if exp=='actpol':
        actpol = load_actpol()
        l = actpol['l']
        dl_bb = actpol['dl_bb']
        sigma_dl_bb = actpol['sigma_dl_bb']
        cl_bb = dl_bb/l/(l+1.)*2.*np.pi
        sigma_cl_bb = sigma_dl_bb/l/(l+1.)*2.*np.pi
        legend_name = 'ACTpol'
        name = 'actpol'
        color = 'darkgreen'
        symbol = 'o'
        capsize=3

    return {'l':l, 'cl':cl_bb, 'dcl':sigma_cl_bb, 
            'name':name, 'legend_name':legend_name, 
            'color':color, 
            'plot':cl_bb*lscaling(l, lpower),
            'dplot':sigma_cl_bb*lscaling(l, lpower), 
            'symbol':symbol, 'capsize':capsize}

def lscaling(l, lpower=1.):
    if (lpower==0) | (lpower==0.5) | (lpower==1): output = l**(lpower)
    if (lpower==1.5): output = l**(0.5)*(l+1.)
    if (lpower==2): output = l*(l+1.)/2./np.pi
    return output


def load_actpol():
    f = open('data/actpol_1405_5524_data.txt','r')
    d = {'l':[], 'dl_tt':[], 'sigma_dl_tt':[], 
         'dl_te':[], 'sigma_dl_te':[], 
         'dl_ee':[], 'sigma_dl_ee':[], 
         'dl_bb':[], 'sigma_dl_bb':[]}
    for line in f:
        tmp = line.split('&')
        d['l'].append( np.float(tmp[0]) )
        d['dl_tt'].append( np.float(tmp[2]) )
        d['sigma_dl_tt'].append( np.float(tmp[3]) )
        d['dl_te'].append( np.float(tmp[4]) )
        d['sigma_dl_te'].append( np.float(tmp[5]) )
        d['dl_ee'].append( np.float(tmp[6]) )
        d['sigma_dl_ee'].append( np.float(tmp[7]) )
        d['dl_bb'].append( np.float(tmp[8]) )
        d['sigma_dl_bb'].append( np.float(tmp[9]) )

    for k in d.keys():
        d[k] = np.array(d[k])
    return d

if (__name__=='__main__'): make_baseline()
