import numpy as np
import matplotlib.pylab as pl
from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica'], 'size':13})
pl.ion()
datadir='./data/'


def make_baseline():
    '''
    Run the make_one() function with the baseline parameters, and 
    save the output image as ee.png.
    '''
    make_one(savename='ee.png', lpower=2, include_wmap=True)



def make_many():
    make_baseline()
    make_one(lpower=2, include_wmap=True)
    make_one(lpower=2, include_wmap=False)



def make_one(lpower=1.5, include_wmap=True, force_crop=True, 
             logy=False, filetype='png', 
             savepath='./', savename=None):
    '''
    The main program for making this figure.
    LPOWER [1.5] controls the power of the multipole scaling.
    INCLUDE_WMAP [True] determines whether to include the WMAP data.
    FORCE_CROP [True] determines whether to use the "BICEP2+QUaD" vertical range.
    FILETYPE ['png'] is the type of image produced.
    SAVEPATH ['./'] is the save path.
    SAVENAME [None] is available if you want to force a particular filename.
    '''

    # load theory
    l, cl_ee_lens = load_theory()

    # load data
    wmap = load_data('wmap', lpower)
    bi = load_data('bicep2', lpower)
    quad = load_data('quad', lpower)
    
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
    fplot(l, cl_ee_lens*lscal, 'k-', linewidth=lw_theory)


    # Plot the data.
    exp2plot = [bi, quad]
    if include_wmap: exp2plot.append(wmap)
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
        #if e['name']=='sptpol': 
        #    pl.plot(e['l'], e['plot']+e['dplot'], e['symbol'], color=e['color'], ms=7*ms_scale)
        #    pl.plot(e['l'], e['plot']-e['dplot'], e['symbol'], color=e['color'], ms=7*ms_scale)


    if not(logy):
        theory_max = np.max((cl_ee_lens)*lscal)
        if force_crop: pl.ylim(-10,50)
        else: pl.ylim(-10,50)

    if logy:
        theory_max = np.max((cl_ee_lens)*lscal)
        if force_crop: pl.ylim(theory_max/150., theory_max*1.)
        else: pl.ylim(theory_max/150., theory_max*1.)

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
        savename = savepath+'ee_l%0.1f_bicep2_quad'%lpower
        if include_wmap: savename += '_wmap'
        if not(force_crop): savename += '_nocrop'
        if logy: savename += '_logy'
        savename += '.'+filetype
    print 'making '+savename
    pl.savefig(savename, dpi=300)



def get_ylabel(lpower):
    if lpower==0: o=r'$C_\ell^{EE}$'
    if lpower==0.5: o=r'$\ell^{0.5} C_\ell^{EE}$'
    if lpower==1: o=r'$\ell C_\ell^{EE}$'
    if lpower==1.5: o=r'$\ell^{0.5} ( \ell + 1)  C_\ell^{EE}$'
    if lpower==2: o=r'$\ell ( \ell + 1) C_\ell^{EE} /(2 \pi)$'
    o += '  '
    o += '$[\mu K^2]$'
    return o


def load_theory():
    # Read in the theory spectra.
    cambfile = datadir+'planck2013_TableIICol4_lensedCls.dat'
    tmp = np.loadtxt(cambfile)
    l = tmp[:,0]
    nl = len(l)
    dl_ee_lens = tmp[:,2]
    cl_ee_lens = dl_ee_lens/l/(l+1.)*2.*np.pi
    return l, cl_ee_lens


def load_data(exp, lpower):
    if exp=='wmap':
        # Load WMAP EE bandpowers.
        tmp=np.loadtxt(datadir+'wmap_binned_ee_spectrum_9yr_v5.txt')
        wh_keep = np.where(tmp[:,0] <= 200)[0] # only show L<LMAX WMAP data.
        tmp = tmp[wh_keep, :]
        l = tmp[:,0]
        dl_ee = tmp[:,3]
        sigma_dl_ee = tmp[:,4]
        cl_ee = dl_ee/l/(l+1.)*2.*np.pi
        sigma_cl_ee = sigma_dl_ee/l/(l+1.)*2.*np.pi
        legend_name = 'WMAP'
        name = 'wmap'
        color = 'blue'
        symbol = 'o'
        capsize=0
        
    if exp=='bicep2':
        # Load BICEP2 bandpowers.
        tmp=np.loadtxt(datadir+'B2_3yr_bandpowers_20140314.txt')
        l = tmp[:,1]
        cl_ee = tmp[:,5]/l/(l+1.)*2.*np.pi
        sigma_cl_ee = tmp[:,11]/l/(l+1.)*2.*np.pi
        legend_name = 'BICEP2'
        name = 'bicep2'
        color = 'darkred'
        symbol = 'o'
        capsize=3

    if exp=='quad':
        # Load QUAD bandpowers.
        tmp=np.loadtxt(datadir+'quad_ee.txt')
        l = tmp[:,0]
        dl_ee = tmp[:,1]
        sigma_dl_ee = tmp[:,2]
        cl_ee = dl_ee/l/(l+1.)*2.*np.pi
        sigma_cl_ee = sigma_dl_ee/l/(l+1.)*2.*np.pi
        legend_name = 'QUaD'
        name = 'quad'
        color = 'green'
        symbol = 'o'
        capsize=3

    return {'l':l, 'cl':cl_ee, 'dcl':sigma_cl_ee, 
            'name':name, 'legend_name':legend_name, 
            'color':color, 
            'plot':cl_ee*lscaling(l, lpower),
            'dplot':sigma_cl_ee*lscaling(l, lpower), 
            'symbol':symbol, 'capsize':capsize}

def lscaling(l, lpower=1.):
    if (lpower==0) | (lpower==0.5) | (lpower==1): output = l**(lpower)
    if (lpower==1.5): output = l**(0.5)*(l+1.)
    if (lpower==2): output = l*(l+1.)/2./np.pi
    return output


if (__name__=='__main__'): make_baseline()