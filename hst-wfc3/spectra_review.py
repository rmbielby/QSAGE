def matching(ra1,dec1,ra2,dec2):
    import numpy as np
    match  = np.zeros(len(ra1))-1
    for i in np.arange(len(ra1)):
        mcoff  = ((ra1[i]-ra2)**2+(dec1[i]-dec2)**2)**0.5
        print 'mcoff = ',np.min(mcoff)*3600.
        if np.min(mcoff) < 1.6/3600.:
            match[i] = np.argmin(mcoff)
    return match.astype(np.int)


import numpy as np
import matplotlib.pyplot as plt
import astropy.io.fits as fits
import glob
from scipy.ndimage.filters import gaussian_filter1d
import seaborn as sns
from matplotlib.patches import Ellipse
import os,sys
import idlsave
from shutil import copyfile
from optparse import OptionParser

parser = OptionParser()

parser.add_option('-q','--quasar',dest='Quasar',
       help='quasar/field name',type='string',
       default='HB890232-042')
parser.add_option('-w','--wfc3dir',dest='wfc3dir',
       help='quasar/field name',type='string',
       default=None)
parser.add_option('-p','--phzfile',dest='phzfile',
       help='Photo-z file',type='string',
       default=None)
parser.add_option('-m','--musefile',dest='musefile',
       help='MUSE MARZ output file',type='string',
       default=None)

try:
    options,args = parser.parse_args(sys.argv[1:])
except:
    print "Error ... check usage with -h ";
    sys.exit(1)


# Set field Parameters
field   = options.Quasar
homedir = os.getenv("HOME")



# Set plit Parameters
xtitle   = r'Wavelength ($\AA$)'
ytitle   = r'Flux 10$^{-17}$ erg s$^{-1}$ cm$^{-2}$ $\AA^{-1}$'
fillcol  = sns.xkcd_rgb['pale orange']
errors   = 1
smfctr   = 0.64
linecols = ['dark green','denim blue','pale red','deep purple','burnt orange','slate','pink']

linelam  = np.array([1216.,1400.,1549.,1909.,2798.,3727.,4101.7,4340.5,4861.,5007.,6300.,6563.,6723.,7135.,8237.,9069.,9531.])
linenm   = np.array(['Lya','SiIV','CIV','CIII]','MgII','[OII]',r'H$\delta$',r'H$\gamma$',r'H$\beta$','[OIII]','[OI]','Ha','[SII]','[ArIII]','HeII','[SIII]','[SIII]'])
abslam   = np.array([4000.,5892.])
absnm    = np.array(['D4000','Na'])

keylines = np.array([3727.,4000.,4101.7,4340.5,4861.,5007.,6562.])
keyname  = ['[OII]','D4000',r'H$\delta$',r'H$\gamma$',r'H$\beta$','[OIII]',r'H$\alpha$']


# Read in F140W image & catalogue
fitsfile = '{1}/Dropbox/QSAGEdata/hst-wfc3/{0}/MasterCatalogs/{0}_F140W.fits'.format(field,homedir)
print fitsfile
f140w_hdulist = fits.open(fitsfile)
stack = f140w_hdulist[0].data
catfile = '{1}/Dropbox/QSAGEdata/hst-wfc3/{0}/MasterCatalogs/{0}_F140W_v2.cat'.format(field,homedir)
print catfile
gid,gx,gy,gra,gdec,gtheta,gaimage,gbimage,gmag = np.genfromtxt(catfile,
    unpack=True,usecols=(0,1,2,3,4,7,8,9,18))

# Check for MUSE data
if options.musefile:
    mzfile = options.musefile
else:
    mzfile = '{1}/Dropbox/QSAGE/data/vlt-muse/{0}/linecombine/COMBINED_IMAGE_F140W_marz_RMB.mz'.format(field,homedir)


if os.path.isfile(mzfile):
    data    = np.genfromtxt(mzfile,unpack=True,delimiter=',',dtype=None)
    mus_sel = np.where((data['f13'] >= -10.) & (data['f12'] >= -10.0) & (data['f12'] <= 16.))[0]
    mus_ra  = data['f2'][mus_sel].astype(np.float) #*180./np.pi  -0.000644
    mus_dec = data['f3'][mus_sel].astype(np.float) #*180./np.pi  -0.000397
    mus_z   = data['f12'][mus_sel].astype(np.float)
    mus_mag = data['f12'][mus_sel].astype(np.float)*0.
    mqop     = data['f13'][mus_sel]

# Check for photo-z's
if options.phzfile:
    phzfile = options.phzfile
else:
    phzfile = '{1}/Dropbox/QSAGE/data/hst-wfc3/{0}/MatchedData/{0}_CFHTLS_ugrizF140WF160W_redshifts_1802.cat'.format(field,homedir)

if os.path.isfile(phzfile):
    phzdata  = np.genfromtxt(phzfile,unpack=True,dtype=None)
    phzid    = phzdata['f0']
    pra      = phzdata['f1'].astype(np.float)
    pdec     = phzdata['f2'].astype(np.float)
    imag     = phzdata['f11'].astype(np.float)
    phz      = phzdata['f19'].astype(np.float)

# Read in fitted parameters & get spectra files
name1,zlist   = np.genfromtxt('redshifts.cat',unpack=True,skip_header=1)
name2,linewllist,SNlist = np.genfromtxt('lines.cat',unpack=True,skip_header=1)
speclist = np.array(glob.glob('ExtractedSpectra/*fits'))
nspec    = len(speclist)

Qualarr  = np.zeros(nspec)
zarray   = np.zeros(nspec)
zarrayalt= np.zeros(nspec)
namearr  = np.zeros(nspec)
startind = 0
endind   = len(speclist)-1

carryon = raw_input('Resume from zchoicepartway.dat?')
zchfile = 'zchoicepartway.dat'

if carryon == 'y':
    if os.path.isfile(zchfile):
        zchtable = np.genfromtxt(zchfile,unpack=True,dtype=None)
        namearr,Qualarr,zarray,zarrayalt,notearr = zchtable['f0'].astype(np.int),zchtable['f1'].astype(np.int),zchtable['f2'].astype(np.float),zchtable['f3'].astype(np.float),zchtable['f4']
        try:
            startind = np.where(Qualarr == 0)[0][0]
        except:
            startind = 0
        print startind
        copyfile(zchfile,'zchoicepartway_backup.dat')
elif carryon == 'n':
    if os.path.isfile(zchfile):
        copyfile(zchfile,'zchoicepartway_backup.dat')
    startind = 0
else:
    if os.path.isfile(zchfile):
        zchtable = np.genfromtxt(zchfile,unpack=True,dtype=None)
        namearr,Qualarr,zarray,zarrayalt,notearr = zchtable['f0'].astype(np.int),zchtable['f1'].astype(np.int),zchtable['f2'].astype(np.float),zchtable['f3'].astype(np.float),zchtable['f4']
        startind = np.where(namearr == np.float(carryon))[0]
        print startind
        copyfile(zchfile,'zchoicepartway_backup.dat')

for q in np.arange(nspec):
    split  = speclist[q].split('ID')
    split2 = split[1].split('_1D')
    namearr[q]  = np.float(split2[0])
spsort = np.argsort(namearr)
speclist = speclist[spsort]
namearr  = namearr[spsort]


#****start main loop
for q in np.arange(startind,len(namearr)):
    print q,speclist[q],namearr[q]
    name = namearr[q]
    # Read in spectrum from fits file
    hdulist = fits.open(speclist[q])
    wl       = hdulist[-1].data['wl']
    wl       = wl + (wl[1]-wl[0])
    flux     = hdulist[-1].data['flux_med']
    flux_err = hdulist[-1].data['ferr']
    fl_sm    = gaussian_filter1d(flux,smfctr)
    fler_sm  = gaussian_filter1d(flux_err,smfctr)
    # Identify with IDs from fitting files
    z        = zlist[np.where(name1 == name)[0]]
    range    = np.where((wl > 1.1e4) & (wl < 1.68e4))[0]
    x1       = np.min(wl)
    x2       = np.max(wl)
    y2       = 1.2 * np.nanmax(flux[range])
    y1       = np.nanmin([-0.1*np.abs(y2),np.nanpercentile(flux[range],2.4)])

# Check for matches in other catalogues
    csel = np.where(gid == name)[0]
    mzmatch = -1
    if os.path.isfile(mzfile):
        mzmatch    = matching(gra[csel],gdec[csel],mus_ra,mus_dec)[0]

    # Plot the spectrum
#    f,ax = plt.subplots(1,1,figsize=(12,4.8),num=1)
    f = plt.figure(num=1,figsize=(14.,5.6))
    plt.subplots_adjust(left=0.08,right=0.96,top=0.95,hspace=0.24)
    plt.pause(0.1)
    ax1 = plt.subplot2grid((2,5), (0, 0), colspan=3)
    ax2 = plt.subplot2grid((2,5), (1, 0), colspan=3)
    ax3 = plt.subplot2grid((2,5), (0, 3), rowspan=2,colspan=2)
    if errors == 1:
        ax1.plot(wl[range],fl_sm[range]+fler_sm[range],color=fillcol)
        ax1.plot(wl[range],fl_sm[range]-fler_sm[range],color=fillcol)
        ax1.fill_between(wl[range],fl_sm[range]+fler_sm[range],fl_sm[range]-fler_sm[range],color=fillcol)
    ax1.plot(wl[range],fl_sm[range],lw=3,color='k',zorder=12)
    ax1.plot(wl[range],fl_sm[range]*0.,lw=1.6,color='k',zorder=10,linestyle=':')
    if Qualarr[q] > 1:
        ax1.text(x1,y1,'{0},{1}'.format(zarray[q],Qualarr[q]))
        for l_index,line in enumerate(keylines):
            ax1.plot(np.array([line,line])*(1.+zarray[q]),np.array([0.5,0.6])*y2,color=sns.xkcd_rgb['blood red'])
            ax1.text(line*(1.+zarray[q]),0.5*y2,r'{0}$_W$'.format(keyname[l_index]),color=sns.xkcd_rgb['blood red'],size=10)
        # ax1.scatter((zarray[q]+1)*np.array([6562.,5007.]),np.array([y2,y2])/2.,marker='o',s=96,facecolor='none')
    print 'Number of z solns = ',len(z)
    if os.path.isfile(phzfile):
        pzmatch    = matching(gra[csel],gdec[csel],pra,pdec)[0]
        print pzmatch,phz[pzmatch]
        print 'z_p = {0}'.format(phz[pzmatch],pzmatch)
        for l_index,line in enumerate(keylines):
            ax1.plot(np.array([line,line])*(1.+phz[pzmatch]),np.array([0.2,0.25])*y2,color=sns.xkcd_rgb['denim blue'])
            ax1.text(line*(1.+phz[pzmatch]),0.26*y2,r'{0}$_p$'.format(keyname[l_index]),color=sns.xkcd_rgb['denim blue'],size=10)
            # ax1.text(x2,y1,r'$z_p=${0:5.3f}'.format(phz[pzmatch]),color='k',size=12,horizontalalignment='right')
            # ax1.scatter((1.+phz[pzmatch])*6562.,y2*0.6,marker='s',s=72,color=sns.xkcd_rgb['denim blue'])
        ax1.text(x2,y1+0.1*(y2-y1),r'$z_p={0}$'.format(phz[pzmatch]),horizontalalignment='right',color=sns.xkcd_rgb['denim blue'])
    for isoln in np.arange(len(z)):
        linez =linelam*(1+z[isoln])
        lsel  = np.where((linez>x1) & (linez<x2))[0]
        ytop = y2*(1.-(0.1*isoln))
        for j in lsel:
            ax1.plot([linez[j],linez[j]],[0., ytop],color=sns.xkcd_rgb[linecols[isoln]],linestyle='--')
            ax1.text(linez[j],ytop,linenm[j],color=sns.xkcd_rgb[linecols[isoln]])
        absz  = abslam*(1+z[isoln])
        asel  = np.where((absz>x1) & (absz<x2))[0]
        ax1.scatter(absz[asel],asel*0.+ytop,color=sns.xkcd_rgb[linecols[isoln]],marker='v',s=96)
        for j in asel:
            ax1.text(linez[j],ytop,absnm[j],color=sns.xkcd_rgb[linecols[isoln]])
        ax1.text(x1*(1.01),ytop,'{0}'.format(isoln),color=sns.xkcd_rgb[linecols[isoln]])
    ax1.set_xlim(x1,x2)
    ax1.set_xlabel(xtitle)
    ax1.set_ylim(y1,y2)
    ax1.set_ylabel(ytitle)
    # Get the galaxy from the primary catalogue.

    ax2.plot(wl[range],fl_sm[range]*0.,lw=1.6,color='k',zorder=1,linestyle=':')
    n_ext = (len(hdulist)-1)/4
    print 'Number of rotation spectra found: ',n_ext
    for i in np.arange(n_ext):
        flux = 0
        for j in np.arange(4):
            ext = i*4+j+1
            try:
                wlind     = hdulist[ext].data['wl']
                flind     = hdulist[ext].data['flux']
                flux     += np.interp(wl,wlind,flind)/4.
            except:
                print "Couldn't find spectrum in extension {0}".format(ext)
        fl_sm    = gaussian_filter1d(flux,smfctr)
        ax2.plot(wl,flux)
    ax2.set_xlim(x1,x2)
    ax2.set_xlabel(xtitle)
    ax2.set_ylim(1.6*y1,1.6*y2)
    ax2.set_ylabel(ytitle)
    plt.pause(0.1)


    oy = np.int(gx[csel])
    ox = np.int(gy[csel])
    thwdth = 36
    print ox,oy,np.shape(stack)
    thumb = stack[ox-thwdth:ox+thwdth,oy-thwdth:oy+thwdth]
    ax3.imshow(thumb,vmin=np.nanpercentile(thumb,16),vmax=np.nanpercentile(thumb,96),origin='lower',cmap = sns.cubehelix_palette(256,start=0,rot=-5, dark=0.24, light=1.0, as_cmap=True),extent=(-thwdth,+thwdth,-thwdth,+thwdth),interpolation='nearest')
    dx,dy = gx-gx[csel],gy-gy[csel]
    for j in np.where((np.abs(dx)<thwdth) & (np.abs(dy)<thwdth))[0]:
        print dx[j],dy[j]
        arte = Ellipse(xy=[dx[j],dy[j]], width=gaimage[j]*5.2, height=gbimage[j]*5.2, angle=gtheta[j])
        # print gaimage[j],gbimage[j],gtheta[j]
        ax3.add_artist(arte)
        arte.set_clip_box(ax3.bbox)
        arte.set_alpha(0.8)
        arte.set_facecolor('none')
        arte.set_edgecolor('k')
        arte.set_lw(1.2)
        ax3.text(dx[j]+2,dy[j]+2,'{0:3.0f}'.format(gid[j]),color='k',size=16)
    if mzmatch > -1:
        ax1.text(x2,y1,r'$z_m=${0:5.3f} ({1:1.0f})'.format(mus_z[mzmatch],mqop[mzmatch]),color='k',size=12,horizontalalignment='right')
        ax1.scatter((1.+mus_z[mzmatch])*6562.,y2*0.5,marker='*',s=96,color=sns.xkcd_rgb['blood red'])

    plt.pause(0.1)


    print name
    again='y'
    while again == 'y':
        again = 'n'
        Qual = ''
        # Read input from the terminal:
        while (Qual > 6) or (Qual < 0):
            Qual      = raw_input('is this a good (4), medium (3), hmmm (2), crap (1) emission line spec? continuum(5) z? a star (6)? Press 0 to skip without editing entry.')
            if Qual == 'a':
                notearr[i] = raw_input('Enter note for output file (no spaces): ')
                Qual = -99
            else:
                Qual = np.float(Qual)
        print 'you selected: ',Qual
        if Qual >= 1:
            Qualarr[q]= Qual
        if (Qual < 5) & (Qual >1):
            redshift = -1
           # Read input from the terminal:
            redshift = np.int(raw_input('which redshift solution is best? 0, 1, 2, 3... type -1 if all bad or 10 if single line'))
            print 'you selected: ',redshift
            if (redshift > -1) & (redshift < 10):
                zarray[q]=z[redshift]
                altredshift = -1
                altredshift = np.int(raw_input('is there a secondary solution? type -1 if rest are bad  '))
                print 'you selected:  ',altredshift
                if altredshift > -1:
                    zarrayalt[q]=z[altredshift]
            elif (np.float(redshift) == 10):
                single_line_wavelength = np.float(raw_input('Enter wavength for single emission line (in Angstroms): '))
                singlams = np.array([3727.,5007.,6562.])
                for isoln in np.arange(3):
                    redsh = single_line_wavelength/singlams[isoln]-1.
                    linez = linelam*(1+redsh)
                    lsel  = np.where((linez>x1) & (linez<x2))[0]
                    ytop  = y1+(y2-y1)*0.5
                    ybot  = y1+(y2-y1)*((0.1*isoln))
                    for j in lsel:
                        ax1.plot([linez[j],linez[j]],[ybot, ytop],color=sns.xkcd_rgb[linecols[isoln]],linestyle='--')
                        ax1.text(linez[j],ybot,linenm[j],color=sns.xkcd_rgb[linecols[isoln]])
                    absz  = abslam*(1+redsh)
                    asel  = np.where((absz>x1) & (absz<x2))[0]
                    ax1.scatter(absz[asel],asel*0.+ytop,color=sns.xkcd_rgb[linecols[isoln]],marker='v',s=96)
                    for j in asel:
                        ax1.text(absz[j],ybot,absnm[j],color=sns.xkcd_rgb[linecols[isoln]])
                    ax1.text(x1*(1.01),ytop,'{0}'.format(isoln),color=sns.xkcd_rgb[linecols[isoln]])
                    ax1.text(x1*(1.01),ybot,'{0}'.format(isoln),color=sns.xkcd_rgb[linecols[isoln]])
                plt.pause(1.)
                single_line_name       = raw_input('Enter line name (OII, OIII, or Ha) or R for redo: ')
                if single_line_name == 'R':
                    again = 'y'
                elif single_line_name == 'OII':
                    zarray[q],zarrayalt[q] = single_line_wavelength/3727.-1.,single_line_wavelength/5007.-1.
                elif single_line_name == 'OIII':
                    zarray[q],zarrayalt[q] = single_line_wavelength/5007.-1.,single_line_wavelength/6562.-1.
                else:
                    zarray[q],zarrayalt[q] = single_line_wavelength/6562.-1.,single_line_wavelength/5007.-1.

        # again = raw_input('Would you like to redo that one?(y/n) ')
    breakind = q
    zchout = open(zchfile,'w')
    for i in np.arange(len(Qualarr)):
        zchout.write('{0:5.0f} {1:3.0f} {2:10.5f} {3:10.5f} {4}\n'.format(namearr[i],Qualarr[i],zarray[i],zarrayalt[i],notearr[i]))
    zchout.close()
    # idlsave.save(breakind,Qualarr,zarray,zarrayalt,namearr,filename='zchoicepartway.idl')
    plt.clf()
