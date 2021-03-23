# cd /home/timaeus/Desktop/Projects/Plants/UMN/ACME_trait_distribution/scripts

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib import cm
from scipy.signal import detrend
from scipy.stats import gaussian_kde, normaltest
from statsmodels.distributions.empirical_distribution import ECDF

comb_runs = pd.read_pickle('../script_input/flx_comb.p')
# manually add in NaN values for bad years
comb_runs[1].loc['2007','GPP_NT_CUT_REF'] = np.nan
comb_runs[10].loc['1995','GPP_NT_CUT_REF'] = np.nan
# Latin Hypercube trait values
tr = pd.read_csv('../script_input/lhd_100.csv')
# defaults from 1000 sample random 
tr_d = pd.read_csv('../script_input/subPFT_1000.csv')
# append defaults
tr = tr.append(tr_d.iloc[1001,:])
tr.index = range(0,102)

t2m = [y.T2m.mean() for y in comb_runs]

# evaluate after omitting large leaf lifespan values
# 90th percentile, data max, and data max + 1 sd
lls90 = [12.6,8.1,0.5,4.5,3.4,4.6,0.8,0.5,3.0,0.8,0.5,0.2,1.5,0.5]
lls_max = [24.0,8.2,0.5,7.9,6.8,15.5,0.8,0.5,3.4,0.9,0.5,0.2,3.0,0.5]
lls_omax = [26.4,10.1,1.1,10.0,8.5,18.0,2.1,2.0,5.1,2.3,1.8,1.5,5.2,2.2]

incl = []
# set up omission dataframe 
lls_lim = lls_max
i_0 = tr.loc[:,'lls'].values < lls_lim[0]
i_13 = [tr.loc[:,'lls.'+str(i)].values < lls_lim[i] for i in range(1,14)]
for i in range(14):
    if i == 0:
	    incl.append(i_0)
    else:
	    incl.append(i_13[i-1])

incl = pd.DataFrame(incl)

incl = incl.transpose()

site = ['AU-Tum','BR-Sa1','FI-Hyy','FR-Pue','GF-Guy','RU-Cok','RU-Fyo',
	'RU-Sam','US-Blo','US-Ha1','US-PFa','US-SRM','US-UMB','US-Wkg',
	'ZA-Kru']

PFT = ['BET-te','BET-tr','NET-b','BDT(80)-NET(20)-te','BET-tr','BDS-b',
       'NET-b','BDS-b','NET-te','BDT(85)-NET(15)-te','BDT(50)-NET(50)-te',
       'BDS-te','BDT-te','C4','C4(70)-BDT(30)-tr']

handle = ['Tree-Te$_2$','Tree-Tr$_2$','Tree-Bo$_1$','Tree*-Te$_3$',
	  'Tree-Tr$_1$','Shrub-Bo$_1$','Tree-Bo$_2$','Shrub-Bo$_2$',
	  'Tree-Te$_3$','Tree*-Te$_2$','Tree*-Te$_1$','Shrub-Ar',
	  'Tree-Te$_1$','Grass-Ar','Grass*-Ar']

norm = mpl.colors.Normalize(vmin=np.min(t2m),vmax=np.max(t2m))
BuRd = cm.get_cmap('RdBu_r',15)

site_PFT_num = [5,4,2,7,4,11,2,12,1,7,1,10,7,14,14]

def tr_loc_assign(pft_num):
    out = tr.iloc[:,(pft_num-1)*3:pft_num*3].copy()
    out.columns = ['sla','lcn','lls']
    return out

def var_plot(x):
    x = x.upper()
    # locations in dataframe for extracting different outputs
    if x=='GPP':
        start = 7
        end = 415
    elif x=='NPP':
        start = 8
        end = 416
    elif x=='MR':
        start = 9
        end = 417
    else:
        start = 10
        end = 418
    
    norm = mpl.colors.Normalize(vmin=-10,vmax=25)
    BuRd = cm.get_cmap('coolwarm',15)

    for i,s in enumerate(site):
        vals = comb_runs[i].iloc[:,range(start,end,4)]
        # pick out zero values
        pl = [iv != 0 for iv in vals.mean()[0:100]]
	    # remove values outside of trait data range
        pl = [pl[iv] and incl.iloc[iv,site_PFT_num[i]-1] for iv in range(len(pl))]	
        
	    # Histograms
	    # mean
        ym = vals.mean()
	    
        with plt.rc_context({'axes.edgecolor':BuRd(norm(t2m[i])),'axes.linewidth':3}):
            plt.figure(figsize=(1.5,1.5)) 
            n, bins, patches = plt.hist(ym[0:100][pl],20)
            plt.setp(patches,'edgecolor','k','facecolor','grey')
            ax = ym[0:100][pl].mean()
            ay = n.max()/3.0
            if x == 'GPP':
                if PFT[i][4:6] == 'tr':
                    plt.xlim([250,5500])
                elif handle[i][-2:] == 'Ar':
                    plt.xlim([0,750])
                else:
                    plt.xlim([0,2550])
            elif x == 'NPP':
                if PFT[i][4:6] == 'tr':	
                    plt.xlim([500,1500])
                elif handle[i][-2:] == 'Ar':
                    plt.xlim([0,400])
                else:
                    plt.xlim([0,850])
            elif x == 'MR':
                if PFT[i][4:6] == 'tr':
                    plt.xlim([500,2200])
                elif handle[i][-2:] == 'Ar':
                    plt.xlim([0,400])
                else:
                    plt.xlim([0,750])
            plt.annotate("", xy=(ax, 0), xytext=(ax, ay),
            arrowprops=dict(arrowstyle='simple',color='lightgreen'))
            plt.tick_params(axis='both', which='major', labelsize=8) 
            plt.title(handle[i], fontdict={'fontsize':10})
            #plt.xlabel(x+' [gC m$^{-2}$ yr$^{-1}$]') 
            #plt.ylabel('Counts')
            plt.savefig(s+'_'+x.lower()+'_lhd_mhist_T_nol.pdf',bbox_inches='tight')
            plt.close()


var_plot('GPP')
var_plot('NPP')
var_plot('MR')

# colorbar plot
fig = plt.figure(figsize=(1.5, 0.75))
ax = fig.add_axes([0.05, 0.5, 0.9, 0.15])
cmap = mpl.cm.coolwarm
norm = mpl.colors.Normalize(vmin=-10,vmax=25)
cb = mpl.colorbar.ColorbarBase(ax, cmap=cmap,norm=norm,orientation='horizontal')
cb.ax.tick_params(labelsize=8)
cb.set_label('$^{\circ}$C',size=8)
plt.savefig('colorbar.pdf')
plt.close('all')

# gpp range density plot
exi = [1,8,11]
plt.figure(figsize=(3,3))
for ix,i in enumerate(exi):
    valsGPP = comb_runs[i].iloc[:,range(7,415,4)]
    pl = [iv != 0 for iv in valsGPP.mean()[0:100]]
    pl = [pl[iv] and incl.iloc[iv,site_PFT_num[i]-1] for iv in range(len(pl))]
    mGPP = valsGPP.mean()[0:100][pl]
    densGPP = gaussian_kde(mGPP)
    xGPP = np.linspace(mGPP.min(),mGPP.max(),500)
    c = ['b','g','brown']
    label = ['Trop','Te/Bo','Arid']
    plt.plot(xGPP,densGPP(xGPP)/densGPP(xGPP).max(),color=c[ix],label=label[ix])
    plt.legend(loc='upper right',fontsize=10)
    plt.tick_params(axis='both', which='major', labelsize=16)
    #plt.xlabel('GPP [gC m$^{-2}$ yr$^{-1}$]',fontsize=16)
    plt.ylabel('Norm. Prob. Density',fontsize=16)
    plt.savefig('gpp_dens_ex.pdf',bbox_inches='tight')

# Density Plots
for i,s in enumerate(site):
    valsGPP = comb_runs[i].iloc[:,range(7,415,4)]
    valsNPP = comb_runs[i].iloc[:,range(8,416,4)]
    valsMR = comb_runs[i].iloc[:,range(9,417,4)]
    # pick out zero values
    pl = [iv != 0 for iv in valsGPP.mean()[0:100]]
    # remove values outside of trait data range
    pl = [pl[iv] and incl.iloc[iv,site_PFT_num[i]-1] for iv in range(len(pl))]	

    # density plots
    mGPP = valsGPP.mean()[0:100]
    mNPP = valsNPP.mean()[0:100]
    mMR = valsMR.mean()[0:100]
    yGPP = (mGPP[pl]-mGPP[pl].mean())/mGPP[pl].std()
    densGPP = gaussian_kde(yGPP)
    xGPP = np.linspace(yGPP.min(),yGPP.max(),100)
    yNPP = (mNPP[pl]-mNPP[pl].mean())/mNPP[pl].std()
    densNPP = gaussian_kde(yNPP)
    xNPP = np.linspace(yNPP.min(),yNPP.max(),100)
    yMR = (mMR[pl]-mMR[pl].mean())/mMR[pl].std()
    densMR = gaussian_kde(yMR)
    xMR = np.linspace(yMR.min(),yMR.max(),100)

    norm = mpl.colors.Normalize(vmin=-10,vmax=25)
    BuRd = cm.get_cmap('coolwarm',15)
    with plt.rc_context({'axes.edgecolor':BuRd(norm(t2m[i])),'axes.linewidth':2}):
	k2, p = normaltest(mGPP[pl])
	if p>0.05:
	    plt.plot(xGPP,densGPP(xGPP),'k')
	else:
	    plt.plot(xGPP,densGPP(xGPP),'k--')
	k2, p = normaltest(mNPP[pl])
	if p>0.05:
	    plt.plot(xNPP,densNPP(xNPP),'g')
	else:
	    plt.plot(xNPP,densNPP(xNPP),'g--')
	k2, p = normaltest(mMR[pl])
	if p>0.05:
	    plt.plot(xMR,densMR(xMR),'r')
	else:
	    plt.plot(xMR,densMR(xMR),'r--')
	plt.title(handle[i])
        plt.xlabel('Standardized Value')
        plt.ylabel('Probability Density')
        plt.savefig(s+'_lhd_dens_T.pdf',bbox_inches='tight')
        plt.close()

def site_info(x):
    x = x.upper()
    # locations in dataframe for extracting different outputs
    if x=='GPP':
	start = 7
	end = 415
    elif x=='NPP':
	start = 8
	end = 416
    elif x=='MR':
	start = 9
	end = 417
    else:
	start = 10
	end = 418
    
    norm = mpl.colors.Normalize(vmin=-10,vmax=25)
    BuRd = cm.get_cmap('coolwarm',15)

    for i,s in enumerate(site):
	vals = comb_runs[i].iloc[:,range(start,end,4)]
	# pick out zero values
	pl = [iv != 0 for iv in vals.mean()[0:100]]
	# remove values outside of trait data range
	pl = [pl[iv] and incl.iloc[iv,site_PFT_num[i]-1] for iv in range(len(pl))]	
    
	# Histograms
	# mean
    ym = vals.mean()
	ymc = ym[0:100][pl]
	print handle[i],ymc.min(),ymc.max()

site_info('gpp')

site_PFT2_num = [0,0,0,1,0,0,0,0,0,1,7,0,0,0,6]
site_PFT1_mean = list()
site_PFT2_mean = list()
for i,s in enumerate(site):
    valsGPP = comb_runs[i].iloc[:,range(7,415,4)]
    # pick out zero values
    pl = [iv != 0 for iv in valsGPP.mean()[0:100]]
    # remove values outside of trait data range
    pl = [pl[iv] and incl.iloc[iv,site_PFT_num[i]-1] for iv in range(len(pl))]
    
    site_trait1 = tr_loc_assign(site_PFT_num[i])
    site_PFT1_mean.append(site_trait1.iloc[0:100,:][pl].mean())
    if site_PFT2_num[i] == 0:
	site_PFT2_mean.append(np.nan)
    else:
	site_trait2 = tr_loc_assign(site_PFT2_num[i])
	site_PFT2_mean.append(site_trait2.iloc[0:100,:][pl].mean())
