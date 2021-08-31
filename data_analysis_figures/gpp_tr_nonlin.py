# cd /Users/timaeus/Desktop/Projects/Plants/UMN/ACME_trait_distribution/scripts/supp_figs

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

fol = '/home/timaeus/Desktop/Projects/Plants/UMN/ACME_trait_distribution'

f_comb = pd.read_pickle(fol+'/script_input/flx_comb.p')

# manually add in NaN values for bad years
f_comb[1].loc['2007','GPP_NT_CUT_REF'] = np.nan
f_comb[10].loc['1995','GPP_NT_CUT_REF'] = np.nan
# append constrained mean output
cm_out = pd.read_pickle(fol+'/script_input/yrly_cm.p')
comb_runs = list()
for i in range(15):
    comb_runs.append(f_comb[i].join(cm_out[i],rsuffix='cm'))
# Latin Hypercube trait values
tr = pd.read_csv(fol+'/script_input/lhd_100.csv')
# defaults from 1000 sample random 
tr_d = pd.read_csv(fol+'/script_input/subPFT_1000.csv')
# append defaults
tr = tr.append(tr_d.iloc[1001,:])
tr.index = range(0,102)

t2m = [y.T2m.mean() for y in comb_runs]
pre = [y.PRE.mean() for y in comb_runs]
eva = [y.Evap.mean() for y in comb_runs]
tra = [y.Tran.mean() for y in comb_runs]

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

# constrained mean trait values
sla = [0.020531,0.027451, 0.013063, 0.035984, 0.027451, 0.038787, 0.013063, 0.044086, 0.011963, 0.035984, 0.011891, 0.040858, 0.035984, 0.035190, 0.036441]
lcn = [30.047955, 27.810269, 39.846400, 22.405209, 27.810269, 20.259220, 39.846400, 19.033626, 47.813814, 22.405209, 46.732404, 19.575897, 22.405209, 31.865744, 31.757485]
lls = [2.233203, 2.253321, 3.997375, 0.381248, 2.253321, 0.253523, 3.997375, 0.107192, 5.770472, 0.381248, 4.913209, 0.342350, 0.381248, 0.239287, 0.291413]

site_PFT_num = [5,4,2,7,4,11,2,12,1,7,1,10,7,14,14]

def tr_loc_assign(pft_num):
    out = tr.iloc[:,(pft_num-1)*3:pft_num*3].copy()
    out.columns = ['sla','lcn','lls']
    return out

# nonlinear function to fit to trait-gpp relationships
def nfunc(x,b0,b1,p0):
    return b0+b1*x**p0

nparams = list()
tr_r2 = list()
gpp_pdiff = list()
sla_tr_z = list()
lcn_tr_z = list()
gpp_z = list()

for i,s in enumerate(site):
    # look at GPP for 0 values
    vals = comb_runs[i].iloc[:,range(7,419,4)]
    # pick out zero values 
    pl = [iv != 0 for iv in vals.mean()[0:100]]
    # make the line below into a for loop...
    pl = [pl[iv] and incl.iloc[iv,site_PFT_num[i]-1] for iv in range(len(pl))]
    # if plotting something other than GPP (NPP = 8, MR = 9)
    vals = comb_runs[i].iloc[:,range(7,419,4)]
    gpp = vals.mean()[0:100]
    gppcm = vals.GPPcm.mean()
    sit_tr = tr_loc_assign(site_PFT_num[i]).iloc[0:100,:]
    sit_trd = tr_loc_assign(site_PFT_num[i]).iloc[101,:]
    gppd = vals.GPP101.mean()
    gpp_z.append((gpp[pl]-np.mean(gpp[pl]))/np.std(gpp[pl]))
    sla_tr_z.append((sit_tr.sla[pl]-sit_tr.sla[pl].mean())/sit_tr.sla[pl].std())
    lcn_tr_z.append((sit_tr.lcn[pl]-sit_tr.lcn[pl].mean())/sit_tr.lcn[pl].std())
    gpp_pdiff.append((np.mean(gpp[pl])-gppcm)/gppcm)
    # save values for Nonlinear example
    if s=='RU-Sam':
        slaRu = sit_tr.sla[pl]
        lcnRu = sit_tr.lcn[pl]
        gppRu = gpp[pl]
        gppcmRu = vals.GPPcm.mean()
    elif s=='BR-Sa1':
        slaBr = sit_tr.sla[pl]
        lcnBr = sit_tr.lcn[pl]
        gppBr = gpp[pl]
        gppcmBr = vals.GPPcm.mean()
    elif s=='FR-Pue':
        slaFr = sit_tr.sla[pl]
        lcnFr = sit_tr.lcn[pl]
        gppFr = gpp[pl]
        gppcmFr = vals.GPPcm.mean() 
    elif s=='US-Blo':
        slaUs = sit_tr.sla[pl]
        lcnUs = sit_tr.lcn[pl]
        gppUs = gpp[pl]
        gppcmUs = vals.GPPcm.mean() 

evergreen = [0,1,2,4,6,8,10]
decid = [3,5,7,9,11,12,13,14]

ever_g = np.concatenate([gpp_z[e].values for e in evergreen])
ever_sla = np.concatenate([sla_tr_z[e].values for e in evergreen])
ever_lcn = np.concatenate([lcn_tr_z[e].values for e in evergreen])
deci_g = np.concatenate([gpp_z[e].values for e in decid])
deci_sla = np.concatenate([sla_tr_z[e].values for e in decid])
deci_lcn = np.concatenate([lcn_tr_z[e].values for e in decid])

# common colorbar from 0 to 45
plt.figure(figsize=[7.5,7.5])
plt.rcParams.update({'font.size':14})
plt.subplot(221)                                
plt.tight_layout()
h = plt.hist2d(deci_lcn,deci_g,cmin=0,vmin=0,vmax=45)
plt.text(-2.75,2.4,'a)')
plt.xlabel('LCN [z-score]')
plt.ylabel('Deciduous \n GPP [z-score]')
plt.subplot(222)                                
plt.tight_layout()
h = plt.hist2d(deci_sla,deci_g,vmin=0,vmax=45)
plt.text(-2.4,2.4,'b)')
plt.xlabel('SLA [z-score]')
plt.ylabel('GPP [z-score]')
plt.subplot(223)                                
plt.tight_layout()
h = plt.hist2d(ever_lcn,ever_g,vmin=0,vmax=45)
plt.text(-2.8,2.55,'c)')
plt.xlabel('LCN [z-score]')
plt.ylabel('Evergreen \n GPP [z-score]')
plt.subplot(224)                                
plt.tight_layout()
h = plt.hist2d(ever_sla,ever_g,vmin=0,vmax=45)
plt.text(-2.4,2.55,'d)')
plt.xlabel('SLA [z-score]')
plt.ylabel('GPP [z-score]')
# colorbar
plt.subplots_adjust(bottom=0.1, right=0.8, top=0.9)
cax = plt.axes([0.85, 0.1, 0.075, 0.8])
clb = plt.colorbar(cax=cax)
clb.ax.set_title('counts')
plt.savefig('phen_tr.pdf')


# maybe just use curve fits?
plt.figure(figsize=[7.5,7.5])
plt.rcParams.update({'font.size':16})
ax = plt.subplot(221)                                
plt.tight_layout()
# initial values
p0 = (np.mean(gppRu),np.mean(gppRu)/np.mean(lcnRu),1)
p_ru, pcov = curve_fit(nfunc,lcnRu,gppRu,p0=p0)
x = np.linspace(np.min(lcnRu),np.max(lcnRu),50)
plt.plot(x,nfunc(x,*p_ru),'b-')
plt.plot(lcnRu,gppRu,'k.',label='trait values')
plt.plot(np.mean(lcnRu),np.mean(gppRu),'ko',ms=10,label='dist. mean')
plt.plot(np.mean(lcnRu),gppcmRu,'r*',ms=10,label='updated mean')
plt.text(-11,1300,'a)')
plt.xlabel('Leaf CN [gC/gN]')    
plt.ylabel('GPP [gC m$^{-2}$ yr$^{-1}$]')
ax.legend(loc='upper left',bbox_to_anchor=[0.21,1.01],fontsize='x-small',frameon=False,shadow=False,fancybox=False)
plt.subplot(222)                                
plt.tight_layout()
p0 = (np.mean(gppFr),np.mean(gppFr)/np.mean(lcnFr),1)
p_fr, pcov = curve_fit(nfunc,lcnFr,gppFr,p0=p0)
x = np.linspace(np.min(lcnFr),np.max(lcnFr),50)
plt.plot(x,nfunc(x,*p_fr),'b-')
plt.plot(lcnFr,gppFr,'k.',label='trait values')
plt.plot(np.mean(lcnFr),np.mean(gppFr),'ko',ms=10,label='distribution mean')
plt.plot(np.mean(lcnFr),gppcmFr,'r*',ms=10,label='updated mean')
plt.text(-12,1430,'b)')
plt.xlabel('Leaf CN [gC/gN]')    
plt.ylabel('GPP [gC m$^{-2}$ yr$^{-1}$]')
plt.subplot(223)                                
plt.tight_layout()
p0 = (np.mean(gppBr),np.mean(gppBr)/np.mean(slaBr),1)
p_br, pcov = curve_fit(nfunc,slaBr,gppBr,p0=p0)
x = np.linspace(np.min(slaBr),np.max(slaBr),50)
plt.plot(x,nfunc(x,*p_br),'b-')
plt.plot(slaBr,gppBr,'k.')
plt.plot(np.mean(slaBr),np.mean(gppBr),'ko',ms=10)
plt.plot(np.mean(slaBr),gppcmBr,'r*',ms=10)
plt.text(-0.038,5400,'c)')
plt.xticks([0.01,0.04,0.07])
plt.xlabel('SLA [m$^2$/gC]')
plt.ylabel('GPP [gC m$^{-2}$ yr$^{-1}$]')
plt.subplot(224)                                
plt.tight_layout()
p0 = (np.mean(gppUs),np.mean(gppUs)/np.mean(slaUs),1)
p_us, pcov = curve_fit(nfunc,slaUs,gppUs,p0=p0)
x = np.linspace(np.min(slaUs),np.max(slaUs),50)
plt.plot(x,nfunc(x,*p_us),'b-')
plt.plot(slaUs,gppUs,'k.')
plt.plot(np.mean(slaUs),np.mean(gppUs),'ko',ms=10)
plt.plot(np.mean(slaUs),gppcmUs,'r*',ms=10)
plt.text(-0.017,2080,'d)')
plt.xlabel('SLA [m$^2$/gC]')
plt.ylabel('GPP [gC m$^{-2}$ yr$^{-1}$]')
plt.savefig('mdiff_nl_example.pdf')
plt.close('all')

#plt.figure(figsize=[7.5,7.5])
#plt.rcParams.update({'font.size':16})
#plt.subplot(221)                                
#plt.tight_layout()
#plt.plot(slaRu,gppRu,'k.')
#plt.text(-0.035,720,'a)')
#plt.xlabel('SLA [m$^2$/gC]')
#plt.ylabel('GPP [gC m$^{-2}$ yr$^{-1}$]')
#plt.subplot(222)                                
#plt.tight_layout()
#plt.plot(lcnRu,gppRu,'k.')
#plt.text(-4.5,720,'b)')
#plt.xlabel('Leaf CN [gC/gN]')    
#plt.ylabel('GPP [gC m$^{-2}$ yr$^{-1}$]')
#plt.subplot(223)                                
#plt.tight_layout()
#plt.plot(slaBr,gppBr,'k.')
#plt.text(-0.035,5400,'c)')
#plt.xlabel('SLA [m$^2$/gC]')
#plt.ylabel('GPP [gC m$^{-2}$ yr$^{-1}$]')
#plt.subplot(224)                                
#plt.tight_layout()
#plt.plot(lcnBr,gppBr,'k.')
#plt.text(-20,5400,'d)')
#plt.xlabel('Leaf CN [gC/gN]')    
#plt.ylabel('GPP [gC m$^{-2}$ yr$^{-1}$]')
#plt.savefig('nonlin_example.pdf')
#plt.close('all')

# is there a stronger signal when we divide into evergreen and deciduous groups?
# not really - too few points and boreal shrub outlier is too strong
#evergreen = [0,1,2,4,6,8,10]
#decid = [3,5,7,9,11,12,13,14]
#plt.plot(np.array(tra)[evergreen],np.array(gpp_pdiff)[evergreen],'g.')
#plt.plot(np.array(tra)[decid],np.array(gpp_pdiff)[decid],'k.')
#plt.show()

# consolidate into a small subset (4-6) of "types" of response based on dominant trait
# and nonlinearity
for i,s in enumerate(site):
    # look at GPP for 0 values
    vals = comb_runs[i].iloc[:,range(7,419,4)]
    # pick out zero values 
    pl = [iv != 0 for iv in vals.mean()[0:100]]
    # make the line below into a for loop...
    pl = [pl[iv] and incl.iloc[iv,site_PFT_num[i]-1] for iv in range(len(pl))]
    # if plotting something other than GPP (NPP = 8, MR = 9)
    vals = comb_runs[i].iloc[:,range(7,419,4)]
    gpp = vals.mean()[0:100]
    gppcm = vals.GPPcm.mean()
    sit_tr = tr_loc_assign(site_PFT_num[i]).iloc[0:100,:]
    sit_trd = tr_loc_assign(site_PFT_num[i]).iloc[101,:]
    gppd = vals.GPP101.mean()
    # SLA(LLS)-LAI-GPP
    plt.figure(figsize=[10,4])
    # SLA-LLS
    plt.subplot(131)
    plt.tight_layout()
    sla_sit = sit_tr.sla[pl]
    lls_sit = sit_tr.lls[pl]
    plt.plot(sla_sit,lls_sit,'k.')
    plt.plot(sla[i],lls[i],'ro')
    plt.plot(sit_trd.sla,sit_trd.lls,'g*')
    plt.xlabel('SLA [m$^2$/gC]')
    plt.ylabel('LLS [years]')
    # SLA-LAI
    plt.subplot(132)
    plt.tight_layout()
    lai = comb_runs[i].iloc[:,range(10,419,4)]
    lai_s = lai.mean()[0:100]
    laicm = lai.LAIcm.mean()
    laid = lai.LAI101.mean()
    plt.plot(sla_sit,lai_s[pl],'k.')
    plt.plot(sla[i],np.mean(lai_s[pl]),'ro')
    plt.plot(sla[i],laicm,'b*')
    plt.plot(sit_trd.sla,laid,'g*')
    plt.xlabel('SLA [m$^2$/gC]')
    plt.ylabel('LAI [m$^{2}$ m$^{-2}$]')
    plt.title(handle[i]+' '+s)
    # LAI-GPP
    plt.subplot(133)
    plt.tight_layout()
    plt.plot(lai_s[pl],gpp[pl],'k.')
    plt.plot(np.mean(lai_s[pl]),np.mean(gpp[pl]),'ro')
    plt.plot(laicm,gppcm,'b*')
    plt.plot(laid,gppd,'g*')
    plt.xlabel('LAI [m$^{2}$ m$^{-2}$]')
    plt.ylabel('GPP [gC m$^{-2}$ yr$^{-1}$]')
    plt.savefig(s+'_tr_lai_gpp.pdf')
    plt.close('all')
    # GPP-trait relationships
    plt.figure(figsize=[10,4])
    # SLA
    plt.subplot(131)
    plt.tight_layout()
    tra = sit_tr.sla[pl]
    # initial sla values
    # p0 = (np.mean(gpp[pl]),np.mean(gpp[pl])/np.mean(tra),1)
    p0 = (0,np.mean(gpp[pl])/np.mean(tra),1)
    p_sla, pcov = curve_fit(nfunc,tra,gpp[pl],p0=p0,ftol=1e-4)
    x = np.linspace(np.min(tra),np.max(tra),50)
    plt.plot(tra,gpp[pl],'k.')
    plt.plot(x,nfunc(x,*p_sla),'r-')
    plt.plot(sla[i],np.mean(gpp[pl]),'ro')
    plt.plot(sla[i],gppcm,'b*')
    plt.plot(sit_trd.sla,gppd,'g*')
    plt.xlabel('SLA [m$^2$/gC]')
    plt.ylabel('GPP [gC m$^{-2}$ yr$^{-1}$]')
    sla_r2 = np.corrcoef(nfunc(tra,*p_sla),gpp[pl])[0,1]**2
    plt.title('r$^2$= '+str(np.round(sla_r2,2))) 
    # LCN
    plt.subplot(132)
    plt.tight_layout()
    tra = sit_tr.lcn[pl]
    # initial lcn values
    # p0 = (np.mean(gpp[pl]),np.mean(gpp[pl])/np.mean(tra),1)
    p0 = (0,np.mean(gpp[pl])/np.mean(tra),1)
    p_lcn, pcov = curve_fit(nfunc,tra,gpp[pl],p0=p0,ftol=1e-3)
    x = np.linspace(np.min(tra),np.max(tra),50)
    plt.plot(tra,gpp[pl],'k.')
    plt.plot(x,nfunc(x,*p_lcn),'r-')
    plt.plot(lcn[i],np.mean(gpp[pl]),'ro')
    plt.plot(lcn[i],gppcm,'b*')
    plt.plot(sit_trd.lcn,gppd,'g*')
    plt.xlabel('Leaf CN [gC/gN]')    
    plt.ylabel('GPP [gC m$^{-2}$ yr$^{-1}$]')
    lcn_r2 = np.corrcoef(nfunc(tra,*p_lcn),gpp[pl])[0,1]**2
    plt.title(handle[i]+' r$^2$= '+str(np.round(lcn_r2,2)))
     # LLS
    plt.subplot(133)
    plt.tight_layout()
    tra = sit_tr.lls[pl]
    # initial lls values
    # p0 = (np.mean(gpp[pl]),np.mean(gpp[pl])/np.mean(tra),1)
    p0 = (0,np.mean(gpp[pl])/np.mean(tra),1)
    p_lls, pcov = curve_fit(nfunc,tra,gpp[pl],p0=p0,ftol=1e-4)
    x = np.linspace(np.min(tra),np.max(tra),50)
    plt.plot(tra,gpp[pl],'k.')
    plt.plot(x,nfunc(x,*p_lls),'r-')
    plt.plot(lls[i],np.mean(gpp[pl]),'ro')
    plt.plot(lls[i],gppcm,'b*')
    plt.plot(sit_trd.lls,gppd,'g*')
    plt.xlabel('Leaf Lifespan [years]')    
    plt.ylabel('GPP [gC m$^{-2}$ yr$^{-1}$]')
    lls_r2 = np.corrcoef(nfunc(tra,*p_lls),gpp[pl])[0,1]**2
    plt.title('r$^2$= '+str(np.round(lls_r2,2)))
    plt.savefig(s+'_trait_gpp_nl.pdf')
    plt.close('all')
    nparams.append((p_sla,p_lcn,p_lls))
    tr_r2.append((sla_r2,lcn_r2,lls_r2))
