# cd Desktop/Projects/Plants/UMN/ACME_trait_distribution/scripts
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import scipy.stats as stats
from mpl_toolkits.mplot3d import Axes3D

f_comb = pd.read_pickle('../script_input/flx_comb.p')

t2m = [y.T2m.mean() for y in f_comb]
pre = [y.PRE.mean() for y in f_comb]
ev = [y.Evap.mean() for y in f_comb]
tr = [y.Tran.mean() for y in f_comb]

pet = np.array(ev)+np.array(tr)

#comb_runs = pd.read_pickle('flx_comb.p')
# manually add in NaN values for bad years
f_comb[1].loc['2007','GPP_NT_CUT_REF'] = np.nan
f_comb[10].loc['1995','GPP_NT_CUT_REF'] = np.nan
# append constrained mean output
cm_out = pd.read_pickle('../script_input/yrly_cm.p')
yr_comb = list()
for i in range(15):
    yr_comb.append(f_comb[i].join(cm_out[i],rsuffix='cm'))
# trait values
# Latin Hypercube trait values
tr = pd.read_csv('../script_input/lhd_100.csv')
# defaults from 1000 sample random 
tr_d = pd.read_csv('../script_input/subPFT_1000.csv')
# append defaults
tr = tr.append(tr_d.iloc[1001,:])
tr.index = range(0,102)
# cannot track traits in this format for constrained input mean (unique to site)
# this may mess with some of the code below (omit from this file)

# evaluate after omitting large leaf lifespan values
# 90th percentile, data max, and data max + 1 sd
lls90 = [12.6,8.1,0.5,4.5,3.4,4.6,0.8,0.5,3.0,0.8,0.5,0.2,1.5,0.5]
lls_max = [24.0,8.2,0.5,7.9,6.8,15.5,0.8,0.5,3.4,0.9,0.5,0.2,3.0,0.5]
lls_omax = [26.4,10.1,1.1,10.0,8.5,18.0,2.1,2.0,5.1,2.3,1.8,1.5,5.2,2.2]

incl = []
# set up omission dataframe
lls_lim = lls_max
o_0 = tr.loc[:,'lls'].values < lls_lim[0]
o_13 = [tr.loc[:,'lls.'+str(i)].values < lls_lim[i] for i in range(1,14)]
for i in range(14):
    if i == 0:
	incl.append(o_0)
    else:
	incl.append(o_13[i-1])

incl = pd.DataFrame(incl)

incl = incl.transpose()
# add False for constrained mean, only to keep consistent size of ND array
incl = incl.append([np.repeat(False,14)])
incl.index = range(103)

site_PFT_num = [5,4,2,7,4,11,2,12,1,7,1,10,7,14,14]

# difference in means for GPP between f(E[x]) and E[f(X)]
# run with range 0:-2, strip off mean trait and default
# mean values
ymm = list() # uncensored
ymmst = list() # " std
ymc = list() # censored
ymcst = list() # " std
ymd = list() # default value
ymdst = list() # " std
myc = list() # output, censored
mycst = list() # ", " std
my = list() # output, uncensored
myst = list() # ", " std
fm = list() # flux tower
fmst = list() # ", std
flux_diff = list() # test diversity influence on divergence from flux tower
gpp_diff = list()
gpp_max = list()
cnum = list()
for i,y in enumerate(yr_comb):
    ym = y.iloc[:,range(7,419,4)].mean()
    ymst = y.iloc[:,range(7,419,4)].std()
    pl = [yi != 0 for yi in ym]
    pl = np.array([pl[iv] and incl.iloc[iv,site_PFT_num[i]-1] for iv in range(len(pl))])
    cnum.append(pl.sum())
    myc.append(np.array(ym)[pl][0:-3].mean())
    mycst.append(np.array(ymst)[pl][0:-3].mean())
    my.append(np.array(ym)[0:-3].mean())
    myst.append(np.array(ymst)[0:-3].mean())
    ymm.append(ym[100])
    ymmst.append(ymst[100])
    ymd.append(ym[101])
    ymdst.append(ymst[101])
    ymc.append(ym[102])
    ymcst.append(ymst[102])
    flux_mean = np.nanmean(y.GPP_NT_CUT_REF)
    flux_std = np.nanstd(y.GPP_NT_CUT_REF)
    fm.append(flux_mean)
    fmst.append(flux_std)
    gpp_diff.append(np.argmin(np.abs(ym.values-flux_mean)))
    gpp_max.append(np.argmax(ym.values))
    flux_diff.append(flux_mean - myc[i])

site = ['AU-Tum','BR-Sa1','FI-Hyy','FR-Pue','GF-Guy','RU-Cok','RU-Fyo',
	'RU-Sam','US-Blo','US-Ha1','US-PFa','US-SRM','US-UMB','US-Wkg',
	'ZA-Kru']

PFT = ['BET','BET','NET','BDT*','BET','BDS',
       'NET','BDS','NET','BDT*','BDT*',
       'BDS','BDT','C4','C4*']

handle = ['Tree-Te$_2$','Tree-Tr$_2$','Tree-Bo$_1$','Tree*-Te$_3$',
	  'Tree-Tr$_1$','Shrub-Bo$_1$','Tree-Bo$_2$','Shrub-Bo$_2$',
	  'Tree-Te$_3$','Tree*-Te$_2$','Tree*-Te$_1$','Shrub-Ar',
	  'Tree-Te$_1$','Grass-Ar','Grass*-Ar']

# broadleaf deciduous tree indices
bdt = [3,9,12]
# broadleaf evergreen tree
bet = [0,1,4]
# needleleaf evergreen tree
net = [2,6,8,10]
# shrubs
shr = [5,7,11]
# grass
gra = [13,14]
# mixed
mix = [3,9,10,14]

# joint list
pft_l = [bdt,bet,net,shr,gra,mix]

# bar plot
p_ord = np.argsort(t2m)

ymd = [ymd[p] for p in p_ord]
ymdst = [ymdst[p] for p in p_ord]
ymm = [ymm[p] for p in p_ord]
ymmst = [ymmst[p] for p in p_ord]
ymc = [ymc[p] for p in p_ord]
ymcst = [ymcst[p] for p in p_ord]
myc = [myc[p] for p in p_ord]
mycst = [mycst[p] for p in p_ord]
my = [my[p] for p in p_ord]
myst = [myst[p] for p in p_ord]
fm = [fm[p] for p in p_ord]
fmst = [fmst[p] for p in p_ord]
handlep = [handle[p] for p in p_ord]

width = 0.3
x = np.array(range(15))
# add yerr to bar plots - temporal deviation
p0 = plt.bar(x-width,ymc,width,color='black',label='Updated Mean',yerr=ymcst)
p2 = plt.bar(x,myc,width,tick_label=handlep,
	     color='green',label='Distribution Mean',yerr=mycst)
p4 = plt.bar(x+width,ymd,width,color='red',label='Default Mean',yerr=ymdst)
for i in x:
    if i==0:
	p5 = plt.plot([i-0.35,i+0.45],[fm[i],fm[i]],'b-',lw=2,label='Flux Tower Data')
	plt.plot([i-0.35,i+0.45],[fm[i]+fmst[i],fm[i]+fmst[i]],'b:',lw=2)
	plt.plot([i-0.35,i+0.45],[fm[i]-fmst[i],fm[i]-fmst[i]],'b:',lw=2)
    else:
	plt.plot([i-0.35,i+0.45],[fm[i],fm[i]],'b-',lw=2)
	plt.plot([i-0.35,i+0.45],[fm[i]+fmst[i],fm[i]+fmst[i]],'b:',lw=2)
	plt.plot([i-0.35,i+0.45],[fm[i]-fmst[i],fm[i]-fmst[i]],'b:',lw=2)
lh = [p0,p2,p4,p5[0]]
plt.legend(lh,[H.get_label() for H in lh])
plt.xticks(rotation=60)
plt.ylabel('GPP [gC m$^{-2}$ yr$^{-1}$]')
plt.subplots_adjust(bottom=0.2)
#plt.savefig('siteGPPcm_only.pdf')
plt.close()

## including uncensored values
#width = 0.18
#x = np.array(range(15))
## add yerr to bar plots - temporal deviation
#p0 = plt.bar(x-width*2,ymc,width,color='magenta',label='Censored Input Mean',yerr=ymcst)
#p1 = plt.bar(x,my,width,color='cyan',label='Output Mean',yerr=myst)
#p2 = plt.bar(x-width,myc,width,tick_label=handlep,
#	     color='green',label='Censored Output Mean',yerr=mycst)
#p3 = plt.bar(x+width,ymm,width,color='black',label='Input Mean',yerr=ymmst)
#p4 = plt.bar(x+width*2,ymd,width,color='red',label='Default',yerr=ymdst)
#for i in x:
#    if i==0:
#	p5 = plt.plot([i-0.4,i+0.45],[fm[i],fm[i]],'b-',lw=2,label='Flux Tower Data')
#	plt.plot([i-0.4,i+0.45],[fm[i]+fmst[i],fm[i]+fmst[i]],'b:',lw=2)
#	plt.plot([i-0.4,i+0.45],[fm[i]-fmst[i],fm[i]-fmst[i]],'b:',lw=2)
#    else:
#	plt.plot([i-0.4,i+0.45],[fm[i],fm[i]],'b-',lw=2)
#	plt.plot([i-0.4,i+0.45],[fm[i]+fmst[i],fm[i]+fmst[i]],'b:',lw=2)
#	plt.plot([i-0.4,i+0.45],[fm[i]-fmst[i],fm[i]-fmst[i]],'b:',lw=2)
#lh = [p0,p1,p2,p3,p4,p5[0]]
#plt.legend(lh,[H.get_label() for H in lh])
#plt.xticks(rotation=60)
#plt.ylabel('GPP [gC m$^{-2}$ yr$^{-1}$]')
#plt.subplots_adjust(bottom=0.2)
#plt.savefig('siteGPPcm.pdf')
#plt.close()

# check difference in means
fmd = np.array(fm) - np.array(ymm)
fm_err = np.sqrt((fmd**2).mean())
fm_errs = np.sqrt((fmd**2).std())
fm_errm = np.median(np.sqrt((fmd**2)))

# default
fdd = np.array(fm) - np.array(ymd)
fd_err = np.sqrt((fdd**2).mean())
fd_errs = np.sqrt((fdd**2).std())
fd_errm = np.median(np.sqrt((fdd**2)))

fmdd = np.array(fm) - np.array(my)
fmd_err = np.sqrt((fmdd**2).mean())
fmd_errs = np.sqrt((fmdd**2).std())
fmd_errm = np.median(np.sqrt((fmdd**2)))

# censored output mean
fcd = np.array(fm) - np.array(myc)
fc_err = np.sqrt((fcd**2).mean())
fc_errs = np.sqrt((fcd**2).std())
fc_errm = np.median(np.sqrt((fcd**2)))

# censored input mean
fmcd = np.array(fm) - np.array(ymc)
fmc_err = np.sqrt((fmcd**2).mean())
fmc_errs = np.sqrt((fmcd**2).std())
fmc_errm = np.median(np.sqrt((fmcd**2)))

# kruskal-wallis test
stats.kruskal(fdd,fcd,fmcd)

np.array([fd_err, fc_err,fmc_err]).round()
np.array([fd_errs, fc_errs,fmc_errs]).round()
np.array([fm_errm, fd_errm, fmd_errm, fc_errm]).round()

# abs(pd.DataFrame([fdd,fcd,fmcd])).idxmin()

tower_DF = pd.DataFrame([fdd,fcd,fmcd]).transpose()
tower_DF.index = handlep
tower_DF.columns = ['Default','Output','Input']
