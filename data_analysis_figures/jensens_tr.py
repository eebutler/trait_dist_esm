# cd Desktop/Projects/Plants/UMN/ACME_trait_distribution/scripts

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

f_comb = pd.read_pickle('../../script_input/flx_comb.p')

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
cm_out = pd.read_pickle('../../script_input/yrly_cm.p')
yr_comb = list()
for i in range(15):
    yr_comb.append(f_comb[i].join(cm_out[i],rsuffix='cm'))
# trait values
# Latin Hypercube trait values
tr = pd.read_csv('../../script_input/lhd_100.csv')
# defaults from 1000 sample random 
tr_d = pd.read_csv('../../script_input/subPFT_1000.csv')
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
# re-run with 0:-2, strip off mean trait and default
mm = 102 # 100 for distribution mean, 101 for model default,102 for constrained mean 
tmp = list()
mag_diff = list()
for i,y in enumerate(yr_comb):
    ym = y.iloc[:,range(7,419,4)].mean()
    pl = [yi != 0 for yi in ym]
    pl = np.array([pl[iv] and incl.iloc[iv,site_PFT_num[i]-1] for iv in range(len(pl))])
    my = np.array(ym)[pl][0:-2].mean()
    tmp.append((my-ym[mm])/ym[mm])
    mag_diff.append(my-ym[mm])

per_diff = np.round(tmp,decimals=3)*100

site = ['AU-Tum','BR-Sa1','FI-Hyy','FR-Pue','GF-Guy','RU-Cok','RU-Fyo',
	'RU-Sam','US-Blo','US-Ha1','US-PFa','US-SRM','US-UMB','US-Wkg',
	'ZA-Kru']


PFT = ['BET','BET','NET','BDT,80:NET,20','BET','BDS',
       'NET','BDS','NET','BDT,85:NET,15','BDT,50:NET,50',
       'BDS','BDT','C4','C4,70:BDT,30']

handle = ['Tree-Te$_2$','Tree-Tr$_2$','Tree-Bo$_1$','Tree*-Te$_3$',
	  'Tree-Tr$_1$','Shrub-Bo$_1$','Tree-Bo$_2$','Shrub-Bo$_2$',
	  'Tree-Te$_3$','Tree*-Te$_2$','Tree*-Te$_1$','Shrub-Ar',
	  'Tree-Te$_1$','Grass-Ar','Grass*-Ar']

# quick look
# perdiff_df = pd.DataFrame([list(per_diff),PFT]).transpose()
# perdiff_df.index = site

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

p_ord = np.argsort(t2m)

per_diff = [per_diff[p] for p in p_ord]
mag_diff = [mag_diff[p] for p in p_ord]
PFT = [PFT[p] for p in p_ord]
t2m = [t2m[p] for p in p_ord]
pre = [pre[p] for p in p_ord]
pet = [pet[p] for p in p_ord]

t2m = np.array(t2m).round(decimals=1)
t2m_s = map(str,t2m)
handlep = [handle[p] for p in p_ord]

# percent difference
plt.bar(range(15),per_diff,0.5,tick_label=handlep)
#plt.text(-2.5,55,'a)')
plt.ylim(-10,40) 
plt.xticks(rotation=60)
#plt.ylabel('Percent Deviation of E[f(X)]-f(E[X])')
plt.ylabel('Distribution Mean - Updated Mean [%]')
plt.subplots_adjust(bottom=0.2)
plt.savefig('per_diff_Tcm.pdf')
plt.close()

# magnitude difference
plt.bar(range(15),mag_diff,0.5,tick_label=handlep)
plt.text(-2,375,'b)')
plt.ylim(-120,375) 
plt.xticks(rotation=60)
#plt.ylabel('E[f(X)]-f(E[X]) [gC m$^{-2}$ yr$^{-1}$]')
plt.ylabel('Distribution Mean - Updated Mean [gC m$^{-2}$ yr$^{-1}$]')
plt.subplots_adjust(bottom=0.2)
plt.savefig('mag_diff_Tcm.pdf')
plt.close()
