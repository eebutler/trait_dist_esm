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
plt.savefig('siteGPPcm_only.pdf')
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

fdd = np.array(fm) - np.array(ymd)
fd_err = np.sqrt((fdd**2).mean())
fd_errs = np.sqrt((fdd**2).std())
fd_errm = np.median(np.sqrt((fdd**2)))


fmdd = np.array(fm) - np.array(my)
fmd_err = np.sqrt((fmdd**2).mean())
fmd_errs = np.sqrt((fmdd**2).std())
fmd_errm = np.median(np.sqrt((fmdd**2)))

fcd = np.array(fm) - np.array(myc)
fc_err = np.sqrt((fcd**2).mean())
fc_errs = np.sqrt((fcd**2).std())
fc_errm = np.median(np.sqrt((fcd**2)))

fmcd = np.array(fm) - np.array(ymc)
fmc_err = np.sqrt((fmcd**2).mean())
fmc_errs = np.sqrt((fmcd**2).std())
fmc_errm = np.median(np.sqrt((fmcd**2)))

# kruskal-wallis test
stats.kruskal(fdd,fcd,fmcd)

np.array([fd_err, fc_err,fmc_err]).round()
np.array([fd_errs, fc_errs,fmc_errs]).round()
np.array([fm_errm, fd_errm, fmd_errm, fc_errm]).round()

# kruskal-wallis across years
# nans are a problem
#y = yr_comb[2]
pk = list()

for i,y in enumerate(yr_comb):
    ym = y.iloc[:,range(7,419,4)].mean()
    pl = [yi != 0 for yi in ym]
    pl = np.array([pl[iv] and incl.iloc[iv,site_PFT_num[i]-1] for iv in range(len(pl))])
     
    yg = y.iloc[:,range(7,419,4)]
    yc = yg.iloc[:,pl[0:100]] # censored
    yu = yg.iloc[:,0:100] # uncensored
    ym = yg.iloc[:,100] # input mean
    yd = yg.iloc[:,101] # default
    yic = yg.iloc[:,102] # censored input mean
    yf = yr_comb[i].GPP_NT_CUT_REF # flux tower
    kstat,p = stats.kruskal(yc.values.ravel(),yu.values.ravel(),yd.values,ym.values,yic.values,
			    yf.values,nan_policy = 'omit')
    pk.append(p)




# trait combinations that gave rise to GPP values
tr_close = list()
tr_max = list()
sit_cl_rnk = list()
sit_mx_rnk = list()
for i in range(15):
    site_loc = (site_PFT_num[i]-1)*3
    tr_close.append(tr.iloc[gpp_diff[i],site_loc:site_loc+3].values)
    tr_max.append(tr.iloc[gpp_max[i],site_loc:site_loc+3].values)
    cl_rnk = list()
    mx_rnk = list()
    for j in range(3):
	tmp = tr.iloc[:,site_loc+j].sort_values(kind='mergesort')
	cl_rnk.append(tmp.index.get_loc(gpp_diff[i]))
	mx_rnk.append(tmp.index.get_loc(gpp_max[i]))
    sit_cl_rnk.append(cl_rnk)
    sit_mx_rnk.append(mx_rnk)

# https://stackoverflow.com/questions/40425840/least-squares-fit-in-python-for-3d-surface

def tr3dplot(tr,name,label,rank=False,optim=False):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    
    tr = np.array(tr)
    
    fc = ['g','b','c','m','brown']
    labels = ['Broadleaf Deciduous','Broadleaf Evergreen',
	      'Needleleaf Evergreen','Shrub','Grass','Mixed']
    for i in range(5):
	x = tr[pft_l[i]][:,0]
	y = tr[pft_l[i]][:,1]
	z = tr[pft_l[i]][:,2]
	for j in range(len(x)): 
		ax.plot3D([x[j],x[j]],[y[j],y[j]],[0,z[j]],color=fc[i])
	ax.scatter(tr[pft_l[i]][:,0],tr[pft_l[i]][:,1],tr[pft_l[i]][:,2],c=fc[i],
		   label=labels[i],depthshade=False)
     
    ax.scatter(tr[pft_l[5]][:,0],tr[pft_l[5]][:,1],tr[pft_l[5]][:,2],
	       facecolor='None',edgecolors='black',depthshade=False,s=45,lw=1.5,label=labels[5])
    
    ax.legend(loc='upper right')

    ax.set_xlabel('SLA [m$^2$ kgC$^{-1}$]')
    ax.set_ylabel('LCN [gC gN$^{-1}$]')
    ax.set_zlabel('LLS [years]')
    # max GPP limits
    ax.set_xlim(0,0.05)
    ax.set_ylim(0,40)
    ax.set_zlim(0,10)
    if optim:
	ax.set_xlim(0,0.06)   
    	ax.set_ylim(0,70)
    	ax.set_zlim(0,10)
    if rank:
	ax.set_xlabel('SLA')
	ax.set_ylabel('LCN')
    	ax.set_zlabel('LLS')
	ax.set_xlim(0,100)
    	ax.set_ylim(0,100)
    	ax.set_zlim(0,100)
    	ax.text(0,-5,120,label,fontsize=12)
    else:
	ax.text(0,-1,12,label,fontsize=12)
    plt.savefig(name+'.pdf')
    plt.close()

#tr3dplot(tr_close,'closest_trait','a)',False,True)
tr3dplot(tr_max,'gppmax_trait','a)')
#tr3dplot(sit_cl_rnk,'closest_trait_rank','b)',True,True)
##jitter
#tr3dplot(np.array(sit_mx_rnk)+np.random.normal(0,0.75,[15,3]),'gppmax_trait_rank','b)',True)

# full trait space
x = tr.iloc[:,range(0,42,3)].values
y = tr.iloc[:,range(1,42,3)].values
z = tr.iloc[:,range(2,42,3)].values

pl = list()
for i in range(15):
    ym = yr_comb[i].iloc[:,range(7,415,4)].mean()
    m0 = [yi != 0 for yi in ym]
    pl.append(np.array([m0[iv] and incl.iloc[iv,site_PFT_num[i]-1] for iv in range(len(m0))]))

xp = x[:,np.array(site_PFT_num)-1]
yp = y[:,np.array(site_PFT_num)-1]
zp = z[:,np.array(site_PFT_num)-1]

xpc = np.concatenate([xp[pl[i],i] for i in range(15)])
ypc = np.concatenate([yp[pl[i],i] for i in range(15)])
zpc = np.concatenate([zp[pl[i],i] for i in range(15)])

log_tr = np.log10(np.array([xpc,ypc,zpc])).T
log_trm = np.min(log_tr,0)
log_trx = np.max(log_tr,0)
# gpp max
log_trmax = np.log10(np.array(tr_max))
log_trxm = np.min(log_trmax,0)
log_trxx = np.max(log_trmax,0)

# set up figure
fig = plt.figure()
ax = fig.gca(projection='3d')
ax.set_aspect("equal")

# from https://codereview.stackexchange.com/questions/155585/plotting-a-rectangular-prism
# draw cube
def rect_prism(x_range, y_range, z_range,c):
    xx, yy = np.meshgrid(x_range, y_range)
    ax.plot_wireframe(xx, yy, z_range[0], color=c)
    ax.plot_surface(xx, yy, z_range[0], color=c, alpha=0.2)
    ax.plot_wireframe(xx, yy, z_range[1], color=c)
    ax.plot_surface(xx, yy, z_range[1], color=c, alpha=0.2)

    yy, zz = np.meshgrid(y_range, z_range)
    ax.plot_wireframe(x_range[0], yy, zz, color=c)
    ax.plot_surface(x_range[0], yy, zz, color=c, alpha=0.2)
    ax.plot_wireframe(x_range[1], yy, zz, color=c)
    ax.plot_surface(x_range[1], yy, zz, color=c, alpha=0.2)

    xx, zz = np.meshgrid(x_range, z_range)
    ax.plot_wireframe(xx, y_range[0], zz, color=c)
    ax.plot_surface(xx, y_range[0], zz, color=c, alpha=0.2)
    ax.plot_wireframe(xx, y_range[1], zz, color=c)
    ax.plot_surface(xx, y_range[1], zz, color=c, alpha=0.2)

    ax.set_xlabel('log$_{10}$(SLA [m$^2$ kgC$^{-1}$])') 
    ax.set_ylabel('log$_{10}$(LCN [gC gN$^{-1}$])')
    ax.set_zlabel('log$_{10}$(LLS [years])')

    #ax.set_xlim([-1,1])
    #ax.set_ylim([-1.5,1.6])
    #ax.set_zlim([-2,1])

rect_prism(np.array([log_tr[:,0].min(), log_tr[:,0].max()]), 
	   np.array([log_tr[:,1].min(), log_tr[:,1].max()]), 
	   np.array([log_tr[:,2].min(), log_tr[:,2].max()]),"r")
rect_prism(np.array([log_trmax[:,0].min(), log_trmax[:,0].max()]), 
	   np.array([log_trmax[:,1].min(), log_trmax[:,1].max()]), 
	   np.array([log_trmax[:,2].min(), log_trmax[:,2].max()]),"b")

xmin = np.repeat(log_tr[:,0].min(),log_tr.shape[0])
ymax = np.repeat(log_tr[:,1].max(),log_tr.shape[0])
zmin = np.repeat(log_tr[:,2].min(),log_tr.shape[0])

ax.scatter3D(xmin,log_tr[:,1],log_tr[:,2],color='grey',depthshade=False,zorder=1)
ax.scatter3D(log_tr[:,0],ymax,log_tr[:,2],color='grey',depthshade=False,zorder=1)
#ax.scatter3D(log_tr[:,0],log_tr[:,1],zmin,color='grey',depthshade=False,zorder=1)

trmx = np.log10(np.array(tr_max))
x = trmx[:,0]
y = trmx[:,1]
z = trmx[:,2]
for j in range(len(x)): 
	ax.plot3D([x[j],x[j]],[y[j],y[j]],[z.min(),z[j]],color='k',zorder=2)
ax.scatter(trmx[:,0],trmx[:,1],trmx[:,2],c='k',depthshade=False,zorder=2)

ax.text(-2.5,0.6,2.5,'b)')
#plt.show()
plt.savefig('log_niches.pdf')

# trait space volumes
all_dims = 10**(log_trx)-10**(log_trm) 
all_dimsx = 10**(log_trxx)-10**(log_trxm)
vol_ratio = (all_dimsx[0]*all_dimsx[1]*all_dimsx[2])/(all_dims[0]*all_dims[1]*all_dims[2])


# regression fits
tr_fit = list()
trmax_fit = list()
ltr = list()
ltxr = list()
mx_x0 = list()
mi_x0 = list()
xmx_x0 = list()
xmi_x0 = list()
resid_tr = list()
resid_trx = list()
for i in range(3):
    x = [0,0,1]
    y = [1,2,2]
    tr_fit.append(np.polyfit(log_tr[:,x[i]],log_tr[:,y[i]],1))
    ltr.append(np.linspace(log_trm[i],log_trx[i]))
    res_tr = log_tr[:,y[i]]-(tr_fit[i][1]+tr_fit[i][0]*log_tr[:,x[i]])
    rmx = np.argmax(res_tr)
    rmi = np.argmin(res_tr)
    resid_tr.append([res_tr.min(),res_tr.max()])
    mx_x0.append(-1*(tr_fit[i][0]*log_tr[rmx,x[i]]-log_tr[rmx,y[i]]))
    mi_x0.append(-1*(tr_fit[i][0]*log_tr[rmi,x[i]]-log_tr[rmi,y[i]]))
    
    trmax_fit.append(np.polyfit(log_trmax[:,x[i]],log_trmax[:,y[i]],1))
    ltxr.append(np.linspace(log_trxm[i],log_trxx[i]))
    res_trx = log_trmax[:,y[i]]-(trmax_fit[i][1]+trmax_fit[i][0]*log_trmax[:,x[i]])
    rmx = np.argmax(res_trx)
    rmi = np.argmin(res_trx)
    resid_trx.append([res_trx.min(),res_trx.max()])
    xmx_x0.append(-1*(trmax_fit[i][0]*log_trmax[rmx,x[i]]-log_trmax[rmx,y[i]]))
    xmi_x0.append(-1*(trmax_fit[i][0]*log_trmax[rmi,x[i]]-log_trmax[rmi,y[i]]))

# set up figure
fig = plt.figure()
ax = fig.gca(projection='3d')
ax.set_aspect("equal")

rect_prism(np.array(resid_tr[0]), 
	   np.array(resid_tr[1]), 
	   np.array(resid_tr[2]),"r")
rect_prism(np.array(resid_trx[0]), 
	   np.array(resid_trx[1]), 
	   np.array(resid_trx[2]),"b")

xmin = np.repeat(resid_tr[:,0].min(),resid_tr.shape[0])
ymax = np.repeat(resid_tr[:,1].max(),resid_tr.shape[0])
zmin = np.repeat(resid_tr[:,2].min(),resid_tr.shape[0])

ax.scatter3D(xmin,resid_tr[:,1],resid_tr[:,2],color='grey')
ax.scatter3D(resid_tr[:,0],ymax,resid_tr[:,2],color='grey')
ax.scatter3D(resid_tr[:,0],resid_tr[:,1],zmin,color='grey')

# 2D plots
# sla-lcn
plt.plot(log_tr[:,0],log_tr[:,1],'.',color='grey')
plt.plot(log_trmax[:,0],log_trmax[:,1],'ko')
plt.plot(ltr[0],tr_fit[0][1]+tr_fit[0][0]*ltr[0],'r')
plt.plot(ltr[0],mx_x0[0]+tr_fit[0][0]*ltr[0],'r--')
plt.plot(ltr[0],mi_x0[0]+tr_fit[0][0]*ltr[0],'r--')
plt.plot(ltxr[0],trmax_fit[0][1]+trmax_fit[0][0]*ltxr[0],'b')
plt.plot(ltxr[0],xmx_x0[0]+trmax_fit[0][0]*ltxr[0],'b--')
plt.plot(ltxr[0],xmi_x0[0]+trmax_fit[0][0]*ltxr[0],'b--')

plt.xlabel('log$_{10}$(SLA [m$^2$ kgC$^{-1}$])')
plt.ylabel('log$_{10}$(LCN [gC gN$^{-1}$])')
plt.savefig('sla_lcn_gppmax_log.pdf')
plt.close()
# sla-lls
plt.plot(log_tr[:,0],log_tr[:,2],'.',color='grey')
plt.plot(log_trmax[:,0],log_trmax[:,2],'ko')
plt.xlabel('log$_{10}$(SLA [m$^2$ kgC$^{-1}$])')
plt.ylabel('log$_{10}$(LLS [years])')
plt.savefig('sla_lls_gppmax_log.pdf')
plt.close()
# lcn-lls
plt.plot(log_tr[:,1],log_tr[:,2],'.',color='grey')
plt.plot(log_trmax[:,1],log_trmax[:,2],'ko')
plt.xlabel('log$_{10}$(LCN [gC gN$^{-1}$])')
plt.ylabel('log$_{10}$(LLS [years])')
plt.savefig('lcn_lls_gppmax_log.pdf')
plt.close()



