# cd /Users/timaeus/Desktop/Projects/Plants/UMN/ACME_trait_distribution/scripts/supp_figs

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors
from mpl_toolkits.basemap import Basemap

fol = '~/Desktop/Projects/Papers_InProgress/trait-dist_esm/site_information/'

sites = pd.read_excel(fol+'15_sites_PFT_years.xlsx',skip_footer=16)

pft = pd.read_excel(fol+'15_sites_PFT_years.xlsx',sheetname=5)

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

# narrow pft species in TRY subset
pft_num = [59,53,5,1883,540,298,323,172,157,91,62,81,3823,267]
pnum = np.log(np.array(pft_num))
pnum = pnum/pnum.mean()
msize = [pnum[4],pnum[3],pnum[1],(pnum[6]*0.8+pnum[0]*0.2),pnum[3],
	 pnum[10],pnum[1],pnum[10],pnum[0],(pnum[6]*0.85+pnum[0]*0.15),
	 (pnum[0]*0.5+pnum[6]*0.5),pnum[9],pnum[6],pnum[13],
	 (pnum[13]*0.7+pnum[5]*0.3)]


# need to jitter the lat of one of the AZ sites to make it more visible
#sites.loc[11,'lat'] = sites.loc[11,'lat']+2
#sites.loc[11,'lon'] = sites.loc[11,'lon']-2

m = Basemap(llcrnrlon=-130,llcrnrlat=-45,urcrnrlon=155,urcrnrlat=77,\
                   projection='cyl',resolution='c')

m.drawcoastlines(linewidth=0.7)
m.drawparallels(np.arange(-90.,120.,30.),labels=[1,0,0,0],linewidth=0.3)
m.drawmeridians(np.arange(0.,360.,60.),labels=[0,0,0,1],linewidth=0.3)

# Broadleaf Deciduous
m.plot(sites.lon[3],sites.lat[3],'g.',ms=24*msize[3])
m.plot(sites.lon[9],sites.lat[9],'g.',ms=24*msize[9])
m.plot(sites.lon[12],sites.lat[12],'g.',ms=24*msize[12])
m.plot(sites.lon[3],sites.lat[3],'.',ms=24*msize[3],mfc='None',mec='white',lw=1)
m.plot(sites.lon[9],sites.lat[9],'.',ms=24*msize[9],mfc='None',mec='white',lw=1)
m.plot(sites.lon[12],sites.lat[12],'.',ms=24*msize[12],mfc='None',mec='white',lw=1)
# Broadleaf Evergreen
m.plot(sites.lon[0],sites.lat[0],'b.',ms=24*msize[0])
m.plot(sites.lon[1],sites.lat[1],'b.',ms=24*msize[1])
m.plot(sites.lon[4],sites.lat[4],'b.',ms=24*msize[4])
m.plot(sites.lon[0],sites.lat[0],'b.',ms=24*msize[0],mfc='None',mec='white',lw=1)
m.plot(sites.lon[1],sites.lat[1],'b.',ms=24*msize[1],mfc='None',mec='white',lw=1)
m.plot(sites.lon[4],sites.lat[4],'b.',ms=24*msize[4],mfc='None',mec='white',lw=1)
# Needleleaf Evergreen
m.plot(sites.lon[2],sites.lat[2],'c.',ms=24*msize[2])
m.plot(sites.lon[6],sites.lat[6],'c.',ms=24*msize[6])
m.plot(sites.lon[8],sites.lat[8],'c.',ms=24*msize[8])
m.plot(sites.lon[10],sites.lat[10],'c.',ms=24*msize[10])
m.plot(sites.lon[2],sites.lat[2],'c.',ms=24*msize[2],mfc='None',mec='white',lw=1)
m.plot(sites.lon[6],sites.lat[6],'c.',ms=24*msize[6],mfc='None',mec='white',lw=1)
m.plot(sites.lon[8],sites.lat[8],'c.',ms=24*msize[8],mfc='None',mec='white',lw=1)
m.plot(sites.lon[10],sites.lat[10],'c.',ms=24*msize[10],mfc='None',mec='white',lw=1)
# Grass
m.plot(sites.lon[13],sites.lat[13],'m.',ms=24*msize[13])
m.plot(sites.lon[14],sites.lat[14],'m.',ms=24*msize[14])
m.plot(sites.lon[13],sites.lat[13],'m.',ms=24*msize[13],mfc='None',mec='white',lw=1)
m.plot(sites.lon[14],sites.lat[14],'m.',ms=24*msize[14],mfc='None',mec='white',lw=1)
# Shrub
m.plot(sites.lon[5],sites.lat[5],'.',color='brown',ms=24*msize[5])
m.plot(sites.lon[7],sites.lat[7],'.',color='brown',ms=24*msize[7])
m.plot(sites.lon[11]-2.5,sites.lat[11],'.',color='brown',ms=24*msize[11])
m.plot(sites.lon[5],sites.lat[5],'.',color='brown',ms=24*msize[5],mfc='None',mec='white',lw=1)
m.plot(sites.lon[7],sites.lat[7],'.',color='brown',ms=24*msize[7],mfc='None',mec='white',lw=1)
m.plot(sites.lon[11]-2.5,sites.lat[11],'.',color='brown',ms=24*msize[11],mfc='None',mec='white',lw=1)
# Mixed
m.plot(sites.lon[3],sites.lat[3],'.',ms=24*msize[3],mfc='None',mec='black',lw=1)#,label='Mixed')
m.plot(sites.lon[9],sites.lat[9],'.',ms=24*msize[9],mfc='None',mec='black',lw=1)#,label='Mixed')
m.plot(sites.lon[10],sites.lat[10],'.',ms=24*msize[10],mfc='None',mec='black',lw=1)#,label='Mixed')
m.plot(sites.lon[14],sites.lat[14],'.',ms=24*msize[14],mfc='None',mec='black',lw=1)#,label='Mixed')

m.plot([],[],'g.',ms=24,label='Broadleaf Deciduous')
m.plot([],[],'b.',ms=24,label='Broadleaf Evergreen')
m.plot([],[],'c.',ms=24,label='Needleleaf Evergreen')
m.plot([],[],'.',color='brown',ms=24,label='Shrub')
m.plot([],[],'m.',ms=24,label='C4 Grass')
m.plot([],[],'.',ms=24,mfc='None',mec='black',lw=1,label='Mixed')

# PFTs as background
lat = pd.read_csv('supp_figs/pft_lat.csv').values
lon = pd.read_csv('supp_figs/pft_lon.csv').values
# need to convert PFTs to the five categories present + "other"
dom = pd.read_csv('supp_figs/domPFT.csv').values
# more elegant solution to consolidate PFTs?
dom[np.where(dom==1)] = 1
dom[np.where(dom==2)] = 2
dom[np.where(dom==3)] = 2
dom[np.where(dom==4)] = 99
dom[np.where(dom==5)] = 3
dom[np.where(dom==6)] = 3
dom[np.where(dom==7)] = 4
dom[np.where(dom==8)] = 4
dom[np.where(dom==9)] = 4
dom[np.where(dom==10)] = 5
dom[np.where(dom==11)] = 5
dom[np.where(dom==12)] = 5
dom[np.where(dom==13)] = 99
dom[np.where(dom==14)] = 99
dom[np.where(dom==15)] = 6
dom[np.where(dom==16)] = 99
dom[np.where(dom==99)] = 7

cmap = matplotlib.colors.ListedColormap(['white','c','b','g','brown','m','grey'])

x,y = m(lon,lat)
cs = m.pcolormesh(x,y,dom,cmap=cmap) # need custom cmap to match PFTs

#m.plot(sites.lon[bdt],sites.lat[bdt],'g.',ms=24,label='Broadleaf Deciduous')
#m.plot(sites.lon[bet],sites.lat[bet],'b.',ms=24,label='Broadleaf Evergreen')
#m.plot(sites.lon[net],sites.lat[net],'c.',ms=24,label='Needleleaf Evergreen')
#m.plot(sites.lon[shr],sites.lat[shr],'.',color='brown',ms=24,label='Shrub')
#m.plot(sites.lon[gra],sites.lat[gra],'m.',ms=24,label='C4 Grass')
#m.plot(sites.lon[mix],sites.lat[mix],'.',ms=24,mfc='None',mec='black',lw=1,label='Mixed')

plt.legend(bbox_to_anchor=(0.,-0.5,1.,0.),loc=8,ncol=2)

lon_offset = [0,0,0,0,-2,-2,0,0,0,0,-30,0,-2,2,2]
lat_offset = [0,0,0,0,0,-5,0,0,0,0,0,0,2.5,-5,-5]
for i in range(sites.shape[0]):
    plt.annotate(s=handle[i],xy=(sites.lon[i]+lon_offset[i],sites.lat[i]+lat_offset[i]),
		 fontsize=8,fontweight='bold')
plt.title('Site Locations and Dominant PFT')

plt.savefig('site_loc_pft_lab.pdf',bbox='tight')

#PFT1 
#59
#PFT2
#53 
#PFT3
#5
#PFT4
#1883
#PFT5
#540
#PFT6
#298
#PFT7
#323
#PFT8
#172
#PFT9
#157
#PFT10
#91
#PFT11
#62
#PFT12
#81
#PFT13
#3823
#PFT14
#267
