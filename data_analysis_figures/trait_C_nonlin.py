#cd /Users/timaeus/Desktop/Projects/Plants/UMN/ACME_trait_distribution

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

usu = pd.read_pickle('../script_input/USU_all.p')[0]

for i,u in enumerate(usu):
    usu[i] = u.loc['2000':'2014']

# divide into random and latin hypercube
usur = usu[0:1002]

tr = pd.read_csv('../script_input/subPFT_1000.csv')

tr = tr.loc[:,['sla.6','lnm.6','lls.6']]

ugpp = [u.GPP.mean() for u in usur]
unpp = [u.NPP.mean() for u in usur]
umr = [u.MR.mean() for u in usur]
ulai = [u.LAI.mean() for u in usur]

ugpp = np.array(ugpp)
nz = np.where(ugpp!=0)[0]
plt.plot(tr.loc[0:999,'sla.6'][nz],ugpp[nz],'k.')
plt.rcParams.update({'font.size':20})
#plt.plot(tr.loc[1000,'sla.6'],ugpp[1000],'k.')
#plt.plot(tr.loc[1001,'sla.6'],ugpp[1001],'r.')
plt.xlabel('SLA [m$^2$/gC]')
plt.ylabel('GPP [gC m$^{-2}$ yr$^{-1}$]')
plt.text(-0.05,1275, 'a)')
#plt.title('GPP Non-linearity with SLA')
plt.savefig('gpp_sla.pdf',bbox_inches='tight',pad_inches=0.1)
plt.close()

plt.plot(tr.loc[0:999,'lnm.6'][nz],ugpp[nz],'k.')
plt.rcParams.update({'font.size':20})
#plt.plot(tr.loc[1000,'lnm.6'],ugpp[1000],'k.')
#plt.plot(tr.loc[1001,'lnm.6'],ugpp[1001],'r.')
plt.xlabel('LCN [gC/gN]')
plt.ylabel('GPP [gC m$^{-2}$ yr$^{-1}$]')
plt.text(-14,1275,'b)')
#plt.title('GPP Non-linearity with leaf CN ratio')
plt.savefig('gpp_lcn.pdf',bbox_inches='tight',pad_inches=0.1)
plt.close()


plt.plot(tr.loc[0:999,'sla.6'],unpp[0:1000],'.',color='0.8')
plt.plot(tr.loc[1000,'sla.6'],unpp[1000],'k.')
plt.plot(tr.loc[1001,'sla.6'],unpp[1001],'r.')
plt.xlabel('SLA [m^2/gC]')
plt.ylabel('NPP [gC m^-2 yr^-1]')
plt.title('NPP Non-linearity with SLA')
plt.savefig('npp_sla.pdf')
plt.close()

plt.plot(tr.loc[0:999,'lnm.6'],unpp[0:1000],'.',color='0.8')
plt.plot(tr.loc[1000,'lnm.6'],unpp[1000],'k.')
plt.plot(tr.loc[1001,'lnm.6'],unpp[1001],'r.')
plt.xlabel('LCN [gC/gN]')
plt.ylabel('NPP [gC m^-2 yr^-1]')
plt.title('NPP Non-linearity with leaf CN ratio')
plt.savefig('npp_lcn.pdf')
plt.close()

plt.plot(tr.loc[0:999,'sla.6'],umr[0:1000],'.',color='0.8')
plt.plot(tr.loc[1000,'sla.6'],umr[1000],'k.')
plt.plot(tr.loc[1001,'sla.6'],umr[1001],'r.')
plt.xlabel('SLA [m^2/gC]')
plt.ylabel('MR [gC m^-2 yr^-1]')
plt.title('MR Non-linearity with SLA')
plt.savefig('mr_sla.pdf')
plt.close()

plt.plot(tr.loc[0:999,'lnm.6'],umr[0:1000],'.',color='0.8')
plt.plot(tr.loc[1000,'lnm.6'],umr[1000],'k.')
plt.plot(tr.loc[1001,'lnm.6'],umr[1001],'r.')
plt.xlabel('LCN [gC/gN]')
plt.ylabel('MR [gC m^-2 yr^-1]')
plt.title('MR Non-linearity with leaf CN ratio')
plt.savefig('mr_lcn.pdf')
plt.close()

plt.plot(tr.loc[0:999,'sla.6'],ulai[0:1000],'.',color='0.8')
plt.plot(tr.loc[1000,'sla.6'],ulai[1000],'k.')
plt.plot(tr.loc[1001,'sla.6'],ulai[1001],'r.')
plt.xlabel('SLA [m^2/gC]')
plt.ylabel('LAI [gC m^-2 yr^-1]')
plt.title('LAI Non-linearity with SLA')
plt.savefig('lai_sla.pdf')
plt.close()

plt.plot(tr.loc[0:999,'lnm.6'],ulai[0:1000],'.',color='0.8')
plt.plot(tr.loc[1000,'lnm.6'],ulai[1000],'k.')
plt.plot(tr.loc[1001,'lnm.6'],ulai[1001],'r.')
plt.xlabel('LCN [gC/gN]')
plt.ylabel('LAI [gC m^-2 yr^-1]')
plt.title('LAI Non-linearity with leaf CN ratio')
plt.savefig('lai_lcn.pdf')
plt.close()

### Brazil, Santarem site

brs = pd.read_pickle('../script_input/BRS_all.p')[0]
    
for i,b in enumerate(brs):
    brs[i] = b.loc['2002':'2011']

# divide into random and latin hypercube
brsr = brs[0:1002]

tr = pd.read_csv('../script_input/subPFT_1000.csv')

tr = tr.loc[:,['sla.3','lnm.3','lls.3']]

bgpp = [b.GPP.mean() for b in brsr]
bnpp = [b.NPP.mean() for b in brsr]
bmr = [b.MR.mean() for b in brsr]
blai = [b.LAI.mean() for b in brsr]

bgpp = np.array(bgpp[0:1000])
nz = np.where(bgpp!=0)[0]

plt.plot(tr.loc[0:999,'sla.3'][nz],bgpp[nz],'k.')
plt.rcParams.update({'font.size':20})
#plt.plot(tr.loc[1000,'sla.3'],bgpp[1000],'k.')
#plt.plot(tr.loc[1001,'sla.3'],bgpp[1001],'r.')
plt.xlabel('SLA [m$^2$/gC]')
plt.ylabel('GPP [gC m$^{-2}$ yr${^-1}$]')
plt.text(-0.0275,6475,'c)')
#plt.title('GPP Non-linearity with SLA')
plt.savefig('gpp_sla_brs.pdf',bbox_inches='tight',pad_inches=0.1)
plt.close()

plt.plot(tr.loc[0:999,'lnm.3'][nz],bgpp[nz],'k.')
plt.rcParams.update({'font.size':20})
#plt.plot(tr.loc[1000,'lnm.3'],bgpp[1000],'k.')
#plt.plot(tr.loc[1001,'lnm.3'],bgpp[1001],'r.')
plt.xlabel('LCN [gC/gN]')
plt.ylabel('GPP [gC m$^{-2}$ yr$^{-1}$]')
plt.text(-28,6475,'d)')
#plt.title('GPP Non-linearity with leaf CN ratio')
plt.savefig('gpp_lcn_brs.pdf',bbox_inches='tight',pad_inches=0.1)
plt.close()

plt.plot(tr.loc[0:999,'lls.3'],bgpp[0:1000],'.',color='0.8')
plt.plot(tr.loc[1000,'lls.3'],bgpp[1000],'k.')
plt.plot(tr.loc[1001,'lls.3'],bgpp[1001],'r.')
plt.xlabel('LLS [years]')
plt.ylabel('GPP [gC m^-2 yr^-1]')
plt.title('GPP Non-linearity with leaf lifespan')
plt.savefig('gpp_lls_brs.pdf')
plt.close()

plt.plot(tr.loc[0:999,'sla.3'],bnpp[0:1000],'.',color='0.8')
plt.plot(tr.loc[1000,'sla.3'],bnpp[1000],'k.')
plt.plot(tr.loc[1001,'sla.3'],bnpp[1001],'r.')
plt.xlabel('SLA [m^2/gC]')
plt.ylabel('NPP [gC m^-2 yr^-1]')
plt.title('NPP Non-linearity with SLA')
plt.savefig('npp_sla_brs.pdf')
plt.close()

plt.plot(tr.loc[0:999,'lnm.3'],bnpp[0:1000],'.',color='0.8')
plt.plot(tr.loc[1000,'lnm.3'],bnpp[1000],'k.')
plt.plot(tr.loc[1001,'lnm.3'],bnpp[1001],'r.')
plt.xlabel('LCN [gC/gN]')
plt.ylabel('NPP [gC m^-2 yr^-1]')
plt.title('NPP Non-linearity with leaf CN ratio')
plt.savefig('npp_lcn_brs.pdf')
plt.close()

plt.plot(tr.loc[0:999,'lls.3'],bnpp[0:1000],'.',color='0.8')
plt.plot(tr.loc[1000,'lls.3'],bnpp[1000],'k.')
plt.plot(tr.loc[1001,'lls.3'],bnpp[1001],'r.')
plt.xlabel('LLS [years]')
plt.ylabel('NPP [gC m^-2 yr^-1]')
plt.title('NPP Non-linearity with leaf lifespan')
plt.savefig('npp_lls_brs.pdf')
plt.close()

plt.plot(tr.loc[0:999,'sla.3'],bmr[0:1000],'.',color='0.8')
plt.plot(tr.loc[1000,'sla.3'],bmr[1000],'k.')
plt.plot(tr.loc[1001,'sla.3'],bmr[1001],'r.')
plt.xlabel('SLA [m^2/gC]')
plt.ylabel('MR [gC m^-2 yr^-1]')
plt.title('MR Non-linearity with SLA')
plt.savefig('mr_sla_brs.pdf')
plt.close()

plt.plot(tr.loc[0:999,'lnm.3'],bmr[0:1000],'.',color='0.8')
plt.plot(tr.loc[1000,'lnm.3'],bmr[1000],'k.')
plt.plot(tr.loc[1001,'lnm.3'],bmr[1001],'r.')
plt.xlabel('LCN [gC/gN]')
plt.ylabel('MR [gC m^-2 yr^-1]')
plt.title('MR Non-linearity with leaf CN ratio')
plt.savefig('mr_lcn_brs.pdf')
plt.close()

plt.plot(tr.loc[0:999,'lls.3'],bmr[0:1000],'.',color='0.8')
plt.plot(tr.loc[1000,'lls.3'],bmr[1000],'k.')
plt.plot(tr.loc[1001,'lls.3'],bmr[1001],'r.')
plt.xlabel('LLS [years]')
plt.ylabel('MR [gC m^-2 yr^-1]')
plt.title('MR Non-linearity with leaf lifespan')
plt.savefig('mr_lls_brs.pdf')
plt.close()

plt.plot(tr.loc[0:999,'sla.3'],blai[0:1000],'.',color='0.8')
plt.plot(tr.loc[1000,'sla.3'],blai[1000],'k.')
plt.plot(tr.loc[1001,'sla.3'],blai[1001],'r.')
plt.xlabel('SLA [m^2/gC]')
plt.ylabel('LAI [gC m^-2 yr^-1]')
plt.title('LAI Non-linearity with SLA')
plt.savefig('lai_sla_brs.pdf')
plt.close()

plt.plot(tr.loc[0:999,'lnm.3'],blai[0:1000],'.',color='0.8')
plt.plot(tr.loc[1000,'lnm.3'],blai[1000],'k.')
plt.plot(tr.loc[1001,'lnm.3'],blai[1001],'r.')
plt.xlabel('LCN [gC/gN]')
plt.ylabel('LAI [gC m^-2 yr^-1]')
plt.title('LAI Non-linearity with leaf CN ratio')
plt.savefig('lai_lcn_brs.pdf')
plt.close()

plt.plot(tr.loc[0:999,'lls.3'],blai[0:1000],'.',color='0.8')
plt.plot(tr.loc[1000,'lls.3'],blai[1000],'k.')
plt.plot(tr.loc[1001,'lls.3'],blai[1001],'r.')
plt.xlabel('LLS [years]')
plt.ylabel('LAI [gC m^-2 yr^-1]')
plt.title('LAI Non-linearity with leaf lifespan')
plt.savefig('lai_lls_brs.pdf')
plt.close()
