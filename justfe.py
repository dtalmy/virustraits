from numpy import *
from pylab import *
from scipy import *
from ehuxbasic import *
from basicforshow import *

#####################################################
# load in the data
#####################################################

# sample times
fehtimes = r_[[0,1.3,2.2,2.9,3.96,4.96,5.96]]*24
fevtimes = r_[[0,7]]*24

# hosts
hferep = r_[[541460,916550,1216350,1080625,715825,255125,162225]]
hfelim = r_[[515810,733725,839775 ,847875 ,713300,469600,370550]]

hferepsd = mean(r_[[0,11694.8,30875.1417,97687.62,34796.5,43829.81,9944.5]])
hfelimsd = mean(r_[[0,23867.33,8778.2,11244.2,33339.62,46452.82,30189.18]])

# virus
vferep = r_[[5414600,14427851991]]
vfelim = r_[[5158100,5700802151]]

vferepsd = mean(r_[[0,6432926774]])
vfelimsd = mean(r_[[0,893663897.9]])

#####################################################
# set up figures
#####################################################
f1,ax1 = subplots(1,2,figsize=[9.5,4.5])
f2,ax2 = subplots(4,1,figsize=[8,15])
f1.subplots_adjust(bottom=0.13,wspace=0.3,hspace=0.3)
f2.subplots_adjust(hspace=0.4)

#######################################################
# model fitting
#######################################################

## set up first guess params, MHA parameters, etc. 
phi = 2.97463146661e-09
muh = 1.47141995554
lam = 1.0053557416e-01
beta = 2843.1129

params = r_[[phi,muh,lam,beta]]
stds = zeros(4) + .1
opt = r_[[1,1,1,1]]
names = ['phi','muh','lam','beta']
params = log(params)
npars = params.shape[0]

nits = 100000
pits = 10000
burnin = 10000

##################################
# Fe replete
##################################

# initial conditions
inits = r_[[hferep[0],0,vferep[0]]]
m,n = fehtimes.shape[0], fevtimes.shape[0]
hnt,vnt = zeros(m),zeros(n)

# first run just to get error
ehuxbasic(params,inits,hnt,vnt,fehtimes,fevtimes,True)
chi =   sum((hnt - hferep) ** 2 / (hferepsd ** 2)) \
      + sum((vnt - vferep) ** 2 / (vferepsd ** 2))        

# distribution arrays and acceptance ratios
ar = 0.0
ars = r_[[]]
phis,muhs,lams,betas = r_[[]],r_[[]],r_[[]],r_[[]]
pall = [phis,muhs,lams,betas]

# run mha algorithm
for it in arange(1,nits,1):
        parsnew = params + opt*normal(0,stds,npars)
	sus,inf,vir = hferep[0],0,vferep[0]
	inits = r_[[sus,inf,vir]]
        ehuxbasic(parsnew,inits,hnt,vnt,fehtimes,fevtimes,True)
	chinew = sum((hnt - hferep) ** 2 / (hferepsd ** 2)) \
                +sum((vnt - vferep) ** 2 / (vferepsd ** 2))       
	if exp(chi-chinew) > rand():
                chi = chinew
                params = parsnew
                if it > burnin:
                        pall = append(pall,params[:,None],1)
                ar = ar + 1.0
        if (it % pits == 0):
                print it,chi,ar/pits
                ars = append(ars,ar/pits)
                ar = 0.0

print 'Optimal parameters for Fe replete conditions'
pars = r_[[ mean(p) for p in pall]]
for (p,l) in zip(pars,names):
        print l,'=',exp(p)

print ' '
print 'Standard deviations'
for (p,l) in zip(pall,names):
        print l+'std','=',std(exp(p))
print ' '

# redefine times for nicer looking plots
delt = 900.0 / 86400.0
ftimes = arange(0,1.2*amax(fehtimes)/24.0,delt)*24.0
n = ftimes.shape[0]
hnt,vnt = zeros(n),zeros(n)

# run again
sus,inf,vir = hferep[0],0,vferep[0]
inits = r_[[sus,inf,vir]]
forshow(pars,inits,hnt,vnt,ftimes,delt,True)

# plot 207
ax1[0].errorbar(fehtimes,hferep/1e+6,yerr=hferepsd/1e+6,c='r',fmt='o',label='Fe replete')
ax1[1].errorbar(fevtimes,vferep/1e+10,yerr=vferepsd/1e+10,c='r',fmt='o',label='Fe replete')
ax1[0].plot(ftimes,hnt/1e+6,c='b',lw=1.5,label='model fit')
ax1[1].plot(ftimes,vnt/1e+10,c='b',lw=1.5,label='model fit')

# Fe replete
ax2[0].hist(exp(pall[1]),label='Fe replete')
ax2[1].hist(exp(pall[0]),label='Fe replete')
ax2[2].hist(exp(pall[2]),label='Fe replete')
ax2[3].hist(exp(pall[3]),label='Fe replete')

##################################
# Fe limited
##################################

# initial conditions
inits = r_[[hfelim[0],0,vfelim[0]]]
m,n = fehtimes.shape[0], fevtimes.shape[0]
hnt,vnt = zeros(m),zeros(n)

# first run just to get error
ehuxbasic(params,inits,hnt,vnt,fehtimes,fevtimes,True)
chi =   sum((hnt - hfelim) ** 2 / (hfelimsd ** 2)) \
      + sum((vnt - vfelim) ** 2 / (vfelimsd ** 2))        

# distribution arrays and acceptance ratios
ar = 0.0
ars = r_[[]]
phis,muhs,lams,betas = r_[[]],r_[[]],r_[[]],r_[[]]
pall = [phis,muhs,lams,betas]

# run mha algorithm
for it in arange(1,nits,1):
        parsnew = params + opt*normal(0,stds,npars)
	sus,inf,vir = hfelim[0],0,vfelim[0]
	inits = r_[[sus,inf,vir]]
        ehuxbasic(parsnew,inits,hnt,vnt,fehtimes,fevtimes,True)
	chinew = sum((hnt - hfelim) ** 2 / (hfelimsd ** 2)) \
                +sum((vnt - vfelim) ** 2 / (vfelimsd ** 2))       
	if exp(chi-chinew) > rand():
                chi = chinew
                params = parsnew
                if it > burnin:
                        pall = append(pall,params[:,None],1)
                ar = ar + 1.0
        if (it % pits == 0):
                print it,chi,ar/pits
                ars = append(ars,ar/pits)
                ar = 0.0

print 'Optimal parameters for Fe limited conditions'
pars = r_[[ mean(p) for p in pall]]
for (p,l) in zip(pars,names):
        print l,'=',exp(p)

print ' '
print 'Standard deviations'
for (p,l) in zip(pall,names):
        print l+'std','=',std(exp(p))
print ' '

# redefine times for nicer looking plots
delt = 900.0 / 86400.0
ftimes = arange(0,1.2*amax(fehtimes)/24.0,delt)*24.0
n = ftimes.shape[0]
hnt,vnt = zeros(n),zeros(n)

# run again
sus,inf,vir = hfelim[0],0,vfelim[0]
inits = r_[[sus,inf,vir]]
forshow(pars,inits,hnt,vnt,ftimes,delt,True)

# plot 207
ax1[0].errorbar(fehtimes,hfelim/1e+6,yerr=hfelimsd/1e+6,c='gold',fmt='o',label='Fe limited')
ax1[1].errorbar(fevtimes,vfelim/1e+10,yerr=vfelimsd/1e+10,c='gold',fmt='o',label='Fe limited')
ax1[0].plot(ftimes,hnt/1e+6,c='g',lw=1.5,label='model fit')
ax1[1].plot(ftimes,vnt/1e+10,c='g',lw=1.5,label='model fit')

# Fe replete
ax2[0].hist(exp(pall[1]),label='Fe limited')
ax2[1].hist(exp(pall[0]),label='Fe limited')
ax2[2].hist(exp(pall[2]),label='Fe limited')
ax2[3].hist(exp(pall[3]),label='Fe limited')

# other axes bits
ax1[0].set_ylabel(r'Cells ($\times$10$^6$ml$^{-1}$)')
ax1[1].set_ylabel(r'Particles ($\times$10$^{10}$ ml$^{-1}$)')
ax1[0].set_xlabel('Time (hours)')
ax1[1].set_xlabel('Time (hours)')

ax2[0].set_xlabel(r'Host growth rate (day $^{-1}$)')
ax2[1].set_xlabel(r'Rate of infection (ml day$^{-1}$)')
ax2[2].set_xlabel(r'lysis rate (day $^{-1}$)')
ax2[3].set_xlabel(r'Burst size')

ax2[0].set_ylabel(r'Frequency')
ax2[1].set_ylabel(r'Frequency')
ax2[2].set_ylabel(r'Frequency')
ax2[3].set_ylabel(r'Frequency')

ax1[1].set_ylim([0,2])

for ax in ax1.flatten():
	ax.set_xlim([0,200])

for ax in ax2:
	ax.semilogx()

l1 = ax1[0].legend(loc='upper right',prop={'size':12})
l1.draw_frame(False)
l3 = ax2[0].legend(loc='upper left',prop={'size':12})
l3.draw_frame(False)

ax1[0].text(0.07,0.9,'a',ha='center',va='center',color='k',transform=ax1[0].transAxes)
ax1[1].text(0.07,0.9,'b',ha='center',va='center',color='k',transform=ax1[1].transAxes)

f1.savefig('figures/fe_dynamics')
f2.savefig('figures/fe_params')

show()
