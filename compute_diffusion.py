import numpy as np

##########################################################################
# Computing length in Angstroms
lower  = 1.8
upper  = 2.8
length = upper-lower
##########################################################################

##########################################################################
# Computing residence time in fs
# Reading the survival probability curve
ndata = 12500
t     = np.zeros(ndata)
prob  = np.zeros(ndata)

with open('survival_probability.dat') as input_file:
    for i in range(ndata):
        t[i],prob[i] = [t for t in next(input_file).split()]

# Eq. (13)
# Trapezoidal rule
residence_time = prob[0]
for i in range(1,ndata-1):
    residence_time += 2.0*prob[i]
residence_time = 0.5*(t[1]-t[0])*(residence_time+prob[ndata-1])
##########################################################################

##########################################################################
# Computing alpha
# Reading density 
maxbin = 1000
z      = np.zeros(maxbin)
dens   = np.zeros(maxbin)

with open('density.dat') as input_file:
    for i in range(maxbin):
        z[i],dens[i] = [z for z in next(input_file).split()]

rz    = np.zeros(maxbin)
pk    = np.zeros(maxbin)
k     = 0

for i in range(maxbin):
    if (z[i] >= lower and z[i] <= upper):
        rz[k] = z[i]
        pk[k] = np.log(dens[i])
        k    += 1

omega, b = np.polyfit(rz, pk, 1)

# Eq. (27)
om_len   = omega*length
sum_term = 0.0
for j in range(100):
    sum_term += 1.0/((2.0*j+1.0)**4.0*np.pi**4.0              \
                +0.75*om_len**2.0*(2.0*j+1)**2.0*np.pi**2.0   \
                -0.25*om_len**4.0)
alpha = 0.25/(om_len*(np.exp(om_len)+1.0)/(np.exp(om_len)-1.0)*sum_term)
##########################################################################

##########################################################################
# Computing diffusion coefficient in m²/s
# Eq. (15)
diffusion = 1e-5*length**2.0/residence_time/alpha

print('D = ',diffusion,' m²/s')
##########################################################################

