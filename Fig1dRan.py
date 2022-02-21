#
# Fit beta distribution to a single dataset. Tutorial example.
#

import numpy as np
#from sebcolour import Colour as C
#aclr = C.violetred   # Colour for 'adjacent' data
#rclr = C.dodgerblue  # Colour for 'random' data

# Set plotting font defaults
import matplotlib
fs = 18
fnt = {'family' : 'Arial',
       'weight' : 'regular',
       'size'   : fs}
matplotlib.rc('font', **fnt)

# Important for svg output of text as 'things that can be edited in inkscape'
import pylab as pl
pl.rcParams['svg.fonttype'] = 'none'

nbins = 80       # Number of bins in all histos
binw = 2.0/nbins # Bin width in all histos

import scipy
from scipy.optimize import curve_fit

#
# Here are three possible functions to use in the curve fitting
#

# 'full' because I fit the beta distribution to all the data
#
# Compute the beta distribution after transforming x. x is an array of real
# values in the open range (-1,1). It is first transformed to the range (0,1)
# to compute the beta distribution (which can be asym), which is returned to
# the caller.
def beta_distribution_full(x, a, b):
    x_ = (x+1)*0.5
    y_ = (1./scipy.special.beta(a, b)) * np.power(x_, a-1) * np.power(1-x_, b-1) * 0.5
    return y_

# Symmetric version of the above in which b=a
def beta_distribution_full_sym(x, a):
    x_ = (x+1)*0.5
    y_ = (1./scipy.special.beta(a, a)) * np.power(x_, a-1) * np.power(1-x_, a-1) * 0.5
    return y_

# 'half' because I fit the beta distribution to the RHS of the data and then
# mirror.
#
# An alternative kind of beta distribution fit is to try fitting a
# non-symmetric beta distribution to *half* of the real data (the right hand
# side) and then mirror the fit for the left hand side.
def beta_distribution_half(x, a, b):
    x_ = x[int(len(x)/2):] # right hand side
    betadist_rhs = (1./scipy.special.beta(a, b)) * np.power(x_, a-1) * np.power(1-x_, b-1)
    y_ = np.hstack((np.flip(betadist_rhs),betadist_rhs)) # join rhs and lhs
    return y_

logpath = "./logsMorph"
# Load histogram data
hdata = np.genfromtxt(logpath + '/random_correlate1.data')

# Control overall figure size here
fig = pl.figure(figsize=(8,8))

# Histogram our data to produce histogrammed values, h and bin edges, be.
h, be = np.histogram(hdata, bins=nbins,  density=True)
# Compute an offset to turn bin edges into histogram x values
be0 = 0.5*(be[1]-be[0])
# And apply the offset to create x
x = be[:-1] + be0
print(x)
# h*binw turns h (which is a prob dens. func) into 'probability' - h_prob.
h_prob = h * binw
# Create an axis
ax = fig.add_subplot(1,1,1)
# Plot the histogram data bar graph
ax.bar(x, h_prob, width=binw, color='r')
# This is the curve fitting. Choose form beta_distribution_full_sym, beta_distribution_full or beta_distribution_half here:
# we try twice with different intial parameters then select according to sos
# initial guess for the parameters
iparams1 = np.array([3.0])
popt1, pcov1 = curve_fit(beta_distribution_full_sym, x, h, iparams1,  method='lm')
#popt1[0] = 0.1
# Note I transform the fitted vales from PDF to probability to be compared with h_prob:
h_fit1 = 2.0*beta_distribution_full_sym(x, *popt1) * binw
# Compute fit quality
normh_fit1 = np.sum(np.power(h_fit1,2))
fit_sos1 = np.sum(np.power((h_fit1 - h_prob), 2)*beta_distribution_full_sym(x,iparams1[0]))
# initial guess for the parameters
iparams2 = np.array([1.0])
popt2, pcov2 = curve_fit(beta_distribution_full_sym, x, h, iparams2,  method='lm')
#popt2[0] = 1.5
# Note I transform the fitted vales from PDF to probability to be compared with h_prob:
h_fit2 = 2.0*beta_distribution_full_sym(x, *popt2) * binw * 0.5
# Compute fit quality
normh_fit2 = np.sum(np.power(h_fit2,2))
fit_sos2 = np.sum(np.power((h_fit2 - h_prob), 2)*beta_distribution_full_sym(x,iparams2[0]))

print(fit_sos1, " , " , fit_sos2)
# Plot the fitted values
#if fit_sos1 <= fit_sos2:
#if fit_sos1 <= fit_sos2:
#ax.plot(x, h_fit1, linestyle='--', color='b', label='sym, {0}={1}={2:.3f}; sos={3:.3f}'.format(r'$\alpha$', r'$\beta$', popt1[0], fit_sos1))
#else:
ax.plot(x, h_fit2, linestyle='--', color='b', label='sym, {0}={1}={2:.3f}; sos={3:.3f}'.format(r'$\alpha$', r'$\beta$', popt2[0], fit_sos2))
# Display the legend, which shows the fit params
ax.legend(fontsize=12,loc='upper center',frameon=False)

# Fix the ticks, set the labels
ax.set_xticks([-1,0,1])
ax.set_xticklabels(['-1','0','1'])
ax.set_yticks([0,0.1,0.2])
ax.set_yticklabels(['0','.1','.2'])
ax.set_xlabel('c')
ax.set_ylabel('p')

pl.show()
#pl.savefig(logpath + "/Fig1eRan.png")


