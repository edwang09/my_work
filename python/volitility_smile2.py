"""
Created on Tue Feb 16 23:21:45 2016

@author: jwang67
"""
import numpy as np
import pandas as pd
from math import log, sqrt, exp
from scipy import stats
from scipy.optimize import fsolve
import mpl_toolkits.mplot3d.axes3d as p3
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rcParams['font.family'] = 'serif'
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
from scipy.interpolate import griddata
from scipy.integrate import quad
def BSM_call_value(S0, div, K, T, r, sigma):
    ''' Valuation of European call option in BSM Model.
    --> Analytical Formula.
    Parameters
    ==========
    S0: float
    initial stock/index level
    K: float
    strike price
    T: float
    time-to-maturity (for t=0)
    r: float
    constant risk-free short rate
    sigma: float
    volatility factor in diffusion term
    Returns
    =======
    call_value: float
    European call option present value
    '''
    d1 = (np.log(S0 / K) + (r -div+ 0.5 * sigma ** 2) * T)/ (sigma * np.sqrt(T))
    d2 = (np.log(S0 / K) + (r-div - 0.5 * sigma ** 2) * T)/ (sigma * np.sqrt(T))
    BS_C = (S0 * exp(-div * (T))* stats.norm.cdf(d1, 0.0, 1.0) - K * exp(-r * T) * stats.norm.cdf(d2, 0.0, 1.0))

    return BS_C
def imp_vol(r,div,S0,K,m,p, sigma_est=0.2):
    ''' Return implied volatility given option price. '''
    def difference(sigma):
        ssigma = sigma       
        value =BSM_call_value(S0,div, K, m, r, ssigma)
        return value- p
    iv = fsolve(difference, sigma_est)[0]
    return iv
def BSM_call_value_INT(S0, div, K, T, r, sigma):
    ''' Valuation of European call option in BSM model via Lewis (2001)
    --> Fourier-based approach (integral).
    Parameters
    ==========
    S0: float
    initial stock/index level
    K: float
    strike price
    T: float
    time-to-maturity (for t=0)
    r: float
    constant risk-free short rate
    sigma: float
    volatility factor in diffusion term
    Returns
    =======
    call_value: float
    European call option present value'''
    int_value = quad(lambda u:BSM_integral_function(u, S0, K, T, r-div, sigma), 0, 100)[0]
    call_value = max(0, S0*np.exp(-div * T) - np.exp(-r * T) * np.sqrt(S0 * K)/ np.pi * int_value)
    return call_value
def BSM_integral_function(u, S0, K, T, r, sigma):
    ''' Valuation of European call option in BSM model via Lewis (2001)
    --> Fourier-based approach: integral function. '''
    cf_value = BSM_characteristic_function(u - 1j * 0.5, 0.0, T, r, sigma)
    int_value = 1 / (u ** 2 + 0.25) * (np.exp(1j * u * np.log(S0 / K)) * cf_value).real
    return int_value
def BSM_characteristic_function(v, x0, T, r, sigma):
    ''' Valuation of European call option in BSM model via
    Lewis (2001) and Carr-Madan (1999)
    --> Fourier-based approach: characteristic function. '''
    cf_value = np.exp(((x0 / T + r - 0.5 * sigma ** 2) * 1j * v
    - 0.5 * sigma ** 2 * v ** 2) * T)
    return cf_value
    
#
# define parameters
#
    
    
r=0.001
div=0.02
S0=1972.18









xl = pd.ExcelFile("spxcalls20150831.xlsx")
df = xl.parse('WRDS')
tm=df['expiration']
t=df['date']
m=list()
for i in range(517):
    m.append((tm[i]-t[i]).days/365.0)
k=df['strike']
p=df['price']
v=list()
for i in range(517):
    v.append(imp_vol(r,div,S0,k[i],m[i],p[i]))
m=np.asarray(m)
x1 = np.linspace(k.min(), k.max(), )
y1 = np.linspace(m.min(), m.max(), )
x2, y2 = np.meshgrid(x1, y1)
z2 = griddata((k, m), v, (x2, y2), method='cubic')
fig = plt.figure()
ax = fig.gca(projection='3d')
surf = ax.plot_surface(x2, y2, z2, rstride=1, cstride=1, cmap=cm.jet, linewidth=0, antialiased=False)
ax.set_zlim(0.12, 0.35)
ax.zaxis.set_major_locator(LinearLocator(10))
ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))
fig.colorbar(surf, shrink=0.5, aspect=5)
plewis=list()
pbs=list()
for i in range(517):
    plewis.append(BSM_call_value_INT(S0,div, k[i], m[i], r, v[i]))
df['plewis']=plewis
writer = pd.ExcelWriter('spxcalls20150831.xlsx')
df.to_excel(writer,'WRDS')
writer.save()