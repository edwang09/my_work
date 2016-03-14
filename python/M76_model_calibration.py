
"""
Created on Tue Feb 16 16:30:43 2016

@author: jwang67
"""
#
# This is a code to calibration the Merton Jump diffusion model to Lewis Fourier-based aproach as well as Carr-Madan FFT aproach.
# Option price is calculated in both way and used to calculate implied volitility.
# 2 maturities were choosed to be ploted.
#



import pandas as pd
from math import exp, sqrt, log, pi
from scipy import stats
from scipy.optimize import fsolve
import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl
mpl.rcParams['font.family'] = 'serif'
from scipy.integrate import quad
np.set_printoptions(suppress=True,formatter={'all': lambda x: '%5.3f' % x})
import scipy.optimize as sop
from numpy.fft import fft

def M76_call_value_INT(S0, div, K, T, r, sigma, Lamb, mu, delta):
    ''' Valuation of European call option in BSM model via Lewis (2001)
    --> Fourier-based approach (integral).
    Parameters
    ==========
    S0: float
    initial stock/index level
    www.it-ebooks.info
    116 DERIVATIVES ANALYTICS WITH PYTHON
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
    int_value = quad(lambda u: M76_integral_function(u, S0, K, T, r-div ,sigma, Lamb, mu, delta), 0, 50, limit=250)[0]
    call_value = S0*np.exp(-div * T) - np.exp(-r * T) * sqrt(S0 * K) / pi * int_value
    return call_value
def M76_integral_function(u, S0, K, T, r, sigma, lamb, mu, delta):
    ''' Valuation of European call option in BSM model via Lewis (2001)
    --> Fourier-based approach: integral function. '''
    cf_value = M76_characteristic_function(u - 1j * 0.5, T, r-div, sigma, lamb, mu, delta)
    int_value = 1 / (u ** 2 + 0.25) * (np.exp(1j * u * np.log(S0 / K)) * cf_value).real
    return int_value
def M76_characteristic_function(u, T, r, sigma, lamb, mu, delta):
    ''' Valuation of European call option in M76 model via
    Lewis (2001) Fourier-based approach: characteristic function.
    Parameter definitions see function M76_value_call_INT. '''
    omega = r - 0.5 * sigma ** 2 - lamb * (np.exp(mu + 0.5 * delta ** 2) - 1)
    value = np.exp((1j * u * omega - 0.5 * u ** 2 * sigma ** 2 +lamb * (np.exp(1j * u * mu - u ** 2 * delta ** 2 * 0.5) - 1)) * T)
    return value
def M76_characteristic_function_FFT(u, x0, T, r, sigma, lamb, mu, delta):
    ''' Valuation of European call option in M76 model via
    Lewis (2001) Fourier-based approach: characteristic function.
    Parameter definitions see function M76_value_call_FFT. '''
    omega = x0 / T + r - 0.5 * sigma ** 2 - lamb * (np.exp(mu + 0.5 * delta ** 2) - 1)
    value = np.exp((1j * u * omega - 0.5 * u ** 2 * sigma ** 2 +lamb * (np.exp(1j * u * mu - u ** 2 * delta ** 2 * 0.5) - 1)) * T)
    return value
def M76_value_call_FFT(S0,div, K, T, r, sigma, lamb, mu, delta):
    ''' Valuation of European call option in M76 model via
    Lewis (1999) Fourier-based approach.
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
    lamb: float
    jump intensity
    mu: float
    expected jump size
    delta: float
    standard deviation of jump
    Returns
    =======
    call_value: float
    European call option present value
    '''
    S0=S0*exp(-div*T)
    k = log(K / S0)
    x0 = log(S0 / S0)
    g = 2 # factor to increase accuracy
    N = g * 4096
    eps = (g * 150.) ** -1
    eta = 2 * pi / (N * eps)
    b = 0.5 * N * eps - k
    u = np.arange(1, N + 1, 1)
    vo = eta * (u - 1)
    # Modificatons to Ensure Integrability
    if S0 >= 0.95 * K: # ITM case
        alpha = 1.5
        v = vo - (alpha + 1) * 1j
        mod_char_fun = exp(-(r) * T) * M76_characteristic_function_FFT(v, x0, T, r, sigma, lamb, mu, delta)/ (alpha ** 2 + alpha - vo ** 2 + 1j * (2 * alpha + 1) * vo)
    else: # OTM case
        alpha = 1.1
        v = (vo - 1j * alpha) - 1j
        mod_char_fun_1 = exp(-(r) * T) * (1 / (1 + 1j * (vo - 1j * alpha))- exp((r) * T) / (1j * (vo - 1j * alpha))- M76_characteristic_function_FFT(v, x0, T, r, sigma, lamb, mu, delta)/ ((vo - 1j * alpha) ** 2 - 1j * (vo - 1j * alpha)))
        v = (vo + 1j * alpha) - 1j
        mod_char_fun_2 = exp(-(r) * T) * (1 / (1 + 1j * (vo + 1j * alpha))- exp((r) * T) / (1j * (vo + 1j * alpha))- M76_characteristic_function_FFT(v, x0, T, r, sigma, lamb, mu, delta)/ ((vo + 1j * alpha) ** 2 - 1j * (vo + 1j * alpha)))
    # Numerical FFT Routine
    delt = np.zeros(N, dtype=np.float)
    delt[0] = 1
    j = np.arange(1, N + 1, 1)
    SimpsonW = (3 + (-1) ** j - delt) / 3
    if S0 >= 0.95 * K:
        fft_func = np.exp(1j * b * vo) * mod_char_fun * eta * SimpsonW
        payoff = (fft(fft_func)).real
        call_value_m = np.exp(-alpha * k) / pi * payoff
    else:
        fft_func = (np.exp(1j * b * vo)* (mod_char_fun_1 - mod_char_fun_2)* 0.5 * eta * SimpsonW)
        payoff = (fft(fft_func)).real
        call_value_m = payoff / (np.sinh(alpha * k) * pi)
    pos = int((k + b) / eps)
    call_value = call_value_m[pos]
    return call_value * S0
#
# Error Function
#
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
def M76_error_function_FFT(p0):
    ''' Error Function for parameter calibration in M76 Model via
    Carr-Madan (1999) FFT approach.
    Parameters
    ==========
    sigma: float
    volatility factor in diffusion term
    lamb: float
    jump intensity
    mu: float
    expected jump size
    delta: float
    standard deviation of jump
    Returns
    =======
    RMSE: float
    root mean squared error
    '''
    global i, min_RMSE
    sigma, lamb, mu, delta = p0
    if sigma < 0.0 or delta < 0.0 or lamb < 0.0:
        return 500.0
    se = []
    for row, option in options.iterrows():
        T = (option['expiration'] - option['date']).days / 365.
        model_value = M76_value_call_FFT(S0,div, option['strike'], T,r, sigma, lamb, mu, delta)
        se.append((model_value - option['price']) ** 2)
    RMSE = sqrt(sum(se) / len(se))
    min_RMSE = min(min_RMSE, RMSE)
    if i % 50 == 0:
        print '%4d |' % i, np.array(p0), '| %7.3f | %7.3f' % (RMSE, min_RMSE)
    i += 1
    return RMSE
def M76_error_function_LEWIS(p0):
    ''' Error Function for parameter calibration in M76 Model via
    Carr-Madan (1999) FFT approach.
    Parameters
    ==========
    sigma: float
    volatility factor in diffusion term
    lamb: float
    jump intensity
    mu: float
    expected jump size
    delta: float
    standard deviation of jump
    Returns
    =======
    RMSE: float
    root mean squared error
    '''
    global i, min_RMSE
    sigma, lamb, mu, delta = p0
    if sigma < 0.0 or delta < 0.0 or lamb < 0.0:
        return 500.0
    se = []
    for row, option in options.iterrows():
        T = (option['expiration'] - option['date']).days / 365.
        model_value = M76_call_value_INT(S0, div, option['strike'], T,r, sigma, lamb, mu, delta)
        se.append((model_value - option['price']) ** 2)
    RMSE = sqrt(sum(se) / len(se))
    min_RMSE = min(min_RMSE, RMSE)
    if i % 50 == 0:
        print '%4d |' % i, np.array(p0), '| %7.3f | %7.3f' % (RMSE, min_RMSE)
    i += 1
    return RMSE
def generate_plot(optFFT, optLEWIS, options):
    #
    # Calculating Model Prices
    #
    sigma, lamb, mu, delta = optFFT
    options['ModelFFT'] = 0.0
    options['ModelLEWIS'] = 0.0
    options['ImpVol'] = 0.0    
    options['ImpVolFFT'] = 0.0
    options['ImpVolLEWIS'] = 0.0
    for row, option in options.iterrows():
        T = (option['expiration'] - option['date']).days / 365.
        options.loc[row, 'ImpVol'] = imp_vol(r,div,S0,option['strike'],T,option['price'])
        options.loc[row, 'ModelFFT'] = M76_value_call_FFT(S0,div, option['strike'], T, r, sigma, lamb, mu, delta)
    sigma, lamb, mu, delta = optLEWIS
    for row, option in options.iterrows():
        T = (option['expiration'] - option['date']).days / 365.
        options.loc[row, 'ModelLEWIS'] = M76_call_value_INT(S0,div, option['strike'], T, r, sigma, lamb, mu, delta)
        options.loc[row, 'ImpVolLEWIS'] = imp_vol(r,div,S0,option['strike'],T,option['ModelLEWIS'])
    for row, option in options.iterrows():
        T = (option['expiration'] - option['date']).days / 365.
        options.loc[row, 'ImpVolFFT'] = imp_vol(r,div,S0,option['strike'],T,option['ModelFFT'])
        options.loc[row, 'ImpVolLEWIS'] = imp_vol(r,div,S0,option['strike'],T,option['ModelLEWIS'])

    #
    # Plotting
    #
    mats = sorted(set(options['expiration']))
    options = options.set_index('strike')
    for i, mat in enumerate(mats):
        options[options['expiration'] == mat][['ImpVol', 'ImpVolFFT']].plot(style=['b-', 'ro'], title='%sFFT' % str(mat)[:10])
        options[options['expiration'] == mat][['price', 'ModelFFT']].plot(style=['b-', 'ro'], title='%sFFT' % str(mat)[:10])        
        plt.ylabel('option value')
        #plt.savefig('../images/08_m76/M76_calibration_FFT_%s.pdf' % i)
        if i ==1:
            break
    for i, mat in enumerate(mats):
        options[options['expiration'] == mat][['ImpVol', 'ImpVolLEWIS']].plot(style=['b-', 'ro'], title='%sLEWIS' % str(mat)[:10])
        options[options['expiration'] == mat][['price', 'ModelLEWIS']].plot(style=['b-', 'ro'], title='%sLEWIS' % str(mat)[:10])                
        plt.ylabel('option value')
        #plt.savefig('../images/08_m76/M76_calibration_LEWIS_%s.pdf' % i)
        if i ==1:
            break
#
# determine parameters
#


r=0.001
div=0.02
S0=1972.18


#
# import data
#

xl = pd.ExcelFile("spxcalls20150831.xlsx")
options = xl.parse('WRDS')
#
# calibration by Carr-Madan
#
i = 0 # counter initialization
min_RMSE = 100 # minimal RMSE initialization
p0FFT = sop.brute(M76_error_function_FFT, ((0.075, 0.201, 0.05),(0.10, 0.401, 0.2), (-0.5, 0.01, 0.2),(0.10, 0.301, 0.2)), finish=None)
# Finding global minimization
optFFT = sop.fmin(M76_error_function_FFT, p0FFT,maxiter=1500, maxfun=1750,xtol=0.000001, ftol=0.000001)
#local minimazation


#
# calibration by Lewis
#

i=0
min_RMSE=100
p0LEWIS = sop.brute(M76_error_function_LEWIS, ((0.075, 0.201, 0.05),(0.10, 0.401, 0.2), (-0.5, 0.01, 0.2),(0.10, 0.301, 0.2)), finish=None)
# Finding global minimization
optLEWIS = sop.fmin(M76_error_function_LEWIS, p0LEWIS,maxiter=1500, maxfun=1750,xtol=0.000001, ftol=0.000001)
#local minimazation
generate_plot(optFFT, optFFT, options)
