# -*- coding: utf-8 -*-
"""
Created on Fri Jan 15 21:34:22 2016

@author: Jianan Wang
"""
#
# This program read APPL data, calculate realized volitility and,
# calculate call option price using Black Shole's Model.
# Finally plot option price to equalty price.
#
import math
import numpy as np
import pandas.io.data as web
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy.integrate import quad
mpl.rcParams['font.family'] = 'serif'



def dN(x):
    ''' Probability density function of standard normal random variable x. '''
    return math.exp(-0.5 * x ** 2) / math.sqrt(2 * math.pi)


def N(d):
    ''' Cumulative density function of standard normal random variable x. '''
    return quad(lambda x: dN(x), -20, d, limit=50)[0]


def d1f(St, K, t, T, r, sigma):
    ''' Black-Scholes-Merton d1 function.
        Parameters see e.g. BSM_call_value function. '''
    d1 = (math.log(St / K) + (r + 0.5 * sigma ** 2)
            * (T - t)) / (sigma * math.sqrt(T - t))
    return d1

#
# Valuation Functions
#


def BSM_call_value(St, K, t, T, r, sigma):
    ''' Calculates Black-Scholes-Merton European call option value.

    Parameters
    ==========
    St : float
        stock/index level at time t
    K : float
        strike price
    t : float
        valuation date
    T : float
        date of maturity/time-to-maturity if t = 0; T > t
    r : float
        constant, risk-less short rate
    sigma : float
        volatility

    Returns
    =======
    call_value : float
        European call present value at t
    '''
    d1 = d1f(St, K, t, T, r, sigma)
    d2 = d1 - sigma * math.sqrt(T - t)
    call_value = St * N(d1) - math.exp(-r * (T - t)) * K * N(d2)
    return call_value

 
def read_aapl_data():
    ''' Reads historical AAPL data from Yahoo! Finance, calculates log returns, 
    realized variance and volatility.''' 
    AAPL = web.DataReader('AAPL', data_source='yahoo',
                    start='01-01-2015', end='31-12-2015')
    AAPL.rename(columns={'Adj Close' : 'index'}, inplace=True)
    AAPL['returns'] = np.log(AAPL['index'] / AAPL['index'].shift(1))
    AAPL['std']=np.std(AAPL['returns'])    
    AAPL['rea_var'] = 252 * np.cumsum(AAPL['returns'] ** 2) / np.arange(len(AAPL))
    AAPL['rea_vol'] = np.sqrt(AAPL['rea_var'])
    AAPL = AAPL.dropna()
    return AAPL
    
    
def print_statistics(data) :
    std=np.sqrt(251)*np.std(data['returns'])
    print" RETURN AAPL STATISTICS"
    print"---------------------------------------"
    print" Annual Std %9.6f" %  std
    print"---------------------------------------"
    return
appl=read_aapl_data()
print_statistics(appl)
re=appl.returns
sd=appl.std
rv=appl.rea_vol
K=100
vol=0.35
r=0.0025
T=0.15
t=np.linspace(1,251,251)
# Sample Data Generation
S = np.linspace(80, 120, 200) 							 # vector of index level values
h = np.maximum(S - K, 0)  								# inner value of option
C = [BSM_call_value(S0, K, 0, T, r, vol) for S0 in S]  # calculate call option values
plt.figure()
plt.plot(t,re,'b', lw=2.5, label=' daily return')
plt.grid(True)
plt.legend(loc=0)
plt.xlabel('time trend')
plt.ylabel('daily return')

plt.figure()
plt.plot(t,rv,'y', lw=2.5, label=' realized volatility')
plt.grid(True)
plt.legend(loc=0)
plt.xlabel('time trend')
plt.ylabel('realized volatility')



# Graphical Output
plt.figure()
plt.plot(S, h, 'b', lw=2.5, label='pay off diagram') 	# plot inner value at maturity
plt.plot(S, C, 'r', lw=2.5, label='BSM price')			# plot option present value
plt.grid(True)
plt.legend(loc=0)
plt.xlabel('spot price')
plt.ylabel('pay-off/ price')


