"""
Created on Fri Jan 15 22:27:43 2016

@author: mahone
"""

from math import log, sqrt, exp
from scipy import stats
import matplotlib.pyplot as plt
import numpy as np


class risk_reversal(object):
    ''' Class for European call options in BSM Model.
    
    Attributes
    ==========
    S0 : float
        initial underlying level
    K1 : float
        option1 put strike price
    K2 : float
        option2 call strike price
    t : datetime/Timestamp object
        pricing date
    M : datetime/Timestamp object
        maturity date
    r : float
        constant risk-free short rate
    sigma : float
        volatility factor in diffusion term
    d : float
        constant diviend rate
    CP : charecter(c/p)
        call or put
    C : float
        market price of the option
    Methods
    =======
    value : float
        return present value of call option
    vega : float
        return vega of call option
    imp_vol : float
        return implied volatility given option quote
    '''
    
    def __init__(self, S0, K1, K2, t, M, r, div, sigma):
        self.S0 = float(S0)
        self.K1 = K1
        self.K2 = K2
        self.t = t
        self.M = M
        self.r = r
        self.div=div
        self.sigma = sigma
    def update_ttm(self):
        ''' Updates time-to-maturity self.T. and test for valid strike'''
        if self.t > self.M:
            raise ValueError("Pricing date later than maturity.")
        if self.K1 > self.K2:
            raise ValueError("Put strike larger than call strike.")
        self.T = (self.M - self.t).days / 365.

    def dc1(self):
        ''' Helper function. '''
        self.update_ttm()
        d1 = ((log(self.S0 / self.K2)
            + (self.r-self.div + 0.5 * self.sigma ** 2) * self.T)
            / (self.sigma * sqrt(self.T)))
        return d1
    def dp1(self):
        ''' Helper function. '''
        self.update_ttm()
        d1 = ((log(self.S0 / self.K1)
            + (self.r-self.div + 0.5 * self.sigma ** 2) * self.T)
            / (self.sigma * sqrt(self.T)))
        return d1
        
    def valuec(self):
        ''' Return option value. '''
        self.update_ttm()
        d1 = self.dc1()
        d2 = d1-self.sigma*sqrt(self.T)
        value = (self.S0 * exp(-self.div * (self.T))* stats.norm.cdf(d1, 0.0, 1.0)
                    - self.K2 * exp(-self.r * self.T) * stats.norm.cdf(d2, 0.0, 1.0))
        return value
    def valuep(self):
        ''' Return option value. '''
        self.update_ttm()
        d1 = self.dp1()
        d2 = d1-self.sigma*sqrt(self.T)
        value = -self.S0 * exp(-self.div * (self.T))* stats.norm.cdf(-d1, 0.0, 1.0) + self.K1 * exp(-self.r * self.T) * stats.norm.cdf(-d2, 0.0, 1.0)
        return value
    def value(self):
        value=self.valuec() - self.valuep()
        return value        
    def plot_payoff(self):
        S=np.linspace(self.K1*0.9,self.K2*1.1,1000)
        PO=[]
        P=[]
        for j in range(0,1000):
            self.S0=S[j]
            PO.append( np.maximum(S[j] - self.K2, 0)-np.maximum(self.K1-S[j], 0))
            P.append(self.value())
        plt.figure()
        plt.plot(S, PO, 'b', lw=2.5, label='pay off diagram') 
        plt.plot(S, P, 'r', lw=2.5, label='BSM price')
        plt.grid(True)
        plt.legend(loc=0)
        plt.xlabel('stock price')
        plt.ylabel('price')
        return 0