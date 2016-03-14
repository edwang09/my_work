# -*- coding: utf-8 -*-
"""
Created on Fri Jan 15 22:27:43 2016

@author: mahone
"""
#
# This program define an object called eropean_option, that has attributes S0 (initial underlying level),
# K (option strike price), t (pricing date), M (maturity date), r (constant risk-free short rate),
# sigma (volatility), d (continuous dividend rate), CP (call or put), and an optional C (market price of the option), 
# and methods value (return present value of the option), imp_vol (implied volatility if market price is given, 
# delta (option delta), gamma(option gamma), vega (option vega), theta (option theta), rho (option rho)
#


from math import log, sqrt, exp
from scipy import stats
from scipy.optimize import fsolve

class eropean_option(object):
    ''' Class for European call options in BSM Model.
    
    Attributes
    ==========
    S0 : float
        initial underlying level
    K : float
        option strike price
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
    
    def __init__(self, S0, K, t, M, r, sigma, d, CP, C):
        self.S0 = float(S0)
        self.K = K
        self.t = t
        self.M = M
        self.r = r
        self.sigma = sigma
        self.d=d
        self.CP=CP
        self.C=C

    def update_ttm(self):
        ''' Updates time-to-maturity self.T. '''
        if self.t > self.M:
            raise ValueError("Pricing date later than maturity.")
        self.T = (self.M - self.t).days / 365.

    def d1(self):
        ''' Helper function. '''
        d1 = ((log(self.S0 / self.K)
            + (self.r-self.d + 0.5 * self.sigma ** 2) * self.T)
            / (self.sigma * sqrt(self.T)))
        return d1
        
    def value(self):
        ''' Return option value. '''
        self.update_ttm()
        d1 = self.d1()
        d2 = self.d1()-self.sigma*sqrt(self.T)
            
        if self.CP=='c':
            value = (self.S0 * exp(-self.d * (self.T))* stats.norm.cdf(d1, 0.0, 1.0)
                    - self.K * exp(-self.r * self.T) * stats.norm.cdf(d2, 0.0, 1.0))
        else:
             value = -(self.S0 * exp(-self.d * (self.T))* stats.norm.cdf(-d1, 0.0, 1.0)
                    + self.K * exp(-self.r * self.T) * stats.norm.cdf(-d2, 0.0, 1.0))
        return value
        
 

    def imp_vol(self, sigma_est=0.2):
        ''' Return implied volatility given option price. '''
        option = eropean_option(self.S0, self.K, self.t, self.M,
                             self.r, sigma_est)
        option.update_ttm()
        def difference(sigma):
            option.sigma = sigma
            return option.value() - self.C
        iv = fsolve(difference, sigma_est)[0]
        return iv
    def delta(self):
        self.update_ttm()
        d1=self.d1()
        if self.CP=='c':
            delta=exp(-self.d*self.T)*stats.norm.pdf(d1, 0.0, 1.0)
        else:
            delta=-exp(-self.d*self.T)*stats.norm.pdf(-d1, 0.0, 1.0)
        return delta
        
    def gamma(self):
        self.update_ttm()
        d1=self.d1()
        gamma=exp(-self.d*self.T)*stats.norm.pdf(d1, 0.0, 1.0)/(self.S0*self.sigma*sqrt(self.T))
        return gamma
    
    def vega(self):
        ''' Return Vega of option. '''
        self.update_ttm()
        d1 = self.d1()
        vega = self.S0 * exp(-self.d * (self.T))* stats.norm.pdf(d1, 0.0, 1.0) * sqrt(self.T)
        return vega
        
    def theta(self):
        self.update_ttm()
        d1=self.d1()
        d2=self.d1()-self.sigma*sqrt(self.T)
        temp=-exp(-self.d*self.T)*self.S0*self.sigma*stats.norm.pdf(d1, 0.0, 1.0)/(2*sqrt(self.T))
        if self.CP=='c':
            theta=temp+self.d*(self.S0 * exp(-self.d * (self.T))* stats.norm.cdf(d1, 0.0, 1.0)
                    - self.r*self.K * exp(-self.r * self.T) * stats.norm.cdf(d2, 0.0, 1.0))
        else:
            theta=temp-self.d*(self.S0 * exp(-self.d * (self.T))* stats.norm.cdf(-d1, 0.0, 1.0)
                    + self.r*self.K * exp(-self.r * self.T) * stats.norm.cdf(-d2, 0.0, 1.0))
        return theta
    
    def rho(self):
        self.update_ttm()
        d2=self.d1()-self.sigma*sqrt(self.T)
        if self.CP=='c':
            rho=self.K*self.T*exp(-self.r*self.T)*stats.norm.cdf(d2, 0.0, 1.0)
        else:
            rho=-self.K*self.T*exp(-self.r*self.T)*stats.norm.cdf(-d2, 0.0, 1.0)
        return rho