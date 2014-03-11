from __future__ import division
import  matplotlib.pyplot as plt
import numpy as np
import coreFunctions as cc
import scipy as sc
from scipy.optimize import curve_fit
from scipy import stats
import sys



def campbell(fig, RADn, RQDn, color, xmin, xmax, ymin,
             ymax, ct, numSims, xdata, ydata, dims, fs):
    
    r, c, i = dims
    ax = fig.add_subplot(r, c, i)
    
    plt.scatter(RADn, RQDn, s = 30, c = color, alpha = 0.4, edgecolor = 'none')        
    
    x = [xmin, xmax]
    y = list(x)
    plt.plot(x, y, 'k', lw=1)
    
    if ct == numSims-1:
        
        r2 = cc.one2one_rsquare(np.array(xdata), np.array(ydata))
        plt.scatter([0],[0], s=0, label=r'R$^2$ = '+str(round(r2,2)))
        plt.xlim(xmin, xmax)
        plt.ylim(ymin, ymax)
        
        plt.title(r'R$^2$ = '+str(round(r2,2)), fontsize=fs-1)
        
        plt.xlabel("log(rel. abundance)", fontsize=fs)
        plt.ylabel("log(rel. cell quota)", fontsize=fs)
        plt.tick_params(axis='both', which='major', labelsize=fs-1)
        
    return fig
        
        
        
def MCQvAb(fig, RADn, AvgQs, color, xmin, xmax, ymin, 
           ymax, xdata, ydata, dims, ct, numSims, fs):

    r, c, i = dims
    ax = fig.add_subplot(r, c, i) 
    plt.scatter(RADn, AvgQs, s = 30, c = color, alpha = 0.4, edgecolor = 'none')        
    
    if ct == numSims-1:
        #plt.title('slope= '+str(round(slope,2))+', p='+str(round(pval,4)), fontsize=fs-1)
            
        plt.xlim(xmin, xmax)
        plt.ylim(ymin, ymax)
        plt.tick_params(axis='both', which='major', labelsize=fs-1)
        plt.xlabel("log(abundance)", fontsize=fs)
        plt.ylabel("mean cell quota", fontsize=fs)        

    return fig



def normRADs(fig, RADn, ranks, color, fig4_xdata, fig4_ydata, dims, ct, numSims, fs): 

    r, c, i = dims
    ax = fig.add_subplot(r, c, i)
    RADn.sort()
    RADn.reverse()
    
    for i, rank in enumerate(ranks):    
        plt.scatter(rank/len(RADn), RADn[i], s=20, c = color, alpha = 0.4, edgecolor = 'none')
        
    if ct == numSims-1:
        slope, intercept, rval, pval, stderr = sc.stats.linregress(fig4_xdata, fig4_ydata)
        X = list(fig4_xdata)
        Y = list(fig4_ydata)
        z = np.polyfit(X,Y,1)
        #print 'r-squared and slope for fig4:',rval**2, slope
        p = np.poly1d(z)
        xp = np.linspace(min(X), max(X), 1000)
        plt.plot(xp,p(xp),'-',c='k',lw=1)
        
        plt.title(r'R$^2$ = '+str(round(rval**2,2))+'\nm = '+str(round(slope,2))+', p ='+str(round(pval,4)), fontsize=fs-1)
        
        plt.xlim(0,1.1)    
        plt.tick_params(axis='both', which='major', labelsize=fs-1)
        plt.xlabel("Rank/Richness", fontsize=fs)
        plt.ylabel("log(abundance)", fontsize=fs)                        

    return fig
    
    
    
def AbvTau(fig, Taus, Ns, colors, dims, ct, numSims, fs):
    
    Taus = list(np.log(Taus))
    
    r, c, j = dims
    ax = fig.add_subplot(r, c, j) 
    
    slope, intercept, r_val, p_val, stderr = sc.stats.linregress(Taus, Ns)
    X = list(Taus)
    Y = list(Ns)
    z = np.polyfit(X,Y,2)
    p = np.poly1d(z)
    xp = np.linspace(min(X), max(X), 1000)
    plt.plot(xp,p(xp),'-',c='k',lw=1)
    
    colors = filter(lambda x: x != 'w', colors)
    for i, color in enumerate(colors): 
        plt.scatter(Taus[i], Ns[i], s=20, c = color, alpha = 0.4, edgecolor = 'none')
    
    plt.title(r'R$^2$ = '+str(round(r_val**2,2))+'\nm = '+str(round(slope,2))+', p = '+str(round(p_val,2)), fontsize=fs-1)
    
    plt.xlim(min(Taus),max(Taus))
    plt.tick_params(axis='both', which='major', labelsize=fs-1)
    plt.xlabel("Residence time", fontsize=fs)
    plt.ylabel("Total\nabundance", fontsize=fs)       
    
    return fig
    
    
def TOvTau(fig, Taus, TOs, colors, dims, ct, numSims, fs):
    
    Taus = list(np.log(Taus))
    
    r, c, j = dims
    ax = fig.add_subplot(r, c, j) 
    
    slope, intercept, r_val, p_val, stderr = sc.stats.linregress(Taus, TOs)
    X = list(Taus)
    Y = list(TOs)
    z = np.polyfit(X,Y,2)
    p = np.poly1d(z)
    xp = np.linspace(min(X), max(X), 1000)
    plt.plot(xp,p(xp),'-',c='k',lw=1)
    
    colors = filter(lambda x: x != 'w', colors)
    for i, color in enumerate(colors): 
        plt.scatter(Taus[i], TOs[i], s=20, c = color, alpha = 0.4, edgecolor = 'none')
    
    plt.title(r'R$^2$ = '+str(round(r_val**2,2))+'\nm = '+str(round(slope,2))+', p = '+str(round(p_val,2)), fontsize=fs-1)
    
    for i, color in enumerate(colors): 
        plt.scatter(Taus[i], TOs[i], s=20, c = color, alpha = 0.4, edgecolor = 'none')
    
    plt.xlim(min(Taus),max(Taus))
    plt.tick_params(axis='both', which='major', labelsize=fs-1)
    plt.xlabel("Residence time", fontsize=fs)
    plt.ylabel("Biomass\nturnover", fontsize=fs)       
    
    return fig
    
    
    
def CTvTau(fig, Taus, CTs, colors, dims, ct, numSims, fs):
    
    Taus = list(np.log(Taus))
    
    r, c, j = dims
    ax = fig.add_subplot(r, c, j) 
    
    slope, intercept, r_val, p_val, stderr = sc.stats.linregress(Taus, CTs)
    X = list(Taus)
    Y = list(CTs)
    z = np.polyfit(X,Y,2)
    p = np.poly1d(z)
    xp = np.linspace(min(X), max(X), 1000)
    plt.plot(xp,p(xp),'-',c='k',lw=1)
    
    colors = filter(lambda x: x != 'w', colors)
    for i, color in enumerate(colors): 
        plt.scatter(Taus[i], CTs[i], s=20, c = color, alpha = 0.4, edgecolor = 'none')
    
    plt.title(r'R$^2$ = '+str(round(r_val**2,2))+'\nm = '+str(round(slope,2))+', p = '+str(round(p_val,2)), fontsize=fs-1)
    
    colors = filter(lambda x: x != 'w', colors)
    for i, color in enumerate(colors): 
        plt.scatter(Taus[i], CTs[i], s=20, c = color, alpha = 0.4, edgecolor = 'none')
    
    plt.tick_params(axis='both', which='major', labelsize=fs-1)
    plt.xlim(min(Taus),max(Taus))
    plt.xlabel("Residence time", fontsize=fs)
    plt.ylabel("Compositional\nturnover", fontsize=fs)
    
    return fig
 
       
def SvTau(fig, Taus, Ss, colors, dims, ct, numSims, fs):
    
    Taus = list(np.log(Taus))
    
    r, c, j = dims
    ax = fig.add_subplot(r, c, j) 
    
    slope, intercept, r_val, p_val, stderr = sc.stats.linregress(Taus, Ss)
    X = list(Taus)
    Y = list(Ss)
    z = np.polyfit(X,Y,2)
    p = np.poly1d(z)
    xp = np.linspace(min(X), max(X), 1000)
    plt.plot(xp,p(xp),'-',c='k',lw=1)
    
    colors = filter(lambda x: x != 'w', colors)
    for i, color in enumerate(colors): 
        plt.scatter(Taus[i], Ss[i], s=20, c = color, alpha = 0.4, edgecolor = 'none')
    
    plt.title(r'R$^2$ = '+str(round(r_val**2,2))+'\nm = '+str(round(slope,2))+', p = '+str(round(p_val,2)), fontsize=fs-1)
    
    plt.xlim(min(Taus),max(Taus))
    plt.tick_params(axis='both', which='major', labelsize=fs-1)
    plt.xlabel("Residence time", fontsize=fs)
    plt.ylabel("no. Taxa", fontsize=fs)       
    
    return fig
