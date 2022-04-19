#!/usr/bin/env python
# coding: utf-8

# In[27]:


import pandas as pd
import numpy as np
from scipy.signal import find_peaks, peak_widths
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt


# In[28]:


e=1.60218*10**(-19)
m_e=9.10938356*10**(-31)
h=6.62607015*10**(-31)


# In[29]:


def gaussian_2(x, a, b, c, d, f, g):
    return a*np.exp(-np.power(x - b, 2)/(2*np.power(c, 2))) + d*np.exp(-np.power(x - f, 2)/(2*np.power(g, 2)))
def gaussian(x, a, b, c, d):
    return a*np.exp(-np.power(x - b, 2)/(2*np.power(c, 2)))+d


# In[30]:


sheet_names = ['trans 4']


# In[31]:


referance_frame=pd.DataFrame()

referance_frame_raw = pd.read_excel ('Transverse.xlsx', sheet_name="trans 0")
alpha=referance_frame_raw["alpha"]
pct=referance_frame_raw["Pct"]
plt.plot(alpha, pct)
plt.title('Raw Alpha v Intesity')
plt.xlabel('Alpha(degrees)')
plt.ylabel('Intesity')
plt.show()
peaks, _ = find_peaks(pct, height=20)
range_window = peak_widths(pct, peaks, rel_height=0.5)
plt.plot(alpha[peaks], pct[peaks], "x")
referance_frame["alphas"]=alpha[peaks]
referance_frame["pct"]=pct[peaks]
referance_frame["range"]=range_window[0]
referance_frame=referance_frame.reset_index()
referance_frame["hieght"]=0.0
referance_frame["position"]=0.0
referance_frame["SD"]=0.0
referance_frame["BG"]=0.0
referance_frame["peak1"]=0.0
referance_frame["peak2"]=0.0
for i in range(20):
    integer_location = np.where(referance_frame_raw.index == referance_frame["index"][i])[0][0]
    start = max(0, integer_location - int((referance_frame["range"][i]/2)*2.8))
    end = max(1, integer_location + int((referance_frame["range"][i]/2)*4.0))
    df_window = referance_frame_raw.iloc[start:end]
    peaks, _ = find_peaks(df_window["Pct"], height=20)
    pars, cov = curve_fit(f=gaussian,
                          xdata=df_window["alpha"],
                          ydata=df_window["Pct"],
                          p0=[referance_frame["pct"][i],
                              referance_frame["alphas"][i],
                              0.01,
                              3],
                          bounds=(-np.inf, np.inf))
    a, b, c, d= pars
    referance_frame["hieght"][i]=a
    referance_frame["position"][i]=b
    referance_frame["SD"][i]=abs(c)
    referance_frame["BG"][i]=d
    plt.figure()
    plt.title('Raw Alpha v Intesity')
    plt.plot(df_window["alpha"],df_window["Pct"], "x")
    plt.plot(df_window["alpha"],gaussian(df_window["alpha"],a,b,c,d))
    plt.xlabel('Alpha(degrees)')
    plt.ylabel('Intesity')
    plt.show()
referance_frame=referance_frame.truncate(after=19)
referance_frame = referance_frame.drop(labels=[9,10], axis=0)
referance_frame=referance_frame.reset_index(drop=True)


# In[32]:


results=pd.DataFrame(columns=['i','peak1', 'peak2'])
B_values = pd.DataFrame([[1.0,0.044],
                         [2.0,0.0822],
                         [3.0,0.1221],
                         [4.0,0.1652],
                         [5.0,0.1975],
                         [6.0,0.2280],
                         [7.0,0.253],
                         [8.0,0.269],
                         [9.0,0.285]],
                        columns=['I','B'])


# In[38]:


for sheet_name in sheet_names:
    df = pd.read_excel ('Transverse.xlsx', sheet_name=sheet_name)
    for i in range(18):
        background=referance_frame["BG"][i]
        integer_location = np.where(df.index == referance_frame["index"][i])[0][0]
        start = max(0, integer_location - int((referance_frame["range"][i]/2)*2.2))
        end = max(1, integer_location + int((referance_frame["range"][i]/2)*4.2))
        df_window = df.iloc[start:end]
        peaks, _ = find_peaks(df_window["Pct"], height=20)
        pars, cov = curve_fit(f=gaussian_2, 
                              xdata=df_window["alpha"], 
                              ydata=df_window["Pct"],
                              p0=[referance_frame["hieght"][i],
                                 referance_frame["position"][i], 
                                 referance_frame["SD"][i], 
                                 referance_frame["hieght"][i], 
                                 referance_frame["position"][i], 
                                 referance_frame["SD"][i]],
                              bounds=((10, #a
                                       -np.inf, #b
                                       referance_frame["SD"][i]-0.002, #c
                                       -np.inf, #d
                                       -np.inf,#f
                                       referance_frame["SD"][i]-0.002),#g
                                      (np.inf,#a
                                       np.inf,#b
                                       referance_frame["SD"][i],#c
                                       np.inf,#d
                                       np.inf,#f
                                       referance_frame["SD"][i])))#g
        a, b, c, d, f, g= pars
        plt.figure()
        plt.plot(df_window["alpha"],df_window["Pct"], "x", label="Experimental")
        plt.plot(df_window["alpha"],gaussian_2(df_window["alpha"],a,b,c,d,f,g), label="Full Fit")
        plt.plot(df_window["alpha"],gaussian(df_window["alpha"],a,b,c,0), color="red", label="Gaussian 1")
        plt.plot(df_window["alpha"],gaussian(df_window["alpha"],d,f,g,0), color="green", label="Gaussian 2")
        plt.title('Double Guassian Fitting of Overlapping Peaks ' + sheet_name)
        plt.legend(loc="upper right")
        plt.xlabel('Alpha(degrees)')
        plt.ylabel('Intesity')
        plt.show()
        B=(B_values["B"][float(sheet_name[-1])-1])
        df2 = pd.DataFrame([[i,B,referance_frame["position"][i],b,f]],columns=['i','B','peakRef','peak1', 'peak2'] )
        results=results.append(df2)
        #print(df2)
print("full run") 
print(results)


# In[34]:


results["peakRef"]=(results['peak1']+results['peak2'])/2 
results["betaRef"]=abs(np.arcsin((1/1.457)*np.sin((results['peakRef'])*((np.pi)/180))))
results["beta1"]=abs(np.arcsin((1.457)*np.sin((results['peak1'])*((np.pi)/180))))
results["beta2"]=abs(np.arcsin((1.457)*np.sin((results['peak2'])*((np.pi)/180))))
results["delta_lambda/lambda1"]=(np.cos(results["beta1"])/np.cos(results["betaRef"]))-1.
results["delta_lambda/lambda2"]=(np.cos(results["beta2"])/np.cos(results["betaRef"]))-1.
results["E"]=(h/(4*np.pi))*(e/m_e)*results["B"]
results["deltaE1"]=abs(results["E"]*results["delta_lambda/lambda1"])
results["deltaE2"]=abs(results["E"]*results["delta_lambda/lambda2"])
print(results)


# In[35]:


plt.figure()
plt.plot(results["B"],abs(results["deltaE1"]), 'xb')
plt.plot(results["B"],abs(results["deltaE2"]), 'xc')
plt.show()
print(results["deltaE2"].mean(),results["deltaE1"].mean())
print((results["deltaE2"].mean()+results["deltaE1"].mean())/2)


# In[ ]:





# In[ ]:




