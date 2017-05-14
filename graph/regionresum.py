# -*- coding: utf-8 -*-
"""
Created on Mon May  8 14:10:29 2017

@author: claudio
"""

def regionplot():

    import matplotlib.pyplot as plt;
    import numpy as np;

    x=np.arange(-0.1,100.,0.1);   
    y1=0.1*x/x;
    y0=0.0*x/x;
    y2= (np.sqrt(1.+x*x/(125.09*125.09))-np.sqrt(x*x/(125.09*125.09)))*(np.sqrt(1.+x*x/(125.09*125.09))-np.sqrt(x*x/(125.09*125.09)));
    y3=(np.sqrt(1.+x*x/(125.09*125.09))-np.sqrt(x*x/(125.09*125.09)))*(np.sqrt(1.+x*x/(125.09*125.09))-np.sqrt(x*x/(125.09*125.09)))-0.1;
    
    fig,ax=plt.subplots();
    ax.plot(x,y1,x,y2,x,y3,color='black');
    ax.fill_between(x,y0,y3, where=x <=10., alpha=0.7, facecolor='green');
    ax.fill_between(x,y3,y2, where=x>=10., facecolor= 'red');
    ax.fill_between(x,y0,y1, alpha=0.7,facecolor= 'blue',);
    ax.fill_between(x,y3,y2, where=x<=10.,alpha=1.0,hatch='///', edgecolor='green',facecolor='white'); 
    ax.fill_between(x,y3,y2, where=x<=10.,alpha=0.5, hatch=' - ', edgecolor='red',facecolor='white', color='white');    
    
    
    plt.xlabel(r"$p_T [GeV]$", fontsize=14)
    plt.ylabel(r"$\tau=\frac{m_H^2}{s}$ ", fontsize=14) 
    plt.axis([0.0,100.0,0.,1.0])
    
    plt.savefig("/home/claudio/Joint/graph/regres1.pdf", format='pdf')
    plt.show();
    
    y4=1.0*x/x;
    y5=1.0*x/x-0.1;    
    
    fig2,ax2=plt.subplots();
    ax2.plot(x,y1,x,y4,x,y5,color='black');
    ax2.fill_between(x,y0,y5, where=x <=10., alpha=0.7, facecolor='green');
    ax2.fill_between(x,y5,y4, where=x>=10., facecolor= 'red');
    ax2.fill_between(x,y0,y1, alpha=0.7,facecolor= 'blue',);
    ax2.fill_between(x,y5,y4, where=x<=10.,alpha=1.0,hatch='///', edgecolor='green',facecolor='white'); 
    ax2.fill_between(x,y5,y4, where=x<=10.,alpha=0.5, hatch=' \\\ ', edgecolor='red',facecolor='white', color='white');    
    
    
    plt.xlabel(r"$p_T$ [GeV]", fontsize=14)
    plt.ylabel(r"$\tau'=\frac{\tau}{\tau^{-^{\,}}}=\frac{(\sqrt{m_H^2+p_T^2}+p_T)^2}{s}$", fontsize=14) 
    plt.axis([0.0,100.0,0.,1.0])
    
    plt.savefig("/home/claudio/Joint/graph/regres2.pdf", format='pdf')
    plt.show();
    return 0;
    
def regionplot2():
    import matplotlib.pyplot as plt;
    import numpy as np;

    x=np.arange(-0.1,100.,0.1);   
    y1=0.1*x/x;
    y0=0.0*x/x;
    
    
    y4=1.0*x/x;
    y5=1.0*x/x-0.1;    
    
    fig2,ax2=plt.subplots();
    ax2.plot(x,y1,x,y4,x,y5,color='black');
    ax2.fill_between(x,y0,y5, where=x <=10., alpha=0.7, facecolor='green');
    ax2.fill_between(x,y5,y4, where=x>=10., facecolor= 'red');
    ax2.fill_between(x,y0,y1, alpha=0.7,facecolor= 'blue',);
    ax2.fill_between(x,y5,y4, where=x<=10.,alpha=1.0,hatch='///', edgecolor='green',facecolor='white'); 
    ax2.fill_between(x,y5,y4, where=x<=10.,alpha=0.5, hatch=' \\\ ', edgecolor='red',facecolor='white', color='white');    
    
    
    plt.xlabel(r"$p_T$ [GeV]", fontsize=14)
    plt.ylabel(r"$\tau'=\frac{\tau}{\tau^{-^{\,}}}=\frac{(\sqrt{m_H^2+p_T^2}+p_T)^2}{s}$", fontsize=14) 
    plt.axis([0.0,100.0,0.,1.0])
    
    plt.savefig("/home/claudio/Joint/graph/regres2.pdf", format='pdf')
    plt.show();
    return 0;
    
    
regionplot();    
#regionplot2();
    
        