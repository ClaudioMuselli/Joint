# -*- coding: utf-8 -*-
"""
Created on Sun May 14 16:10:25 2017

@author: claudio
"""

def plotresum():
    import matplotlib.pyplot as plt
    import reading as read
    
    name1="PlotResumFinal_13_LL_tiny.dat"
    name2="PlotResumFinal_13_NLL_tiny.dat"
    name3="PlotResumFinal_13_NNLL_tiny.dat"
    
    name4="../NLLsmallpt.dat" 
    name5="../NNLLsmallpt.dat"
    
    ptLL=read.obtain_x(name1,2)
    ptNLL=read.obtain_x(name2,2)
    ptNNLL=read.obtain_x(name3,2)
    
    JointLL=read.obtain_x(name1,3)
    JointNLL=read.obtain_x(name2,3)
    JointNNLL=read.obtain_x(name3,3)
    
    FixNLL=read.obtain_x(name2,5)
    FixNNLL=read.obtain_x(name3,5)
    
    MatchedLL=read.obtain_x(name1,7)
    MatchedNLL=read.obtain_x(name2,7)
    MatchedNNLL=read.obtain_x(name3,7)
    
    ResptLL=read.obtain_x(name1,9)
    ResptNLL=read.obtain_x(name2,9)
    ResptNNLL=read.obtain_x(name3,9)
    
    ptfreeNLL=read.obtain_x(name4,1)    
    Li2freeNLL=read.obtain_x(name4,2)
    ptfreeNNLL=read.obtain_x(name5,1)
    Li2freeNNLL=read.obtain_x(name5,2)
    
    line1LL, line2LL, line3LL =plt.plot(ptLL,JointLL, ptLL,MatchedLL,ptLL,ResptLL)
    plt.setp(line1LL,color='b',linewidth=2.0,ls='-')
    plt.setp(line2LL,color='r',linewidth=2.0,ls='--')
    plt.setp(line3LL,color='black', linewidth=2.0,ls=':')
    plt.legend([line1LL,line2LL,line3LL],[r"LL Small-$p_T$ Consistent Resummation","LL Combined Resummation",r"LL Small-$p_T$ CSS Resummation"], loc='upper right',fontsize=11, frameon=False)
    plt.title(r"LHC13TeV, PDF4LHC15_nnlo_100 pdf set, $\mu_F=\mu_R=m_H$")    
    plt.axis([0,250,-0.1,1])    
    plt.xlabel(r"$p_T$ [GeV]", fontsize=14)
    plt.ylabel(r"$\frac{d\sigma}{dp_T}$ [pb]", fontsize=14)
    plt.savefig("LL.pdf", format='pdf')
    plt.close() 
    
    line1NLL, line2NLL, line3NLL, line4NLL =plt.plot(ptNLL,JointNLL, ptNLL,MatchedNLL,ptNLL,ResptNLL, ptNLL, FixNLL)
    plt.setp(line1NLL,color='b',linewidth=2.0,ls='-')
    plt.setp(line2NLL,color='r',linewidth=2.0,ls='--')
    plt.setp(line3NLL,color='black', linewidth=2.0,ls=':')
    plt.setp(line4NLL,color='orange', linewidth=1.5,ls='-.')
    plt.legend([line1NLL,line2NLL,line3NLL,line4NLL],[r"NLL Small-$p_T$ Consistent Resummation","NLL Combined Resummation",r"NLL Small-$p_T$ CSS Resummation",r"NLL Fixed-$p_T$ threshold Resummation"], loc='upper right',fontsize=11, frameon=False)
    plt.title(r"LHC13TeV, PDF4LHC15_nnlo_100 pdf set, $\mu_F=\mu_R=m_H$")    
    plt.annotate(r"      $F_{match}=\frac{N^2 \xi_p^3}{1+N^2 \xi_p^3}$",xy=(0,0),xytext=(0,0), horizontalalignment='left')        
    plt.axis([0,250,-0.1,1])
    plt.xlabel(r"$p_T$ [GeV]", fontsize=14)
    plt.ylabel(r"$\frac{d\sigma}{dp_T}$ [pb]", fontsize=14)
    plt.savefig("NLL.pdf", format='pdf')
    plt.close()  
    
    line1NNLL, line2NNLL, line3NNLL, line4NNLL =plt.plot(ptNNLL,JointNNLL, ptNNLL,MatchedNNLL,ptNNLL,ResptNNLL, ptNNLL, FixNNLL)
    plt.setp(line1NNLL,color='b',linewidth=2.0,ls='-')
    plt.setp(line2NNLL,color='r',linewidth=1.0,ls='--')
    plt.setp(line3NNLL,color='black', linewidth=1.0,ls=':')
    plt.setp(line4NNLL,color='orange', linewidth=1.5,ls='-.')
    plt.title(r"LHC13TeV, PDF4LHC15_nnlo_100 pdf set, $\mu_F=\mu_R=m_H$")    
    plt.legend([line1NNLL,line2NNLL,line3NNLL,line4NNLL],[r"NNLL Small-$p_T$ Consistent Resummation","NNLL Combined Resummation",r"NNLL Small-$p_T$ CSS Resummation",r"NNLL Fixed-$p_T$ threshold Resummation"], loc='upper right',fontsize=11, frameon=False)    
    plt.annotate(r"      $F_{match}=\frac{N^2 \xi_p^3}{1+N^2 \xi_p^3}$",xy=(0,0),xytext=(0,0), horizontalalignment='left')    
    plt.axis([0,250,-0.1,1])
    plt.xlabel(r"$p_T$ [GeV]", fontsize=14)
    plt.ylabel(r"$\frac{d\sigma}{dp_T}$ [pb]", fontsize=14)
    plt.savefig("NNLL.pdf", format='pdf')
    plt.close() 
    
    line1pt,line2pt,line3pt =plt.plot(ptLL,ResptLL,ptNLL,ResptNLL,ptNNLL,ResptNNLL)
    plt.setp(line1pt,color='black', linewidth=1.0,ls=':')
    plt.setp(line2pt,color='black', linewidth=1.0,ls='--')
    plt.setp(line3pt,color='black', linewidth=1.0,ls='-')
    plt.title(r"Small-$p_T$ CSS Resummation")
    plt.legend([line1pt,line2pt,line3pt],["LL", "NLL", "NNLL"], loc='upper right',fontsize=11, frameon=False)
    plt.axis([0,250,-0.1,1])
    plt.xlabel(r"$p_T$ [GeV]", fontsize=14)
    plt.ylabel(r"$\frac{d\sigma}{dp_T}$ [pb]", fontsize=14)
    plt.savefig("Smallptseries.pdf", format='pdf')
    plt.close()
    
    line1j,line2j,line3j =plt.plot(ptLL,JointLL,ptNLL,JointNLL,ptNNLL,JointNNLL)
    plt.setp(line1j,color='b', linewidth=1.0,ls=':')
    plt.setp(line2j,color='b', linewidth=1.0,ls='--')
    plt.setp(line3j,color='b', linewidth=1.0,ls='-')
    plt.title(r"Small-$p_T$ Consistent Resummation")
    plt.legend([line1j,line2j,line3j],["LL", "NLL", "NNLL"], loc='upper right',fontsize=11, frameon=False)
    plt.axis([0,250,-0.1,1])
    plt.xlabel(r"$p_T$ [GeV]", fontsize=14)
    plt.ylabel(r"$\frac{d\sigma}{dp_T}$ [pb]", fontsize=14)
    plt.savefig("Jointseries.pdf", format='pdf')
    plt.close()
    
    line1fin,line2fin,line3fin,line4fin,line5fin,line6fin=plt.plot(ptNLL,JointNLL,ptfreeNLL,Li2freeNLL,ptNNLL,JointNNLL,ptfreeNNLL,Li2freeNNLL,ptNLL,ResptNLL,ptNNLL,ResptNNLL)
    plt.setp(line1fin,color='b', linewidth=1.5,ls='--')
    plt.setp(line2fin,color='b', linewidth=1.5,ls=':')
    plt.setp(line3fin,color='r', linewidth=1.5,ls='--')
    plt.setp(line4fin,color='r', linewidth=1.5,ls=':')
    plt.setp(line5fin,color='black', linewidth=1.0, ls='-')
    plt.setp(line6fin,color='black', linewidth=1.0, ls='-.')
    plt.title("Study of the impact of threshold genuine contributions")    
    plt.legend([line1fin,line2fin,line3fin,line4fin,line5fin,line6fin],["NLL Joint", r"NLL Joint - no $Li_2$", "NNLL Joint", r"NNLL Joint - no $Li_2$", "NLL CSS", "NNLL CSS"], loc='lower center',fontsize=11, frameon=False)
    plt.axis([0,40,-0.1,1])
    plt.xlabel(r"$p_T$ [GeV]", fontsize=14)
    plt.ylabel(r"$\frac{d\sigma}{dp_T}$ [pb]", fontsize=14)
    plt.savefig("Li2impact.pdf", format='pdf')
    plt.close()




plotresum()
    