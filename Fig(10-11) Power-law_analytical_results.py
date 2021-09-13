#Infiltration in an initially-dry soil with power-law porosity decay
#Mohammad Afzal Shadab and Marc Hesse
#Date modified: 08/05/21

from supporting_tools import *

#parameters
simulation_name = f'power_law_analytical'
m = 3 #Cozeny-Karman coefficient for numerator K = K0 (1-phi_i)^m
n = 2 #Corey-Brooks coefficient krw = krw0 * sw^n
s_wr = 0.0 #Residual water saturation
s_gr = 0.0 #Residual gas saturation
phi_0 = 0.5#Porosity at the surface
z0 = 1.0   #point specified for power law porosity variation
p  = 7.63  #index specified for power law porosity variation 
R  = 0.16  #dimensionless rainfall rates

#spatial discretization
zbottom = 0.3
Nz   = 1000
zc   = np.linspace(0,zbottom,Nz) #depth array
phi  = np.transpose([phi_0*(1-zc/z0)**p])#porosity vector: upper layer

#temporal discretization
tmax = 0.5   #time scaling with respect to z0/fc
Nt   = 20000 #time steps
t = np.linspace(0,tmax,Nt+1) #time array
#time stamps of interest: [beginning, stage1, stage2: just after saturation (ts), stage2, stage2: ponding time (tp), stage3]
t_interest = [0,0.05,0.11224956732746585+0.0001,0.20,0.3003003003003003,tmax]

############################################################# 
#Functions for two-layered soils
#############################################################

#Flux in the saturated region qs
def qs(zu,zl):
    return ((m*p - 1)*(zu/z0-zl/z0))/((1-zu/z0)**(1-m*p)-(1-zl/z0)**(1-m*p))
    
#Inverting for moisture content at the surface from rainfall
def func_R2phil_U(theta_0):
    return (f_Cm(np.array([phi_0]),m,phi_0)[0]*f_Cn(np.array([theta_0]),np.array([phi_0]),s_gr,s_wr,n))[0] - R 
theta_0 = opt.fsolve(func_R2phil_U,0.01)[0]  
fcbyR = 1/f_Cn(np.array([theta_0]),np.array([phi_0]),s_gr,s_wr,n)[0]
  
#Stage 2: ODEs for upper and lower shock locations y[0] and y[1] respectively: Equations (26) and (27) from paper 
def rhs_stage2(t, y): 
    return [ (qs(y[0],y[1]) -f_Cm(np.array([phi_0]),m,phi_0)*f_Cn(np.array([theta_0]),np.array([phi_0]),s_gr,s_wr,n))[0]/(phi_0*(1-y[0]/z0)**p*(1- s_gr - s_wr)*(1-(1/fcbyR)**(1/n)*(1-y[0]/z0)**(-m*p/n))), \
              qs(y[0],y[1])/(phi_0*(1-y[1]/z0)**p*(1-s_wr-s_gr))]  
    
#Stage 3: ODEs for upper and lower shock locations y[0] and y[1] respectively: Equations (26) and (27) from paper 
def rhs_stage3(t, y): 
    return [ 0, qs(y[0],y[1])/(phi_0*(1-y[1]/z0)**p*(1-s_wr-s_gr))] 
        
    
zs = [z0*(1-(1/fcbyR)**(1/(m*p)))] #dimensionless depth of complete saturation
ts = n/(-m*p + n*p + n)*phi_0*z0*fcbyR*(1 - s_wr - s_gr)*((1/fcbyR)**(1/n) - (1/fcbyR)**((p+1)/(m*p))) #dimensionless time of saturation
tnew = t - ts
tnew = tnew[tnew>0]
res = solve_ivp(rhs_stage2, (0, tnew[-1]), [-1e-15+zs[0],1e-15+zs[0]],t_eval=tnew)
tp = tnew[np.argwhere(res.y[0,:]<0)][0,0] + ts 

print("############################################################# \n")
print("10 (top) Volume fraction with depth")
print("############################################################# \n")

# First set up the figure, the axis
fig,([ax1,ax2,ax3,ax4,ax5,ax6]) = time_sequenced_figure(zbottom,0,0,1,t_interest)

xf = z0 - z0*(1 - (-m*p + n*p +n )/(n*fcbyR*phi_0*z0*(1 - s_wr - s_gr))*fcbyR**(1/n)*  t_interest[0])**(n/(-m*p + n*p + n))

S_w_analy_int = s_wr*np.ones((Nz,1))
S_w_analy_int[zc<=xf] = np.transpose([s_wr + (1 - s_gr - s_wr)*(1/fcbyR)**(1/n)*(1-zc[zc<=xf]/z0)**(-m*p/n)])

S_w_analy_int_combined = []
S_w_analy_int_combined = S_w_analy_int.copy()

ax1.fill_betweenx(zc,1, facecolor=red,label=r'$\phi_g$')
ax1.fill_betweenx(zc,(1-phi+phi*S_w_analy_int)[:,0], facecolor=blue,label=r'$\phi_w$')
ax1.fill_betweenx(zc,(1-phi)[:,0], facecolor=brown,label=r'$\phi_s$')

xf = z0 - z0*(1 - (-m*p + n*p +n )/(n*fcbyR*phi_0*z0*(1 - s_wr - s_gr))*fcbyR**(1/n)*  t_interest[1])**(n/(-m*p + n*p + n))
S_w_analy_int = s_wr*np.ones((Nz,1))
S_w_analy_int[zc<=xf] = np.transpose([s_wr + (1 - s_gr - s_wr)*(1/fcbyR)**(1/n)*(1-zc[zc<=xf]/z0)**(-m*p/n)])

S_w_analy_int_combined = np.hstack([S_w_analy_int_combined,S_w_analy_int])

ax2.fill_betweenx(zc,1, facecolor=red,label=r'$\phi_g$')
ax2.fill_betweenx(zc,(1-phi+phi*S_w_analy_int)[:,0], facecolor=blue,label=r'$\phi_w$')
ax2.fill_betweenx(zc,(1-phi)[:,0], facecolor=brown,label=r'$\phi_s$')

res = solve_ivp(rhs_stage2, (0, t_interest[2] - ts), [-1e-16+zs[0],1e-16+zs[0]],t_eval=[t_interest[2]-ts])

S_w_analy_int = (1-s_gr)*np.ones((Nz,1))
S_w_analy_int[zc<=res.y[0,0]] = np.transpose([s_wr + (1 - s_gr - s_wr)*(1/fcbyR)**(1/n)*(1-zc[zc<=res.y[0,0]]/z0)**(-m*p/n)])
S_w_analy_int[zc>=res.y[1,0]] = s_wr


S_w_analy_int_combined = np.hstack([S_w_analy_int_combined,S_w_analy_int])

ax3.fill_betweenx(zc,1, facecolor=red,label=r'$\phi_{gas}$')
ax3.fill_betweenx(zc,(1-phi+phi*S_w_analy_int)[:,0], facecolor=blue,label=r'$\phi_{water}$')
ax3.fill_betweenx(zc,(1-phi)[:,0], facecolor=brown,label=r'$\phi_{soil}$')


res = solve_ivp(rhs_stage2, (0, t_interest[3] - ts), [-1e-16+zs[0],1e-16+zs[0]],t_eval=[t_interest[3]-ts])
S_w_analy_int = (1-s_gr)*np.ones((Nz,1))
S_w_analy_int[zc<=res.y[0,0]] = np.transpose([s_wr + (1 - s_gr - s_wr)*(1/fcbyR)**(1/n)*(1-zc[zc<=res.y[0,0]]/z0)**(-m*p/n)])
S_w_analy_int[zc>=res.y[1,0]] = s_wr

S_w_analy_int_combined = np.hstack([S_w_analy_int_combined,S_w_analy_int])

ax4.fill_betweenx(zc,1, facecolor=red,label=r'$\phi_g$')
ax4.fill_betweenx(zc,(1-phi+phi*S_w_analy_int)[:,0], facecolor=blue,label=r'$\phi_w$')
ax4.fill_betweenx(zc,(1-phi)[:,0], facecolor=brown,label=r'$\phi_s$')

res = solve_ivp(rhs_stage2, (0, tnew[-1]), [-1e-16+zs[0],1e-16+zs[0]],t_eval=[t_interest[4]-ts])

S_w_analy_int = (1-s_gr)*np.ones((Nz,1))
S_w_analy_int[zc<=res.y[0,0]] = np.transpose([s_wr + (1 - s_gr - s_wr)*(1/fcbyR)**(1/n)*(1-zc[zc<=res.y[0,0]]/z0)**(-m*p/n)])
S_w_analy_int[zc>=res.y[1,0]] = s_wr

S_w_analy_int_combined = np.hstack([S_w_analy_int_combined,S_w_analy_int])

ax5.fill_betweenx(zc,1, facecolor=red,label=r'$\phi_{gas}$')
ax5.fill_betweenx(zc,(1-phi+phi*S_w_analy_int)[:,0], facecolor=blue,label=r'$\phi_{water}$')
ax5.fill_betweenx(zc,(1-phi)[:,0], facecolor=brown,label=r'$\phi_{soil}$')

ytop = res.y[0,0]
ybot = res.y[1,0]
    
res = solve_ivp(rhs_stage3, (0, t_interest[5]-tp), [ytop,ybot],t_eval=[t_interest[5]-tp])
S_w_analy_int = (1-s_gr)*np.ones((Nz,1))
S_w_analy_int[zc<=res.y[0,0]] = np.transpose([s_wr + (1 - s_gr - s_wr)*(1/fcbyR)**(1/n)*(1-zc[zc<=res.y[0,0]]/z0)**(-m*p/n)])
S_w_analy_int[zc>=res.y[1,0]] = s_wr

S_w_analy_int_combined = np.hstack([S_w_analy_int_combined,S_w_analy_int])

ax6.fill_betweenx(zc,1, facecolor=red,label=r'$\phi_{gas}$')
ax6.fill_betweenx(zc,(1-phi+phi*S_w_analy_int)[:,0], facecolor=blue,label=r'$\phi_{water}$')
ax6.fill_betweenx(zc,(1-phi)[:,0], facecolor=brown,label=r'$\phi_{soil}$')


ax3.set_title(r'''%.2f ($t_s'$)'''%t_interest[2], fontsize='medium')
ax5.set_title(r'''%.2f ($t_p'$)'''%tp, fontsize='medium')
ax1.legend(loc='lower left', shadow=False, fontsize='medium')
plt.xlabel("Volume fractions $\phi$", fontsize='medium')
plt.savefig(f"./Figures/{simulation_name}_CL{theta_0}.pdf")

print("############################################################# \n")
print("10 (bottom) Saturation with depth")
print("############################################################# \n")

fig,([ax1,ax2,ax3,ax4,ax5,ax6]) = time_sequenced_figure(zbottom,0,-0.05,1.05,t_interest)
ax1.plot(S_w_analy_int_combined[:,0],zc , c = 'k',linestyle='--')
ax2.plot(S_w_analy_int_combined[:,1],zc , c = 'k',linestyle='--')
ax3.plot(S_w_analy_int_combined[:,2],zc , c = 'k',linestyle='--')
ax4.plot(S_w_analy_int_combined[:,3],zc , c = 'k',linestyle='--')
ax5.plot(S_w_analy_int_combined[:,4],zc , c = 'k',linestyle='--')
ax6.plot(S_w_analy_int_combined[:,5],zc , c = 'k',linestyle='--')

ax3.set_title(r'''%.2f ($t_s'$)'''%(ts+0.005), fontsize='medium')
ax5.set_title(r'''%.2f ($t_p'$)'''%tp, fontsize='medium')
plt.xlabel("Water saturation $s_w$", fontsize='medium')
plt.savefig(f"./Figures/swvsZpanelshock_{simulation_name}_phi0{theta_0}.pdf")

print("############################################################# \n")
print("11 Dimensionless infiltration rate with dimensionless time")
print("############################################################# \n")

print('R/fc:',R,''', tp':''',tp,''', zs':''',z0,''', ts':''',ts,'\n')

fig = plt.figure(figsize=(8,8) , dpi=100)
t0 = np.linspace(0,tp,1000)
qs0 = (f_Cm(np.array([phi_0]),m,phi_0)*f_Cn(np.array([theta_0]),np.array([phi_0]),s_gr,s_wr,n))* np.ones_like(t0)
t1 = np.linspace(tp,5,10000)
res= solve_ivp(rhs_stage3, (t1[0],t1[-1]), [ytop,ybot],t_eval=t1)
y  = res.y
qs1 = qs(y[0],y[1])

T  = np.concatenate([t0,t1])
QS = np.concatenate([qs0,qs1])
plot = plt.plot(T,QS/f_Cm(np.array([phi_0]),m,phi_0),c='r')
plt.ylabel(r'''$I(t')/f_c$''', fontsize='medium')
plt.xlabel(r'''$t'$''', fontsize='medium')
plt.xticks(fontsize='medium')
plt.yticks(fontsize='medium')
plt.xlim([np.min(T), 5])
plt.ylim([0,round(np.max(QS),1)+0.03])
plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
plt.savefig(f"./Figures/ICvsT_phi0array_{simulation_name}_phi0{phi_0}.pdf")
