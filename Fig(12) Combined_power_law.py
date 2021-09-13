#Infiltration in an initially-dry soil with power-law decay (Combined)
#Mohammad Afzal Shadab and Marc Hesse
#Date modified: 08/05/21

from supporting_tools import *

#parameters
simulation_name = f'combined_power_law'
m = 3      #Cozeny-Karman coefficient for numerator K = K0 (1-phi_i)^m
n = 2      #Corey-Brooks coefficient krw = krw0 * se(sw)^n
p  = 7.63  #index specified for power law porosity variation
s_wr = 0.0 #Residual water saturation
s_gr = 0.0 #Residual gas saturation
phi_0 = 0.5#Porosity at the surface
theta_L = s_gr*phi_0 #setting the saturation in domain
z0 = 1.0             #characteristic e-folding depth, surface is at z= 0

#temporal discretization
tmax = 10   #time scaling with respect to z0/fc
Nt   = 60000#time steps
t = np.linspace(0,tmax,Nt+1)

print("############################################################# \n")
print("(a) Varying dimensionless rainfall rate R/fc")
print("############################################################# \n")

R_array =[0.15,0.2,0.3,0.4,0.50,0.60,0.70,0.80] #Array of dimensionless rainfall rates

Blues = cm.get_cmap('Blues', len(R_array)+1) #Color palette

fig = plt.figure(figsize=(8,8) , dpi=100)
#looping over the arrays off rainfall
for count, R in enumerate(R_array):
    #Inverting for moisture content at the surface from rainfall
    def func_R2phil_U(theta_0):
        return (f_Cm(np.array([phi_0]),m,phi_0)[0]*f_Cn(np.array([theta_0]),np.array([phi_0]),s_gr,s_wr,n))[0] - R 
    theta_0 = opt.fsolve(func_R2phil_U,0.01)[0]    

    fcbyR = 1/f_Cn(np.array([theta_0]),np.array([phi_0]),s_gr,s_wr,n)[0]  #Dimensionless ratio of fc/R
    zs = [z0*(1-(1/fcbyR)**(1/(m*p)))] #dimensionless depth of complete saturation
    
    # %% Define derivative function
    ts = n/(-m*p + n*p + n)*phi_0*z0*fcbyR*(1 - s_wr - s_gr)*((1/fcbyR)**(1/n) - (1/fcbyR)**((p+1)/(m*p))) #dimensionless time of saturation
    tnew = t - ts
    tnew = tnew[tnew>0]

    #Stage 2: ODEs for upper and lower shock locations y[0] and y[1] respectively: Equations (26) and (27) from paper 
    def qs(zu,zl):
        return ((m*p - 1)*(zu/z0-zl/z0))/((1-zu/z0)**(1-m*p)-(1-zl/z0)**(1-m*p))
    
    def rhs_stage2(t, y): 
        return [ (qs(y[0],y[1]) -f_Cm(np.array([phi_0]),m,phi_0)*f_Cn(np.array([theta_0]),np.array([phi_0]),s_gr,s_wr,n))[0]/(phi_0*(1-y[0]/z0)**p*(1- s_gr - s_wr)*(1-(1/fcbyR)**(1/n)*(1-y[0]/z0)**(-m*p/n))), \
                  qs(y[0],y[1])/(phi_0*(1-y[1]/z0)**p*(1-s_wr-s_gr))]  
    
    res = solve_ivp(rhs_stage2,(0, tnew[-1]), [-1e-16+zs[0],1e-16+zs[0]],t_eval=tnew)
    y   = res.y
    qs_int = qs(y[0],y[1])
    
    tp = tnew[np.argwhere(y[0,:]<0)][0,0] + ts #dimensionless time of ponding: where upper shock < 0
    print('R/fc:',R,''', tp':''',tp,''', zs':''',zs[0],''', ts':''',ts,'\n')
    res = solve_ivp(rhs_stage2, (0, tp-ts), [-1e-16+zs[0],1e-16+zs[0]],t_eval=[tp-ts])
    ytop = res.y[0,0]    #dimensionless shock location at the time of ponding
    ybot = res.y[1,0]    #dimensionless shock location at the time of ponding
    
    #Stage 3: ODEs for upper and lower shock locations y[0] and y[1] respectively: Equations (26) and (27) from paper 
    def rhs_stage3(t, y): 
        return [ 0, qs(y[0],y[1])/(phi_0*(1-y[1]/z0)**p*(1-s_wr-s_gr))] 
            
    #Dimensionless infiltration rate with dimensionless time
    t0 = np.linspace(0,tp,1000)
    qs0 = (f_Cm(np.array([phi_0]),m,phi_0)*f_Cn(np.array([theta_0]),np.array([phi_0]),s_gr,s_wr,n))* np.ones_like(t0)
    t1 = np.linspace(tp,1,1000)
    res= solve_ivp(rhs_stage3, (t1[0],t1[-1]), [ytop,ybot],t_eval=t1)
    y  = res.y
    qs1 = qs(y[0],y[1])
    
    T  = np.concatenate([t0,t1])
    QS = np.concatenate([qs0,qs1])
    plot = plt.plot(T,QS/f_Cm(np.array([phi_0]),m,phi_0),c=Blues(count+1),label='$R/f_c=$%0.2f'%R)
    plt.plot(t1[0],qs1[0],'k.')

plt.legend(loc='best', shadow=False, fontsize='medium')
plt.ylabel(r'''$I(t')/f_c$''', fontsize='medium')
plt.xlabel(r'''$t'$''', fontsize='medium')
plt.xticks(fontsize='medium')
plt.yticks(fontsize='medium')
plt.xlim([np.min(T), 0.5])
plt.ylim([0,0.85])
plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
plt.savefig(f"ICvsT_Rarray_{simulation_name}_phi0{phi_0}.pdf")

print('\n')
print("############################################################# \n")
print("(b) Varying porosity at the surface: phi_0")
print("############################################################# \n")

phi_0_array =[0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8] #Array of porosity at the surface
R = 0.8 #Dimensionless rainfall rate R/fc

Reds = cm.get_cmap('Reds', len(phi_0_array)+1)

fig = plt.figure(figsize=(8,8) , dpi=100)
#looping over the array of surface porosity
for count, phi_0 in enumerate(phi_0_array):
    #Inverting for moisture content at the surface from rainfall
    def func_R2phil_U(theta_0):
        return (f_Cm(np.array([phi_0]),m,phi_0)[0]*f_Cn(np.array([theta_0]),np.array([phi_0]),s_gr,s_wr,n))[0] - R 
    theta_0 = opt.fsolve(func_R2phil_U,0.01)[0]    

    fcbyR = 1/f_Cn(np.array([theta_0]),np.array([phi_0]),s_gr,s_wr,n)[0]  #Dimensionless ratio of fc/R
    zs = [z0*(1-(1/fcbyR)**(1/(m*p)))] #dimensionless depth of complete saturation
    
    # %% Define derivative function
    ts = n/(-m*p + n*p + n)*phi_0*z0*fcbyR*(1 - s_wr - s_gr)*((1/fcbyR)**(1/n) - (1/fcbyR)**((p+1)/(m*p))) #dimensionless time of saturation
    tnew = t - ts
    tnew = tnew[tnew>0]

    #Stage 2: ODEs for upper and lower shock locations y[0] and y[1] respectively: Equations (26) and (27) from paper 
    def qs(zu,zl):
        return ((m*p - 1)*(zu/z0-zl/z0))/((1-zu/z0)**(1-m*p)-(1-zl/z0)**(1-m*p))
    
    def rhs_stage2(t, y): 
        return [ (qs(y[0],y[1]) -f_Cm(np.array([phi_0]),m,phi_0)*f_Cn(np.array([theta_0]),np.array([phi_0]),s_gr,s_wr,n))[0]/(phi_0*(1-y[0]/z0)**p*(1- s_gr - s_wr)*(1-(1/fcbyR)**(1/n)*(1-y[0]/z0)**(-m*p/n))), \
                  qs(y[0],y[1])/(phi_0*(1-y[1]/z0)**p*(1-s_wr-s_gr))]  
    
    res = solve_ivp(rhs_stage2,(0, tnew[-1]), [-1e-16+zs[0],1e-16+zs[0]],t_eval=tnew)
    y   = res.y
    qs_int = qs(y[0],y[1])
    
    tp = tnew[np.argwhere(y[0,:]<0)][0,0] + ts #dimensionless time of ponding: where upper shock < 0
    print('phi_0:',phi_0,''', tp':''',tp,''', zs':''',zs[0],''', ts':''',ts,'\n')
    res = solve_ivp(rhs_stage2, (0, tp-ts), [-1e-16+zs[0],1e-16+zs[0]],t_eval=[tp-ts])
    ytop = res.y[0,0]    #dimensionless shock location at the time of ponding
    ybot = res.y[1,0]    #dimensionless shock location at the time of ponding
    
    #Stage 3: ODEs for upper and lower shock locations y[0] and y[1] respectively: Equations (26) and (27) from paper 
    def rhs_stage3(t, y): 
        return [ 0, qs(y[0],y[1])/(phi_0*(1-y[1]/z0)**p*(1-s_wr-s_gr))] 
            
    #Dimensionless infiltration rate with dimensionless time
    t0 = np.linspace(0,tp,1000)
    qs0 = (f_Cm(np.array([phi_0]),m,phi_0)*f_Cn(np.array([theta_0]),np.array([phi_0]),s_gr,s_wr,n))* np.ones_like(t0)
    t1 = np.linspace(tp,1,1000)
    res= solve_ivp(rhs_stage3, (t1[0],t1[-1]), [ytop,ybot],t_eval=t1)
    y  = res.y
    qs1 = qs(y[0],y[1])
    
    T  = np.concatenate([t0,t1])
    QS = np.concatenate([qs0,qs1])
    plot = plt.plot(T,QS/f_Cm(np.array([phi_0]),m,phi_0),c=Reds(count+1),label='$\phi_0=$%0.2f'%phi_0)
    plt.plot(t1[0],qs1[0],'k.')

plt.legend(loc='best', shadow=False, fontsize='medium')
plt.ylabel(r'''$I(t')/f_c$''', fontsize='medium')
plt.xlabel(r'''$t'$''', fontsize='medium')
plt.xticks(fontsize='medium')
plt.yticks(fontsize='medium')
plt.xlim([np.min(T), 0.5])
plt.ylim([0,0.85])
plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
plt.savefig(f"ICvsT_phi0array_{simulation_name}_phi0{phi_0}.pdf")





