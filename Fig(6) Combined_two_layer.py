#Infiltration in a two-layered initially-dry soil (Combined)
#Mohammad Afzal Shadab and Marc Hesse
#Date modified: 08/05/21

from supporting_tools import *

#parameters
simulation_name = f'combined_two_layered'
m = 3      #Cozeny-Karman coefficient for numerator K = K0 (1-phi_i)^m
n = 2      #Corey-Brooks coefficient krw = krw0 * se(sw)^n
s_wr = 0.0 #Residual water saturation
s_gr = 0.0 #Residual gas saturation
phi_U = 0.5; phi_L = 0.2 #Porosity of upper and lower layer
theta_L = s_gr*phi_L #setting the saturation in domain
zsurface = -1    #dimensionless surface height
z0 = 0.0         #dimensionless location of jump

#temporal discretization
tmax = 20   #time scaling with respect to z0/fc
Nt   = 20000#time steps
t = np.linspace(0,tmax,Nt+1)

print("############################################################# \n")
print("6(a) Varying dimensionless rainfall rate R/fc")
print("############################################################# \n")

R_array = np.linspace(0.2,0.9,8)  #Array of dimensionless rainfall rates

Blues = cm.get_cmap('Blues', len(R_array)+1) #Color palette

fig = plt.figure(figsize=(8,8) , dpi=100)
#looping over the arrays off rainfall
for count, R in enumerate(R_array):
    #Inverting for moisture content at the surface from rainfall
    def func_R2phil_U(theta_U):
        return (f_Cm(np.array([phi_U]),m,phi_U)[0]*f_Cn(np.array([theta_U]),np.array([phi_U]),s_gr,s_wr,n))[0] - R 
    theta_U = opt.fsolve(func_R2phil_U,0.01)[0]    

    #Stage 2: ODEs for upper and lower shock locations y[0] and y[1] respectively: Equations (26) and (27) from paper 
    def rhs_stage2(t, y): 
        return [((y[0]-y[1])/(y[0]/f_Cm(np.array([phi_U]),m,phi_U)-y[1]/f_Cm(np.array([phi_L]),m,phi_U))-f_Cm(np.array([phi_U]),m,phi_U)*f_Cn(np.array([theta_U]),np.array([phi_U]),s_gr,s_wr,n))[0]/(phi_U*(1-s_gr)-theta_U), \
                ((y[0]-y[1])/(y[0]/f_Cm(np.array([phi_U]),m,phi_U)-y[1]/f_Cm(np.array([phi_L]),m,phi_U))-f_Cm(np.array([phi_L]),m,phi_U)*f_Cn(np.array([theta_L]),np.array([phi_L]),s_gr,s_wr,n))[0]/(phi_L*(1-s_gr)-theta_L)]
    
    #Flux in the saturated region qs
    def qs(x):
        return (x-1)/(x/f_Cm(np.array([phi_U]),m,phi_U) - 1/f_Cm(np.array([phi_L]),m,phi_U))
    
    #Relation (32) - (35) to find the shock speed ratio k
    aa = phi_U/phi_L*(1-theta_U/phi_U-s_gr)/(1-theta_L/phi_L-s_gr)*(1- f_Cm(np.array([phi_L]),m,phi_U) * f_Cn(np.array([theta_L]),np.array([phi_L]),s_gr,s_wr,n) / f_Cm(np.array([phi_U]),m,phi_U))
    bb =-phi_U/phi_L*(1-theta_U/phi_U-s_gr)/(1-theta_L/phi_L-s_gr)*(1- f_Cm(np.array([phi_L]),m,phi_U) * f_Cn(np.array([theta_L]),np.array([phi_L]),s_gr,s_wr,n) / f_Cm(np.array([phi_L]),m,phi_U)) - (1 - f_Cm(np.array([phi_U]),m,phi_U) * f_Cn(np.array([theta_U]),np.array([phi_U]),s_gr,s_wr,n) / f_Cm(np.array([phi_U]),m,phi_U))
    cc = 1 - f_Cm(np.array([phi_U]),m,phi_U) * f_Cn(np.array([theta_U]),np.array([phi_U]),s_gr,s_wr,n) / f_Cm(np.array([phi_L]),m,phi_U)
    k = (-bb-np.sqrt(bb**2-4*aa*cc))/(2*aa)
    
    #Calculating the shock speeds
    s_U_analy =((qs(k) - f_Cm(np.array([phi_U]),m,phi_U)*f_Cn(np.array([theta_U]),np.array([phi_U]),s_gr,s_wr,n))[0]/(phi_U*(1-s_gr-theta_U/phi_U)))
    s_L_analy =((qs(k) - f_Cm(np.array([phi_L]),m,phi_U)*f_Cn(np.array([theta_L]),np.array([phi_L]),s_gr,s_wr,n))[0]/(phi_L*(1-theta_L/phi_L-s_gr)))
    
    Sf = ((f_Cm(np.array([phi_U]),m,phi_U)[0]*f_Cn(np.array([theta_U]),np.array([phi_U]),s_gr,s_wr,n))[0] / (phi_U*(theta_U/phi_U - s_wr))) #initial front speed
    ts = (z0 - zsurface)/Sf             #dimensionless time of saturation
    tp = ts + (zsurface-z0)/s_U_analy   #dimensionless time of ponding
    print('R/fc:',R,''', tp':''',tp,''', zs':''',z0,''', ts':''',ts,'\n')
    res = solve_ivp(rhs_stage2, (0, tp-ts), [-1e-14,1e-15],t_eval=[tp-ts])
    ytop = res.y[0,0]    #dimensionless shock location at the time of ponding
    ybot = res.y[1,0]    #dimensionless shock location at the time of ponding
    
    #Stage 3: ODEs for upper and lower shock locations y[0] and y[1] respectively: Equations (26) and (27) from paper 
    def rhs_stage3(t, y): 
        return [0, \
                ((y[0]-y[1])/(y[0]/f_Cm(np.array([phi_U]),m,phi_U)-y[1]/f_Cm(np.array([phi_L]),m,phi_U))-f_Cm(np.array([phi_L]),m,phi_U)*f_Cn(np.array([theta_L]),np.array([phi_L]),s_gr,s_wr,n))[0]/(phi_L*(1-s_gr)-theta_L)]
    
    #Dimensionless infiltration rate with dimensionless time
    t0 = np.linspace(0,ts,1000)
    qs0 = (f_Cm(np.array([phi_U]),m,phi_U)*f_Cn(np.array([theta_U]),np.array([phi_U]),s_gr,s_wr,n))* np.ones_like(t0)
    t1 = np.linspace(ts,tp,1000)
    qs1 = (f_Cm(np.array([phi_U]),m,phi_U)*f_Cn(np.array([theta_U]),np.array([phi_U]),s_gr,s_wr,n))*np.ones_like(t1)
    t2 = np.linspace(tp,10,10000)
    res= solve_ivp(rhs_stage3, (t2[0],t2[-1]), [ytop,ybot],t_eval=t2)
    y  = res.y
    qs2 = (y[0]-y[1])/(y[0]/f_Cm(np.array([phi_U]),m,phi_U) - y[1]/f_Cm(np.array([phi_L]),m,phi_U))
    
    T  = np.concatenate([t0,t1,t2])
    QS = np.concatenate([qs0,qs1,qs2])
    plot = plt.plot(T,QS/f_Cm(np.array([phi_U]),m,phi_U),c=Blues(count+1),label='$R/f_c=$%0.2f'%R)
    plt.plot(t2[0],qs2[0],'k.')

plt.hlines(f_Cm(np.array([phi_L]),m,phi_U), np.min(T), np.max(T), colors=gray, linestyles='--')
plt.legend(loc='best', shadow=False, fontsize='medium')
plt.ylabel(r'''$I(t')/f_c$''', fontsize='medium')
plt.xlabel(r'''$t'$''', fontsize='medium')
plt.xticks(fontsize='medium')
plt.yticks(fontsize='medium')
plt.xlim([np.min(T), 5])
plt.ylim([0,round(np.max(QS),1)+0.1])
plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
plt.savefig(f"ICvsT_shock_{simulation_name}_CL{theta_U}CR{theta_L}_phiL{phi_U}phiR{phi_L}.pdf")

print("\n")
print("############################################################# \n")
print("6(b) Varying porosity of lower soil: phi_L")
print("############################################################# \n")

phi_L_array = [0.01,0.05,0.10,0.15,0.20,0.25,0.288,0.35,0.40]

Reds = cm.get_cmap('Reds', len(phi_L_array)+1)

fig = plt.figure(figsize=(8,8) , dpi=100)

#looping over the arrays of lower layer porosity
for count, phi_L in enumerate(phi_L_array):

    # %% Define derivative function
    def rhs_stage2(t, y): 
        return [((y[0]-y[1])/(y[0]/f_Cm(np.array([phi_U]),m,phi_U)-y[1]/f_Cm(np.array([phi_L]),m,phi_U))-f_Cm(np.array([phi_U]),m,phi_U)*f_Cn(np.array([theta_U]),np.array([phi_U]),s_gr,s_wr,n))[0]/(phi_U*(1-s_gr)-theta_U), \
                ((y[0]-y[1])/(y[0]/f_Cm(np.array([phi_U]),m,phi_U)-y[1]/f_Cm(np.array([phi_L]),m,phi_U))-f_Cm(np.array([phi_L]),m,phi_U)*f_Cn(np.array([theta_L]),np.array([phi_L]),s_gr,s_wr,n))[0]/(phi_L*(1-s_gr)-theta_L)]
    
    aa = phi_U/phi_L*(1-theta_U/phi_U-s_gr)/(1-theta_L/phi_L-s_gr)*(1- f_Cm(np.array([phi_L]),m,phi_U) * f_Cn(np.array([theta_L]),np.array([phi_L]),s_gr,s_wr,n) / f_Cm(np.array([phi_U]),m,phi_U))
    bb =-phi_U/phi_L*(1-theta_U/phi_U-s_gr)/(1-theta_L/phi_L-s_gr)*(1- f_Cm(np.array([phi_L]),m,phi_U) * f_Cn(np.array([theta_L]),np.array([phi_L]),s_gr,s_wr,n) / f_Cm(np.array([phi_L]),m,phi_U)) - (1 - f_Cm(np.array([phi_U]),m,phi_U) * f_Cn(np.array([theta_U]),np.array([phi_U]),s_gr,s_wr,n) / f_Cm(np.array([phi_U]),m,phi_U))
    cc = 1 - f_Cm(np.array([phi_U]),m,phi_U) * f_Cn(np.array([theta_U]),np.array([phi_U]),s_gr,s_wr,n) / f_Cm(np.array([phi_L]),m,phi_U)
    
    k = (-bb-np.sqrt(bb**2-4*aa*cc))/(2*aa)
    
    def qs(x):
        return (x-1)/(x/f_Cm(np.array([phi_U]),m,phi_U) - 1/f_Cm(np.array([phi_L]),m,phi_U))
    
    s_U_analy =((qs(k) - f_Cm(np.array([phi_U]),m,phi_U)*f_Cn(np.array([theta_U]),np.array([phi_U]),s_gr,s_wr,n))[0]/(phi_U*(1-s_gr-theta_U/phi_U)))
    s_L_analy =((qs(k) - f_Cm(np.array([phi_L]),m,phi_U)*f_Cn(np.array([theta_L]),np.array([phi_L]),s_gr,s_wr,n))[0]/(phi_L*(1-theta_L/phi_L-s_gr)))
    
    Sf = ((f_Cm(np.array([phi_U]),m,phi_U)[0]*f_Cn(np.array([theta_U]),np.array([phi_U]),s_gr,s_wr,n))[0] / (phi_U*(theta_U/phi_U - s_wr)))
    ts = (z0 - zsurface)/Sf
    tp = ts + (zsurface-z0)/s_U_analy
    print('phi_L:',phi_L,''', tp':''',tp,''', zs':''',z0,''', ts':''',ts,'\n')  
    res = solve_ivp(rhs_stage2, (0, tp-ts), [-1e-14,1e-15],t_eval=[tp-ts])
    
    ytop = res.y[0,0]
    ybot = res.y[1,0]
    
    def rhs_stage3(t, y): 
        return [0, \
                ((y[0]-y[1])/(y[0]/f_Cm(np.array([phi_U]),m,phi_U)-y[1]/f_Cm(np.array([phi_L]),m,phi_U))-f_Cm(np.array([phi_L]),m,phi_U)*f_Cn(np.array([theta_L]),np.array([phi_L]),s_gr,s_wr,n))[0]/(phi_L*(1-s_gr)-theta_L)]
    
    #Infiltration rate with time
    t0 = np.linspace(0,ts,1000)
    qs0 = (f_Cm(np.array([phi_U]),m,phi_U)*f_Cn(np.array([theta_U]),np.array([phi_U]),s_gr,s_wr,n))* np.ones_like(t0)
    
    t1 = np.linspace(ts,tp,1000)
    qs1 = (f_Cm(np.array([phi_U]),m,phi_U)*f_Cn(np.array([theta_U]),np.array([phi_U]),s_gr,s_wr,n))*np.ones_like(t1)

    t2 = np.linspace(tp,50,10000)
    res= solve_ivp(rhs_stage3, (t2[0],t2[-1]), [ytop,ybot],t_eval=t2)
    y  = res.y
    
    qs2 = (y[0]-y[1])/(y[0]/f_Cm(np.array([phi_U]),m,phi_U) - y[1]/f_Cm(np.array([phi_L]),m,phi_U))
    
    T  = np.concatenate([t0,t1,t2])
    QS = np.concatenate([qs0,qs1,qs2])
    plot = plt.plot(T,QS/f_Cm(np.array([phi_U]),m,phi_U),c=Reds(count+1),label='$\phi_l=$%0.2f'%phi_L)
    plt.plot(t2[0],qs2[0],'k.')

plt.legend(loc='best', shadow=False, fontsize='medium')
plt.xlim([np.min(T), 5])
plt.ylim([0,round(np.max(QS),1)+0.1])
plt.ylabel(r'''$I(t')/f_c$''', fontsize='medium')
plt.xlabel(r'''$t'$''', fontsize='medium')
plt.xticks(fontsize='medium')
plt.yticks(fontsize='medium')
plt.xlim([np.min(T), 1.2])
plt.savefig(f"ICvsT_phiR_shock_{simulation_name}_CL{theta_U}CR{theta_L}_phiL{phi_U}phiR{phi_L}.pdf")
