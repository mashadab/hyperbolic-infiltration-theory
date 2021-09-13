from supporting_tools import *

#For more information about the parameters visit
#http://www.pc-progress.com/en/OnlineHelp/HYDRUS3/Hydrus.html?WaterFlowParameters.html

#############################################################
#Modified van-Genuchten model: Used in Hydrus
#############################################################

############ [upper , lower] layer parameters ############
theta_r = [0,0]           #residual moisture content [-]
theta_s = [0.43,0.1]      #saturated moisture content [-]
Alpha = [0.5,0.5]         #modified van-Genuchten model parameter [1/cm]
n = np.array([2.68,2.68]) #modified van-Genuchten model parameter [-]
Ks = [106.1,1.33447]      #saturated hydrualic conductivity [cm/s]
l = [0.5,0.5]             #modified van-Genuchten model parameter [-]
theta_m = [0.43,0.1]      #modified van-Genuchten model parameter [-]
theta_a = [0,0]           #modified van-Genuchten model parameter [-]
theta_k = [0.43,0.1]      #modified van-Genuchten model parameter [-]
Kk = [106.1,1.33447]      #modified van-Genuchten model parameter [cm/s]
m  = 1 - 1/n
theta_array = np.linspace(theta_r,theta_s,1000)
 
#Plotting
plt.figure(figsize=(8,14) , dpi=100)
plt.plot(Se(theta_array[:,0],theta_s[0],theta_r[0]),Kr(Kk[0],Ks[0],theta_array[:,0],theta_r[0],theta_s[0],theta_k[0],theta_m[0],theta_a[0],m[0]),'b-', alpha=0.5,label=r'Modified vG: Upper layer')
plt.plot(Se(theta_array[:,1],theta_s[1],theta_r[1]),Kr(Kk[1],Ks[1],theta_array[:,1],theta_r[1],theta_s[1],theta_k[1],theta_m[1],theta_a[1],m[1])*Ks[1]/Ks[0],'r-', alpha=0.5,label=r'Modified vG: Lower layer')


#############################################################
#Brooks-Corey model: Used in numerical & analytical results
#############################################################

############ Brooks-Corey model parameters ############
phi_U = 0.43
phi_L = 0.1
m = 3      #Cozeny-Karman coefficient for numerator K = K0 (1-phi_i)^m
n = 7.15306#Corey-Brooks coefficient krw = krw0 * sw^n
s_wr = 0.0 #Residual water saturation
s_gr = 0.0 #Residual gas saturation

s_w_new = np.linspace(s_wr,1-s_gr,1000)
theta_U = phi_U*s_w_new
flux_U  = f_Cm(np.array([phi_U]),m,phi_U)*f_Cn(theta_U,np.array([phi_U]),s_gr,s_wr,n)
theta_L = phi_L*s_w_new
flux_L  = f_Cm(np.array([phi_L]),m,phi_U)*f_Cn(theta_L,np.array([phi_L]),s_gr,s_wr,n)                      

#Plotting continued
plt.plot(s_w_new,flux_U,'b--',label=r'Brooks-Corey: Upper layer')
plt.plot(s_w_new,flux_L,'r--',label=r'Brooks-Corey: Lower layer')
plt.legend(loc='best', shadow=False, fontsize='medium')
if s_wr==0 and s_gr==0:
    plt.xlabel(r'${s_w}$', fontsize='medium')
else:
    plt.xlabel(r'$\frac{s_w - s_{wr}}{1 - s_{gr} - s_{wr}}$', fontsize='medium')
plt.ylabel(r'$f/f_c$', fontsize='medium')
plt.xlim([-0.01,1.01])
plt.xticks(fontsize='medium')
plt.yticks(fontsize='medium')
plt.ylim([-0.01,1.01])
plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
plt.axis('equal')
plt.savefig(f"./Figures/hydrus_fluxes.pdf")


