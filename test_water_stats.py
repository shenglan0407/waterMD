from water_stats import WaterStats

import mdtraj as md
import numpy as np
import matplotlib.pyplot as plt

data_path='/Users/shenglanqiao/Documents/GitHub/waterMD/data'
traj = md.load_trr(data_path+'/nvt-pr.trr', top = data_path+'/water-sol.gro')
print ('here is some info about the trajectory we are looking at:')
print traj
test = WaterStats(traj)
# 
# test.radial_dist([0.0,1.0])
# rs, g_R, g_err = test.rdf[0],test.rdf[1],test.rdf[2]
# fig3 = plt.figure()
# plt.errorbar(rs,g_R, yerr=g_err)
# plt.title('gn(r)')
# plt.xlabel('r (nm)')
# plt.ylabel('gn(r)')
# plt.ylim(0,max(g_R)*1.2)
# fig3.savefig('/Users/shenglanqiao/Documents/GitHub/waterMD/output/gn_r.png')
# plt.close(fig3)
# 
# Rs = np.linspace(0.2,0.4,10)
# 
# Sn_0R = []
# for R in Rs:
#     Sn_0R.append(test.struct_factor(0,R,1)[0])
# fig1 = plt.figure()
# plt.plot(Rs, Sn_0R)
# plt.title("Sn(0,R) with dt = 1.0 ps")
# plt.xlabel("R (nm)")
# plt.ylabel("Sn(0,R)") 
# fig1.savefig('/Users/shenglanqiao/Documents/GitHub/waterMD/output/Sn_0R.png')
# plt.close(fig1)

# R_water = 0.3
# R_max = 0.5
# dt = 8.0
# Qs = 2.*np.pi*np.linspace(0.0,1/R_water,60)
# test.estimate_struct_factor(Qs,R_max,dt)
# 
# fig = plt.figure()
# N_bar=[]
# for Q in Qs:
#     if Q == 0:
#         N_bar.append(0)
#     else:
#         
#         N_bar.append(-4./3.*np.pi*test.rho*3./Q**3.*(np.sin(Q*R_max)-Q*R_max*np.cos(Q*R_max)))
# plt.plot(np.array(Qs)*R_water/(2.*np.pi),N_bar)
# plt.xlabel('Q (unit of water radius)')
# fig.savefig('/Users/shenglanqiao/Documents/GitHub/waterMD/output/S_Q_BessTerm.png')
# plt.close(fig)
# 
# fig = plt.figure()

# Sn_QR = []
# for Q in Qs:
#     Sn_QR.append(test.struct_factor(Q,R_max,dt)[0])
# Qs,S_Q,S_Qerr=np.array(test.ssf[0]),test.ssf[1],test.ssf[2]
# plt.errorbar(Qs*R_water,S_Q,yerr=S_Qerr)
# plt.plot(Qs,Sn_QR,'go')
# plt.title("Estimate of S(Q) with dt = %.1f ps and R_inf = %.1f" % (dt,R_max))
# plt.xlabel("Q (rad)")
# plt.ylabel("S(Q)")
# plt.xlim(0,max(Qs*R_water))
# fig.savefig('/Users/shenglanqiao/Documents/GitHub/waterMD/output/S_Q.png')
# plt.close(fig)
# ts = np.linspace(1,10,10)
# 
# In_0tR1 = []
# In_0tR2 = []
# R1 = 0.3 # nm
# R2 = 0.5 # nm
# 
# for tt in ts:
#     In_0tR1.append(test.scat_func(0,R1,tt))
#     In_0tR2.append(test.scat_func(0,R2,tt))
#     
# fig2 = plt.figure()
# plt.plot(ts, In_0tR1-min(In_0tR1),label = "R = %.2f" % R1)
# plt.plot(ts, In_0tR2-min(In_0tR2),label = "R = %.2f" % R2)
# plt.legend()
# plt.title("In(0,t,R)")
# plt.xlabel("t (ps)")
# plt.ylabel("In(0,t,R)") 
# fig2.savefig('/Users/shenglanqiao/Documents/GitHub/waterMD/output/In_0tR.png')
# plt.close(fig2)
# 

# Sn_QR=[]
# for Q in Qs:
#     Sn_QR.append(test.struct_factor(Q,R,1)[0])
# 
# print Qs
# print Sn_QR
# plt.plot(Qs,Sn_QR,'o')
# plt.show()