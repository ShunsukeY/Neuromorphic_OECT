import matplotlib.pyplot as plt
import numpy as np

#parameters -----------------------------------------------------------------------------
#time (based on pulse width = 1)
delta_t = 0.1			#delta t for PPD
tau_G = 0.3				#gate response time const
tau_Din = 0.3 			#drain response injection
tau_Dex = 0.6	        #drain response extaction

#current
i_G0 = 0.7					#gate current amplitude
i_Gspike= 0.8				#gate spike current amplitude
i_Dbase= -1.16				#drain current baseline at VG=0

#--------------------------------------------------------------------------------------------

#setting arrays----------------------------------------------------------------------------
# define result array -- row 0:time 1:VG, 2:IG(ig), 3:ID(id)
result = np.zeros((4,100000))  

# set time row (row0)
tm = np.linspace(0, 1e2, 100000, dtype=float) 
result[0] = tm

#set V_G row (row1)
vg_tmp= [0]*98000 											 # off state
vg_tmp_on= [1] * 1000 										#on state
tm_2in=int(2000+(delta_t*1000)) 						#set the timing of 2nd pulse
vg_tmp=np.insert(vg_tmp, 1000, vg_tmp_on) 		#insert 1st pulse
vg_tmp=np.insert(vg_tmp, tm_2in, vg_tmp_on)		#insert 2nd pulse

result[1] = vg_tmp 													#merge

#--------------------------------------------------------------------------------------------


#calculation of IG---------------------------------------------------------------------------

#1st injection (at row no. 1000)
sw_1in = np.array([0]*1000 + [1]*1000 + [0]*98000) 				#window function
ig_1in = sw_1in * i_G0 * np.exp(-(tm-1)/tau_G)		#main decay
ig_1in[1000]=ig_1in[1000]+i_Gspike						#add spike

#1st extraction (at row no. 2000)
sw_1ex = np.array([0]*2000 + [1]*(tm_2in-2000)+ [0]*(100000-tm_2in))
ig_1ex = sw_1ex * i_G0 * np.exp(-(tm-2)/tau_G)
ig_1ex[2000]=ig_1ex[2000]+i_Gspike

#2nd injection (at tm_2in in column no.)
sw_2in = np.array([0]*tm_2in + [1]*1000 + [0]*(99000-tm_2in))
ig_2in = sw_2in * i_G0 * np.exp(-((tm-(2+delta_t))/tau_G))
ig_2in[tm_2in]=ig_2in[tm_2in]+i_Gspike

#2nd extraction (at tm_2in+1000 in column no.)
sw_2ex = np.array([0]*(tm_2in+1000) + [1]*(99000-tm_2in))
ig_2ex = sw_2ex * i_G0 * np.exp(-(tm-(3+delta_t))/tau_G)
ig_2ex[tm_2in+1000]=ig_2ex[tm_2in+1000]+i_Gspike

ig = ig_1in-ig_1ex+ig_2in-ig_2ex						#merge

result[2] = ig
#--------------------------------------------------------------------------------------------

#calculation of ID---------------------------------------------------------------------------

#set baseline + gate current
id_base = np.array([i_Dbase]*100000)
id_tmp = id_base

#1st injection (Note: window function was shared from IG)
amp_1in = np.array([(0-id_tmp[(1000-1)])]*100000) * sw_1in    #amplitude of decay/rise
id_1in = amp_1in *(1- np.exp(-(tm-1)/tau_Din))	#main decay
id_tmp = id_tmp + id_1in 

#1st extraction 
amp_1ex = np.array([(id_tmp[(2000-1)]-i_Dbase)]*100000) * sw_1ex
id_1ex = amp_1ex * np.exp(-(tm-2)/tau_Dex)
id_tmp = id_tmp + id_1ex

#2nd injection 
amp_2in = np.array([(0-id_tmp[(tm_2in-1)])]*100000) * sw_2in    
id_2in = amp_2in *(1- np.exp(-((tm-(2+delta_t))/tau_Din))) + sw_2in  * (id_tmp[(tm_2in-1)]-i_Dbase)
id_tmp = id_tmp + id_2in 

#2nd extraction 
amp_2ex = np.array([(id_tmp[(tm_2in+1000-1)]-i_Dbase)]*100000) * sw_2ex
id_2ex = amp_2ex * np.exp(-(tm-(3+delta_t))/tau_Dex)
id_tmp = id_tmp + id_2ex

result[3] =    id_tmp - ig

#calc A1 and A2
A1 = -(result[3,1000]-i_Dbase)
A2 = -(result[3,tm_2in]-i_Dbase)
print(1-(A2/A1))
#--------------------------------------------------------------------------------------------




#prepare plots------------------------------------------------------------------------------

#plot settings-----------------------------
plt.rcParams['font.family'] ='Arial'       
plt.rcParams['xtick.direction'] = 'in'              
plt.rcParams['ytick.direction'] = 'in'              
plt.rcParams['xtick.major.width'] = 1.0             
plt.rcParams['ytick.major.width'] = 1.0             
plt.rcParams['font.size'] = 12                          
plt.rcParams['axes.linewidth'] = 1.0   
#------------------------------------------     

fig = plt.figure(figsize=(8,8))

gs = fig.add_gridspec(3, height_ratios=[0.3,1,1])

ax1 = fig.add_subplot(gs[0])
ax2 = fig.add_subplot(gs[1])
ax3 = fig.add_subplot(gs[2])

ax1.plot(result[0], result[1], color='black')
ax1.set_xlim(0.5, 5)
ax1.set_ylim(-0.2, 1.2)
ax1.set_ylabel('$V_{G}$')

ax2.plot(result[0], result[2], color='blue', label = 'model')
ax2.set_xlim(0.5, 5)
ax2.set_ylabel('$I_{G}$')
ax2.legend()

ax3.plot(result[0], result[3], color='red', label = 'model')
ax3.set_xlim(0.5, 5)
ax3.set_xlabel('Time')
ax3.set_ylabel('$I_{D}$')
ax3.legend()

plt.tight_layout()

file_name = 'calcd_v_obsd.tif'
plt.savefig(file_name, format="tif")
        
plt.show()

#--------------------------------------------------------------------------------------------