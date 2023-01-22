# -*- coding: utf-8 -*-
"""
Created on Sat Dec 17 21:22:01 2022

@author: rps207
"""
#TA
marr=np.ma.getdata(meanPlumeAT)
marr=np.ma.getdata(meanNotPlumeAT)
marr=np.ma.getdata(meanAllAT)

month_numbers=np.argwhere(~np.isnan(marr))
month_numbers=np.squeeze(month_numbers)
marr = marr[~np.isnan(marr)]
vb=np.polyfit(month_numbers,marr,1)

fig49 = plt.figure(figsize=(21,24))
gs = fig49.add_gridspec(6, 1, hspace=0.28)
#plt.title("carbonate_timeseries_{0}".format(region))
fig49.tight_layout(pad=3.0)
ax11 = plt.subplot(6,1,1);
ax11.plot(month_numbers, marr, 'r--')
plt.plot(month_numbers,month_numbers*vb[0]+vb[1], 'b--')
plt.show()

#amazon
#TA all trend is 0.04563498757380958 per month
#TA plume trend is -0.11652206319181975 per month
#TA not plume trend is 0.06373828070270533 per month

#congo
#TA all trend is 0.04122747350109547 per month
#TA plume trend is 0.011212307472628542 per month
#TA not plume trend is 0.012127512511196741 per month


#DIC
marr=np.ma.getdata(meanPlumeDIC)
marr=np.ma.getdata(meanNotPlumeDIC)
marr=np.ma.getdata(meanAllDIC)

month_numbers=np.argwhere(~np.isnan(marr))
month_numbers=np.squeeze(month_numbers)
marr = marr[~np.isnan(marr)]
vb=np.polyfit(month_numbers,marr,1)

fig49 = plt.figure(figsize=(21,24))
gs = fig49.add_gridspec(6, 1, hspace=0.28)
#plt.title("carbonate_timeseries_{0}".format(region))
fig49.tight_layout(pad=3.0)
ax11 = plt.subplot(6,1,1);
ax11.plot(month_numbers, marr, 'r--')
plt.plot(month_numbers,month_numbers*vb[0]+vb[1], 'b--')
plt.show()

#amazon
#DIC -all trend is -0.010255446423501333 per month
#DIC -not plume trend is 0.04046963714256433 per month
#DIC -plume trend is -0.20115151712792484 per month

#congo
#DIC -all trend is -0.0013687577860983739 per month
#DIC -not plume trend is -0.007392177629331081 per month
#DIC -plume trend is 0.00047861665691566863 per month


#pH
marr=np.ma.getdata(meanAllpH_at)
month_numbers=np.argwhere(~np.isnan(marr))
month_numbers=np.squeeze(month_numbers)
marr = marr[~np.isnan(marr)]
vb=np.polyfit(month_numbers,marr,1)

fig49 = plt.figure(figsize=(21,24))
gs = fig49.add_gridspec(6, 1, hspace=0.28)
#plt.title("carbonate_timeseries_{0}".format(region))
fig49.tight_layout(pad=3.0)
ax11 = plt.subplot(6,1,1);
ax11.plot(month_numbers, marr, 'r--')
plt.plot(month_numbers,month_numbers*vb[0]+vb[1], 'b--')
plt.show()

#amazon

#pH -dic  plume trend is 1.6992646993276855e-05 per month
#pH -ta  plume trend is 3.490476325322784e-05 per month

#pH -dic nott  plume trend is -1.784295129644643e-05 per month
#pH -ta  not plume trend is -9.75742047632474e-07 per month

#pH -dic all trend is -1.1496168606554863e-05 per month
#pH -ta  all trend is 6.523664739178737e-06 per month

#Congo
#pH -dic all trend is -4.263716069229951e-05 per month
#pH -ta all trend is -3.9622635370333665e-05 per month

#pH -dic not plume trend is -5.539615320291365e-06 per month
#pH -ta not plume trend is -5.100397725194097e-06 per month

#pH -dic  plume trend is -5.6304578965206e-05 per month
#pH -ta  plume trend is -6.25814547532229e-05 per month



#pCO2
marr=np.ma.getdata(meanPlumepco2_at)
month_numbers=np.argwhere(~np.isnan(marr))
month_numbers=np.squeeze(month_numbers)
marr = marr[~np.isnan(marr)]
vb=np.polyfit(month_numbers,marr,1)

fig49 = plt.figure(figsize=(21,24))
gs = fig49.add_gridspec(6, 1, hspace=0.28)
#plt.title("carbonate_timeseries_{0}".format(region))
fig49.tight_layout(pad=3.0)
ax11 = plt.subplot(6,1,1);
ax11.plot(month_numbers, marr, 'r--')
plt.plot(month_numbers,month_numbers*vb[0]+vb[1], 'b--')
plt.show()

#amazon
#pCO2 - dic all trend is 0.015696478604151753 per month
#pCO2 - ta all trend is -0.006554126926569175 per month

#pCO2 - dic plume trend is -0.024018122735742657 per month
#pCO2 - ta plume trend is  -0.04658282490450927 per month

#pCO2 - dic not plume trend is 0.004639948747131351 per month
#pCO2 - ta not plume trend is  -0.024018122735742657 per month


#congo
#pCO2 - dic plume trend is 0.04603363633293212 per month
#pCO2 - ta plume trend is 0.055538444825516765 per month

#pCO2 - dic not plume trend is 0.0019428854745402186 per month
#pCO2 - ta not plume trend is 0.0006007066235220696 per month

#pCO2 - dic all trend is 0.03797113616520264 per month
#pCO2 - ta all trend is  0.03604262011702946 per month


