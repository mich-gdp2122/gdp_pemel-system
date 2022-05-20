% FEEG6013 Group Design Project, 2021-2022
% Group 19
%
% Created by Michael
%
%
% Temperature-dependent thermal property vectors of gases (H2, O2)
% Obtained from NIST Chemistry WebBook: https://webbook.nist.gov/chemistry/
%
%
%% Temperature vector
gasData.T_i = 273.16:5:423.16;  % [K]

%% H2 gas data
% sp. enthalpy [kJ/kg]
gasData.h2_h  = [3575.5	3646.6	3717.8	3789.1	3860.5	3931.9	4003.5	...
				 4075.2	4146.9	4218.7	4290.6	4362.6	4434.5	4506.6	...
				 4578.7	4650.8	4723	4795.2	4867.4	4939.7	5012	...
				 5084.3	5156.6	5229	5301.4	5373.8	5446.2	5518.6	...
				 5591	5663.5	5736];
% cp [kJ/(kg*K)]
gasData.h2_cp = [14.197	14.223	14.246	14.268	14.288	14.306	14.324	...
				 14.339	14.354	14.367	14.380	14.391	14.402	14.411	...
				 14.420	14.428	14.435	14.442	14.448	14.454	14.459	...
				 14.464	14.468	14.472	14.475	14.479	14.482	14.484	...
				 14.487	14.489	14.491];
% viscosity [uPa*s]
gasData.h2_mu = [8.3772	8.4830	8.5882	8.6929	8.7970	8.9005	9.0035	...
				 9.1060	9.2080	9.3094	9.4104	9.5109	9.6110	9.7106	...
				 9.8097	9.9085	10.007	10.105	10.202	10.299	10.396	...
				 10.492	10.588	10.684	10.779	10.874	10.968	11.063	...
				 11.156	11.250	11.343];
% therm. conductivity [mW/(m*K)]
gasData.h2_k  = [173.45	175.98	178.47	180.95	183.39	185.81	188.21	...
				 190.59	192.95	195.28	197.59	199.89	202.17	204.43	...
				 206.67	208.89	211.10	213.30	215.48	217.64	219.80	...
				 221.93	224.06	226.17	228.28	230.37	232.45	234.52	...
				 236.58	238.63	240.67];

%% O2 gas data
% sp. enthalpy [kJ/kg]
gasData.o2_h  = [248.09	252.67	257.25	261.84	266.43	271.02	275.61	...
				 280.21	284.82	289.43	294.04	298.66	303.28	307.91	...
				 312.54	317.18	321.82	326.47	331.13	335.79	340.46	...
				 345.13	349.81	354.50	359.20	363.90	368.60	373.32	...
				 378.04	382.77	387.51];
% cp [kJ/(kg*K)]
gasData.o2_cp = [0.91655	0.91708	0.91765	0.91827	0.91893	0.91963	0.92037	...
				 0.92116	0.92198	0.92285	0.92375	0.92470	0.92568	0.92670	...
				 0.92775	0.92884	0.92997	0.93113	0.93232	0.93353	0.93478	...
				 0.93606	0.93736	0.93869	0.94005	0.94142	0.94282	0.94424	...
				 0.94567	0.94713	0.94860];
% viscosity [uPa*s]
gasData.o2_mu = [19.142	19.428	19.711	19.993	20.273	20.551	20.827	...
				 21.101	21.373	21.644	21.913	22.180	22.446	22.710	...
				 22.972	23.232	23.491	23.749	24.005	24.259	24.512	...
				 24.764	25.014	25.263	25.510	25.756	26		26.244	...
				 26.486	26.726	26.966];
% therm. conductivity [mW/(m*K)]
gasData.o2_k  = [24.346	24.749	25.150	25.549	25.946	26.341	26.735	...
				 27.127	27.517	27.905	28.292	28.677	29.060	29.442	...
				 29.822	30.201	30.578	30.953	31.327	31.700	32.071	...
				 32.441	32.810	33.177	33.542	33.907	34.270	34.631	...
				 34.992	35.351	35.709];