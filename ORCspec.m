function ORCstruct = ORCspec(dTc, dTh_req, mdot_h, Th_in, Th_out, Tc, eff_pmp, eff_tbn, D)
% FEEG6013 Group Design Project, 2021-2022
% Group 19
%
% Created by Michael
%
%
% Calculates ORC performance as well as thermodynamic properties at each
% state
%
%
%% Constant (non-Tmax dependent) properties
% Min cycle temperature, via driving temp [K]
Tmin = Tc + dTc;

% State properties on saturation line
data_c1  = data_r600a_sat(Tmin, 0);  % State 1
data_c4  = data_r600a_sat(Tmin, 1);  % State 4 (sat v)

% Coolant properties and heat rejection required
data_h = data_water(Th_in, Th_out);
Qin = mdot_h*data_h.cp*abs(Th_in - Th_out);  % [W]


%% Calc cycle properties, find Tmax & Tpp via binary search
% Initial guess of max ORC temp
Tmax = Th_out-1;
% Set upper/lower bounds for iteration loop
Tmax_lwr = Tmin;
Tmax_upr = Th_in;
% Iterate to find spec'd dTh
while true
	% State properties on saturation line
	data_c2f = data_r600a_sat(Tmax, 0);  % State 2f
	data_c3v = data_r600a_sat(Tmax, 1);  % State 3v

	% Subcooled liq state properties, via Gibb's eq w/ ds = 0
	h2  = data_c1.h + 1E6*(data_c2f.p - data_c1.p)/data_c1.rho;
	T2 = Tmin + (h2 - data_c1.h)/data_c1.cp;

	% Get x3 from data_c4.s, data_c2f.s, data_c3v.s (isentropic: s3 = data_c4.s), then get state 3 props
	x3 = (data_c4.s - data_c2f.s)/(data_c3v.s - data_c2f.s);
	data_c3 = data_r600a_sat(Tmax, x3);

	% ORC mass flow req'd, via energy balance [kg/s]
	mdot = Qin/(data_c3.h - h2);

	% Calc pinch-point temp, pinch-point T-difference [K]
	Tpp = Th_out + ( mdot*(data_c2f.h - h2) )/(mdot_h*data_h.cp);
	dTh = Tpp - Tmax;


	% Compare actual dTh with required dTh
	ddTh = max(dTh_req, dTh) - min(dTh_req, dTh);
	% Converged!
	if ddTh < 0.0001
		break
	end
	% Not converged...
	if dTh > dTh_req						
		Tmax_lwr = Tmax;				% If dTh too high, set Tmax lower bound,
		Tmax = (Tmax_lwr + Tmax_upr)/2;	% then avg bounds for new Tmax guess
	else
		Tmax_upr = Tmax;				% If dTh too low, set Tmax upper bound,
		Tmax = (Tmax_lwr + Tmax_upr)/2;	% then avg bounds for new Tmax guess
	end
end
clear ddTh Tmax_lwr Tmax_upr

% Heat in to preheater & evaporator sections, respectively [W]
Qin_evp = mdot_h*data_h.cp*(Th_in - Tpp);
Qin_ph  = mdot_h*data_h.cp*(Tpp - Th_out);

% Heat rejection req'd [W]
Qout = mdot*abs(data_c1.h - data_c4.h);


%% Pipe diameter required to maintain incompressible flow (if value not supplied)
if exist('D', 'var') == 0
	Ma_max = 0.28;		% Max incompressible Mach number
	
	% Minimum rho*c
	rhoc_min = ...
		min([data_c1.rho*data_c1.c, ...					% State 1
 			data_c2f.rho*data_c2f.c, ...				% State 2
 			data_c3.rho*data_c3.c, ...					% State 3
 			data_c4.rho*data_c4.c, ...					% State 4
 			]);
	
	% Required diameter [m]
	D  = sqrt( (4*mdot)/(pi*rhoc_min*Ma_max) );
end
% Cross-section area [m^2]
Ac = pi*(D/2)^2;


%% Calculate cycle performance
% Pump and turbine power [W]
pwr_pmp  = (1/eff_pmp)*mdot*abs(h2 - data_c1.h);
pwr_tbne = eff_tbn*mdot*abs(data_c4.h - data_c3.h); 

% Net power out [W]
pwr_net = pwr_tbne - pwr_pmp;

% Thermal efficiency [%]
eff = 100*(pwr_tbne - pwr_pmp)/Qin;


%% ORC T-s cycle points (for plotting)
%     [1,         2,         2f,         3,          4,          1   ]
T_n = [Tmin,      T2,        Tmax,       Tmax,       Tmin,       Tmin];
s_n = [data_c1.s, data_c1.s, data_c2f.s, data_c4.s, data_c4.s, data_c1.s];

% Coolant endpoints
Th_n = [Th_in,     Tpp,         Th_out];
sh_n = [data_c4.s, data_c2f.s,  data_c1.s];
% Heatsink (sea) endpoints
Tc_n = [Tc,        Tc];
sc_n = [data_c1.s, data_c4.s];

%% Put all into specified ORC structure
ORCstruct.dTh      = dTh;
ORCstruct.dTc      = dTc;
ORCstruct.Tpp      = Tpp;
ORCstruct.Tmax     = Tmax;
ORCstruct.Tmin     = Tmin;
ORCstruct.pmax     = data_c2f.p; 
ORCstruct.pmin     = data_c1.p;
ORCstruct.x1	   = 0;			% Obsolete parameter
ORCstruct.x3	   = x3; 
ORCstruct.mdot	   = mdot; 
ORCstruct.Ac	   = Ac; 
ORCstruct.D		   = D; 
ORCstruct.pwr_pmp  = pwr_pmp; 
ORCstruct.pwr_tbne = pwr_tbne; 
ORCstruct.pwr_net  = pwr_net;
ORCstruct.eff	   = eff;
ORCstruct.Qin	   = Qin; 
ORCstruct.Qout     = Qout;
ORCstruct.Qin_ph   = Qin_ph;
ORCstruct.Qin_evp  = Qin_evp;
ORCstruct.T2       = T2;
ORCstruct.h1       = data_c1.h;
ORCstruct.h2       = h2;
ORCstruct.h2f      = data_c2f.h;
ORCstruct.h3       = data_c3.h;
ORCstruct.h4       = data_c4.h;

% For plotting
ORCstruct.plotOpt.Tpp     = Tpp;
ORCstruct.plotOpt.pwr_net = pwr_net;
ORCstruct.plotOpt.x3      = x3;
ORCstruct.plotTs.T        = T_n;
ORCstruct.plotTs.s        = s_n;
ORCstruct.plotTs.Th       = Tpp;
ORCstruct.plotTs.Th_n     = Th_n;
ORCstruct.plotTs.sh_n     = sh_n;
ORCstruct.plotTs.Tc_n     = Tc_n;
ORCstruct.plotTs.sc_n     = sc_n;
end