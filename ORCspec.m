function ORCstruct = ORCspec(dTc, dTh_req, mdot_h, Th_in, Th_out, Tc, eff_pmp, eff_tbn)
% Calculate mass flow and pipe area requirements for ORC

%% Constant (non-Tmax dependent) properties
% Min cycle temperature, via driving temp [K]
Tmin = Tc + dTc;

% State properties on saturation line
[rho1,  ~, cp_c1, ~,~,~,~, pmin,    c1,  h1,  s1]  = data_r600a_sat(Tmin, 0);  % State 1
[rho4,  ~, ~,     ~,~,~,~, ~,       c4,  h4,  s4]  = data_r600a_sat(Tmin, 1);  % State 4 (sat v)

[~, cp_h] = data_water(Th_in, Th_out);
Qin = mdot_h*cp_h*abs(Th_in - Th_out);


%% Calc cycle properties, find Tmax & Tpp via binary search
% Initial guess of max ORC temp
Tmax = Th_out-1;
% Set upper/lower bounds for iteration loop
Tmax_lwr = Tmin;
Tmax_upr = Th_in;
% Iterate to find spec'd dTh
while true
	% State properties on saturation line
	[rho2f, ~,~,~,~,~,~, pmax, c2f, h2f, s2f] = data_r600a_sat(Tmax, 0);  % State 2f
	[~,     ~,~,~,~,~,~, ~,    ~,   ~,   s3v] = data_r600a_sat(Tmax, 1);  % State 3v
	% Subcooled liq state properties, via Gibb's eq w/ ds = 0
	h2  = h1 + 1E6*(pmax - pmin)/rho1;
	T2 = Tmin + (h2 - h1)/cp_c1;

	% Get x3 from s4, s2f, s3v (isentropic: s3 = s4), then get state 3 props
	x3 = (s4 - s2f)/(s3v - s2f);
	[rho3, ~,~,~,~,~,~,~, c3, h3] = data_r600a_sat(Tmax, x3);

	% ORC mass flow req'd, via energy balance [kg/s] (Source: SESM2017 ch9)
	mdot = Qin/(h3 - h2);

	% Calc pinch-point temp, pinch-point T-difference
	Tpp = Th_out + ( mdot*(h2f - h2) )/(mdot_h*cp_h);
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

% Heat rejection req'd [W]
Qout = mdot*abs(h1 - h4);


%% Pipe area required to maintain incompressible flow
Ma_max = 0.27;		% Max incompressible Mach number

% Minimum rho*c
rhoc_min = ...
	min([rho1*c1, ...					% State 1
 		rho2f*c2f, ...					% State 2
 		rho3*c3, ...					% State 3
 		rho4*c4, ...					% State 4
 		]);

% Required diameter [m] and cross-section area [m^2]
D  = sqrt( (4*mdot)/(pi*rhoc_min*Ma_max) );
Ac = pi*(D/2)^2;


%% Calculate cycle performance
% Pump and turbine power [W]
pwr_pmp  = (1/eff_pmp)*mdot*abs(h2 - h1);
pwr_tbne = eff_tbn*mdot*abs(h4 - h3); 

% Net power out
pwr_net = pwr_tbne - pwr_pmp;

% Thermal efficiency [%]
eff = 100*(pwr_tbne - pwr_pmp)/Qin;


%% ORC T-s cycle points
%     [1,    2,      2f,       3,        4,    1   ]
T_n = [Tmin, T2, Tmax, Tmax, Tmin, Tmin];
s_n = [s1,   s1, s2f,  s4,   s4,   s1];

% Coolant endpoints
Th_n = [Th_in, Th_out];
sh_n = [s4,     s1];

%% Put all into specified ORC structure
ORCstruct.dTh      = dTh;
ORCstruct.dTc      = dTc;
ORCstruct.Tpp      = Tpp;
ORCstruct.Tmax     = Tmax;
ORCstruct.Tmin     = Tmin;
ORCstruct.pmax     = pmax; 
ORCstruct.pmin     = pmin;
ORCstruct.x1       = 0;			% Legacy parameter
ORCstruct.x3       = x3; 
ORCstruct.mdot     = mdot; 
ORCstruct.Ac       = Ac; 
ORCstruct.D        = D; 
ORCstruct.pwr_pmp  = pwr_pmp; 
ORCstruct.pwr_tbne = pwr_tbne; 
ORCstruct.pwr_net  = pwr_net;
ORCstruct.Qin      = Qin; 
ORCstruct.Qout     = Qout;
ORCstruct.eff      = eff;
% For plotting
ORCstruct.plotOpt.Tpp     = Tpp;
ORCstruct.plotOpt.pwr_net = pwr_net;
ORCstruct.plotOpt.x3      = x3;
ORCstruct.plotTs.T        = T_n;
ORCstruct.plotTs.s        = s_n;
ORCstruct.plotTs.Th       = Tpp;
ORCstruct.plotTs.Th_n     = Th_n;
ORCstruct.plotTs.sh_n     = sh_n;
end