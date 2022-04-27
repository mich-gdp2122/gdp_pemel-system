function ORCstruct = ORCspec(dTc, dTh, mdot_h, Th_in, Th_out, Tc, eff_pmp, eff_tbn)
% Calculate mass flow and pipe area requirements for ORC

%% Constant (non-Tmax dependent) properties
% Min cycle temperature, via driving temp [K]
Tmin = Tc + dTc;

% State properties on saturation line
[rho1,  ~, cp_c1, ~,~,~,~, pmin,    c1,  h1,  s1]  = data_r600a_sat(Tmin, 0);  % State 1
[rho4,  ~, ~,     ~,~,~,~, ~,       c4,  h4,  s4]  = data_r600a_sat(Tmin, 1);  % State 4 (sat v)

% Coolant cp [J/(kg*K)] and heat transfer req'd [W]
[~, cp_h, ~, ~, ~] = data_water(Th_in, Th_out);
Qin = mdot_h*cp_h*abs(Th_in - Th_out);


%% Pinch-point optimisation
% Iterate over Th range to find optimal pinch-point temp for max net pwr out
Tpp_i = linspace(min(Th_in,Th_out), max(Th_in,Th_out));
for i = 1:length(Tpp_i)
	%% Calculate top-side (Tmax dependent) properties
	% Max cycle temperature, via driving temp [K]
	Tmax(i) = Tpp_i(i) - dTh;
	
	% State properties on saturation line
	[rho2f, ~,~,~,~,~,~, pmax(i), c2a, h2f, s2f(i)] = data_r600a_sat(Tmax(i), 0);  % State 2f
	% Subcooled liq state properties, via Gibb's eq w/ ds = 0
	h2   = h1 + 1E6*(pmax - pmin)/rho1;
	T2(i) = Tmin + (h2 - h1)/cp_c1;
	% State properties in sat region, via isentropic processes
	[rho3, c3, h3, x3(i)] = isentrop(Tmax(i), s4, 1);   % State 3, via process 3->4

	% Sp. enthalpy changes [J/kg]
	dh_12 = h2 - h1;
	dh_23 = h3 - h2;
	dh_34 = h4 - h3;
	dh_41 = h1 - h4;
	
	
	%% Mass flow required to transfer coolant heat
	% ORC mass flow req'd, via energy balance [kg/s]
	mdot(i) = Qin/abs(dh_23);
	
	% Heat rejection req'd [W]
	Qout(i) = mdot(i)*abs(dh_41);
	

	%% Pipe area required to maintain incompressible flow
	Ma_max = 0.27;		% Max incompressible Mach number
	
	% Minimum rho*c
	rhoc_min = ...
		min([rho1*c1, ...					% State 1
	 		rho2f*c2a, ...					% State 2
	 		rho3*c3, ...					% State 3
	 		rho4*c4, ...					% State 4
	 		]);
	
	% Required diameter [m] and cross-section area [m^2]
	D(i)  = sqrt( (4*mdot(i))/(pi*rhoc_min*Ma_max) );
	Ac(i) = pi*(D(i)/2)^2;

	
	%% Calculate cycle performance
	% Pump and turbine power [W]
	pwr_pmp(i)  = (1/eff_pmp)*mdot(i)*abs(dh_12);
	pwr_tbne(i) = eff_tbn*mdot(i)*abs(dh_34); 

	% Net power out
	pwr_net(i) = pwr_tbne(i) - pwr_pmp(i);

	% Thermal efficiency [%]
	eff(i) = 100*(pwr_tbne(i) - pwr_pmp(i))/Qin;
end

%% Find optimal values
[pwr_net_max, i_opt] = max(pwr_net);
opt_Tpp      = Tpp_i(i_opt);
opt_Tmax     = Tmax(i_opt);
opt_pmax     = pmax(i_opt);
opt_x3       = x3(i_opt);
opt_mdot     = mdot(i_opt);
opt_Ac       = Ac(i_opt);
opt_D        = D(i_opt);
opt_pwr_pmp  = pwr_pmp(i_opt);
opt_pwr_tbne = pwr_tbne(i_opt);
opt_eff      = eff(i_opt);
opt_s2f      = s2f(i_opt);
opt_T2       = T2(i_opt);

%% Put all into specified ORC structure
ORCstruct.dTh      = dTh;
ORCstruct.dTc      = dTc;
ORCstruct.Tpp      = opt_Tpp;
ORCstruct.Tmax     = opt_Tmax;
ORCstruct.Tmin     = Tmin;
ORCstruct.pmax     = opt_pmax; 
ORCstruct.pmin     = pmin;
ORCstruct.x1       = 0;			% Legacy parameter
ORCstruct.x3       = opt_x3; 
ORCstruct.mdot     = opt_mdot; 
ORCstruct.Ac       = opt_Ac; 
ORCstruct.D        = opt_D; 
ORCstruct.pwr_pmp  = opt_pwr_pmp; 
ORCstruct.pwr_tbne = opt_pwr_tbne; 
ORCstruct.pwr_net  = pwr_net_max;
ORCstruct.Qin      = Qin; 
ORCstruct.Qout     = Qout;
ORCstruct.eff      = opt_eff;
% For plotting
ORCstruct.plot_opt.Tpp     = Tpp_i;
ORCstruct.plot_opt.pwr_net = pwr_net;
ORCstruct.plot_opt.x3      = x3;
ORCstruct.plot_Ts.T2       = opt_T2;
ORCstruct.plot_Ts.s1       = s1;
ORCstruct.plot_Ts.s2       = s1;
ORCstruct.plot_Ts.s2f      = opt_s2f;
ORCstruct.plot_Ts.s3       = s4;
ORCstruct.plot_Ts.s4       = s4;
end


%%
function [rhoB, cB, hB, xB] = isentrop(TB, sA, xA)
% Calc. outlet properties of isentropic process
	if xA == 1
		x = flip(0:01:1);	% Search from sat v line
	else
		x = 0:01:1;			% Search from sat L line otherwise
	end
	% Initial guess
	[~,~,~,~,~,~,~,~,~,~,sB_prev] = data_r600a_sat(TB, x(1));
	% Cycle thru vapour qualities until sp. entropy out = sp. entropy in
	for i = 2:length(x)
		[~,~,~,~,~,~,~,~,~,~, sB] = data_r600a_sat(TB, x(i));
		if (sB >= sA && sB_prev <= sA) || (sB <= sA && sB_prev >= sA)
			% Interp for xB between found sB vals of interest
			xB = interp1([sB_prev, sB], [x(i-1), x(i)], sA, 'linear');
			% Get output data
			[rhoB, ~,~,~,~,~,~,~, cB, hB, sB] = data_r600a_sat(TB, xB);
			break
		else
			% Otherwise set current sB val as prev then continue loop
			sB_prev = sB;
		end
	end
end