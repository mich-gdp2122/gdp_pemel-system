%% Derived Parameters
% Current density -> current
pemel.I       = pemel.i*pemel.A_cel;	% i -> I [A/cm^2 -> A/m^2 -> A]

% Single cooling, process channel geometries
clch.Vol = clch.Ac*clch.L;	% Cooling channel volume [m^3]
clch.As  = clch.Prm*clch.L;	% Cooling channel surface area [m^2]
prch.Vol = prch.Ac*prch.L;  % Process channel volume [m^3]

% Specs for total no. cells overall 
pemel.totN_cel = pemel.N_cel*pemel.N_stk;  % Total no. cells overall
clch.N_tot     = pemel.totN_cel*clch.N;    % Total no. cool tubes overall
prch.N_tot     = pemel.totN_cel*prch.N;    % Total no. anode/cathode channels (each) overall
clch.Ac_tot    = clch.N_tot*clch.Ac;	   % Cool tube cross-section area [m^2]
clch.Prm_tot   = clch.N_tot*clch.Prm;	   % Cool tube cross-section perimeter [m]
clch.Vol_tot   = clch.N_tot*clch.Vol;	   % Cool tube volume [m^3]
clch.As_tot    = clch.N_tot*clch.As;	   % Cool tube surface area [m^2]
prch.Ac_tot    = prch.N_tot*prch.Ac;	   % Process channel cross-section area [m^2]
%prch.Prm_tot   = prch.N_tot*prch.Prm;	   % Process channel cross-section perimeter [m]

% Stack channel hydraulic diameters [m]
clch.Dh_stk = 4*clch.Ac/clch.Prm;  % Cooling channels
prch.Dh_stk = 4*prch.Ac/prch.Prm;  % Process channels

% External pipe diameters (circular pipe assumed) [m]
prch.D = sqrt(4*prch.Ac_tot/pi);	% Process water
clch.D = sqrt(4*clch.Ac_tot/pi);	% Coolant
%Con.D  = sqrt(4*Con.PortA_A/pi);	% ORC condenser

% Total mass flow rates (for total no. cells overall) [kg/s]
clnt.mdot_tot = 1000*clnt.Vldot_ch*clch.N_tot;
clnt.mdot_stk = clnt.mdot_tot/pemel.N_stk;

% Process fluid total mass flow rates [kg/s]
h2o.mdot_reac_tot = pemel.totN_cel*const.M_h2o*pemel.I/(2*const.F);  % Total h2o consumed in reaction [kg/s]
h2o.mdot_in_tot   = h2o.stoich*h2o.mdot_reac_tot;					 % Total h2o inlet mass flow rate [kg/s]
h2o.mdot_out_tot  = h2o.mdot_in_tot - h2o.mdot_reac_tot;			 % Total h2o outlet mass flow rate [kg/s]
h2.mdot_reac_tot  = pemel.totN_cel*const.M_h2*pemel.I/(2*const.F);   % Total h2 mass produced [kg/s]

% Closed-loop feedwater base temperature
h2o.T0 = calc_T0(h2o.mdot_in_tot, h2o.mdot_out_tot, amb.T, h2o.T_stk_out);


%% Heat exchanger sizing & performance (NO preheat)
% ORC properties
ORC = ORCspec(ORC.dTc, ORC.dTh, clnt.mdot_tot, clnt.T_stk_out, clnt.T_stk_in, amb.T_sea, ...
	BoP.eff_pmp, BoP.eff_tbne);

% Heat rejector
[HX_rj.L, HX_rj.Rt, HX_rj.As, HX_rj.U] = ...
	HXsizer_rjct(clnt.mdot_tot, clch.D, amb.T_sea, clnt.T_stk_out, clnt.T_stk_in);

% Coolant/ORC fluid HX
[HX_ORC.L_h, HX_ORC.L_c, HX_ORC.Rt, HX_ORC.As, HX_ORC.U] = ...
	HXsizer_ORC(clnt.mdot_tot, ORC.mdot, ORC.D, clch.D, ...
	clnt.T_stk_out, clnt.T_stk_in, ORC.Tpp, ORC.Tmax, ORC.Tmin, ORC.x3);

% ORC condenser
[HX_cond.L, HX_cond.Rt, HX_cond.As, HX_cond.U] = ...
	HXsizer_cond(ORC.mdot, ORC.Qout, ORC.D, amb.T_sea, ORC.Tmin, ORC.x1);


%% Heat exchanger sizing & performance (preheat)
% Preheater
[HX_ph.L_h2o, HX_ph.L_clnt, HX_ph.Rt, HX_ph.As, HX_ph.U, clnt.Tout_PhCL] = ...
	HXsizer_PH(h2o.mdot_in_tot, clnt.mdot_tot, prch.D, clch.D, ...
	h2o.T0, h2o.T_stk_in, clnt.T_stk_out);

% ORC properties
ORCph = ORCspec(ORC.dTc, ORC.dTh, clnt.mdot_tot, clnt.Tout_PhCL, clnt.T_stk_in, amb.T_sea, ...
	BoP.eff_pmp, BoP.eff_tbne);

% Heat rejector
[HX_rjPh.L, HX_rjPh.Rt, HX_rjPh.As, HX_rjPh.U] = ...
	HXsizer_rjct(clnt.mdot_tot, clch.D, ...
	amb.T_sea, clnt.Tout_PhCL, clnt.T_stk_in);

% Coolant/ORC fluid HX
[HX_ORCph.L_h, HX_ORCph.L_c, HX_ORCph.Rt, HX_ORCph.As, HX_ORCph.U] = ...
	HXsizer_ORC(clnt.mdot_tot, ORCph.mdot, ORCph.D, clch.D, ...
	clnt.Tout_PhCL, clnt.T_stk_in, ORCph.Tpp, ORCph.Tmax, ORCph.Tmin, ORC.x3);

% ORC condenser length [m], thermal resistance [K/W]
[HX_condPh.L, HX_condPh.Rt, HX_condPh.As, HX_condPh.U] = ...
	HXsizer_cond(ORCph.mdot, ORCph.Qout, ORCph.D, amb.T_sea, ORC.Tmin, ORC.x1);

% Electric heater length (assumed equal to feedwater-side preheater length) [m]
htr.L = HX_ph.L_h2o;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function T0_h2o = calc_T0(mdot_in, mdot_out, T_amb, T_out)
% Calculate post-mixed feedwater temperature for closed-loop config	
	data_amb = data_water(T_amb, T_amb);
	data_out = data_water(T_out, T_out);
	T0_h2o  = ( (mdot_in - mdot_out)*data_amb.cp*T_amb + mdot_out*data_out.cp*T_out )...
			   /( (mdot_in - mdot_out)*data_amb.cp + mdot_out*data_out.cp);
end