%% PLACE INPUT PARAMETERS IN THIS FILE %%

%%% Load constants and gas properties
const_file;
data_gas_file;

%% Ambient conditions
amb.T     = 20 + 273.15;       % Ambient temperature  [K]
amb.p     = 1  * 101325;       % Ambient air pressure [Pa]
amb.phi   = 0.60;              % Ambient humidity     [0-1]
amb.T_sea = 10 + 273.15;	   % Ocean background temperature [K]


%% Stack specifications
pemel.N_stk = 8;		   % No. stacks				  []
pemel.N_cel = 30;          % Number of cells in stack []
pemel.A_cel = 0.1*100^2;   % Cell active area         [cm^2]

% Cooling channels
clch.N   = 60;            % No. cooling tubes per BP plate       []
clch.L   = 0.25;           % Cooling tube length                  [m]
clch.W   = 0.001;		  % Cooling tube width					 [m]
clch.Ac  = clch.W^2;      % Cooling tube cross-section area      [m^2]
clch.Prm = 4*clch.W;      % Cooling tube cross-section perimeter [m]

% Bipolar plate
bp.L   = 0.5;			% BP plate length [m]
bp.W   = 0.2;			% BP plate width [m]
bp.thk = 0.002;			% BP plate thickness [m]
bp.rho = 4494;			% BP density [kg/m^3]
bp.cp  = 534;			% BP plate Sp. heat (Cp)               [J/(kg*K)]
bp.k   = 20.2;        % BP plate heat conductivity           [W/(m*K)]

% Process channels
prch.N   = 4;			% No. anode/cathode channels each per cell []
prch.W   = 0.8E-3;		% Process channel width
prch.Ac  = prch.W^2;    % Process channel cross-section area       [m^2]
prch.Prm = 4*prch.W;    % Process channel cross-section perimeter  [m]
prch.L   = 0.025;       % Process channel length                   [m]
%prch.vol = 0.1;        % Process channel volume                   [m^3]


%% Stack operating conditions
% Inputs
pemel.i		  = 2.2219;	      % Nominal current density [A/cm^2]
pemel.V_stk   = 60;			  % Nominal stack voltage [V]
pemel.T_stk   = 80 + 273.15;  % Nominal stack temperature [K]
pemel.p_ca    = 20 * 100000;  % Nominal cathode pressure [bar -> Pa]
pemel.p_an    = amb.p;	      % Nominal anode pressure   [Pa]
pemel.eff_stk = 0.8;		  % Stack efficiency

% Process water
h2o.stoich    = 1.2;		% Stoichiometric ratio
h2o.T_stk_in  = 340;		% Stack inlet temperature [K]
h2o.T_stk_out = 341;		% Stack outlet temperature [K]
%h2o.mdot_stk  = 12;		% Nominal mass flow rate per stack [kg/s]

% Hydrogen
%h2.mdot_stk = 6.84E-4;		% Nominal mass flow rate per stack [kg/s]

% Coolant
clnt.mdot_stk  = 0.9;	  % Nominal mass flow rate per stack [kg/s]
clnt.T_stk_in  = 339.15;  % Coolant inlet temperature @ stack [K] 
clnt.T_stk_out = 345.15;  % Coolant outlet temperature [K]
%clnt.p_stk_in  = 2742.5 + amb.p; % Coolant stack inlet pressure [Pa]
clnt.dp_stk    = 2681; % Coolant stack pressure drop [Pa]
%clnt.p_stk_out = clnt.p_stk_in - clnt.dp_stk;  % Coolant outlet pressure [Pa]


%% Balance-of-Plant
% Rated powers
BoP.PwrRt_htr = 2;	 % Elec. heater rated power [kW]

% Efficiencies
BoP.eff_pmp = 0.90;  % Coolant pump efficiency []

% Preheater
%ph.L = 1;	 %% Pipe length [m]

% Heat rejector
HX_rj.h_c = 10000;	 % Sea heat transf coeff. [W/(m^2*K)]

% Copper piping
%pipe.k   = 400;     % Thermal conductivity [W/(m*K)]
%pipe.thk = 0.5E-3;  % Pipe thickness [m]
%pipe.cp  = 385;     % Cp heat capacity [J/(kg*K)]
%pipe.rho = 8960;    % Density [kg/m^3]


%% Organic Rankine Cycle
%%%Cooling system model
Cooli.T_initial = 345.15; %Input temperature [K]
ORC.mdot        = 1.7;    % Nominal mass flow rate [kg/s]

% Heat exchanger
HX.L           = 10;      % Coolant-ORC Heat exchanger length (m)

% Condenser
Con.Fluid_V    = 0.5;      %Condenser fluid volume [m^3]
Con.PortA_A    = 0.02;   %Condenser port A area [m^2]
Con.PortB_A    = 0.02;   %Condenser port A area [m^2]

% Thermal 
Con.Area       = 16;   %Condenser thermal resistance area [m^2]
Con.Heat_Coeff = 3500;     %Condenser thermal resistance heat coefficient [W/(m^2*K)]
Con.Coolant_T  = 283.15;     %Condenser cooling temperature [K]

% Power system
Shaft.speed    = 3600;   %Shaft speed [rpm]


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%  DO NOT PUT INPUT PARAMETERS BELOW HERE!  (put them in above section)  %%%%

%% Derived Parameters

%pemel.Q_clt  = pemel.q*(clch.Prm*clch.L);   % Heat flux transfer to single tube [W]
clnt.dT_stk    = 5.98;

% Current density -> current
pemel.I       = pemel.i*pemel.A_cel;	% i -> I [A/cm^2 -> A/m^2 -> A]
%pemel.I_i     = pemel.i_i*pemel.A_cel;  % i -> I [A/cm^2 -> A/m^2 -> A]

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
%pemel.V_tot    = pemel.N_stk*pemel.V_stk;  % Overall voltage [V]

% Stack channel hydraulic diameters [m]
clch.Dh_stk = 4*clch.Ac/clch.Prm;  % Cooling channels
prch.Dh_stk = 4*prch.Ac/prch.Prm;  % Process channels

% External pipe diameters (circular pipe assumed) [m]
prch.D = sqrt(4*prch.Ac_tot/pi);	% Process water
clch.D = sqrt(4*clch.Ac_tot/pi);	% Coolant
Con.D  = sqrt(4*Con.PortA_A/pi);	% ORC condenser

% Total mass flow rates (for total no. cells overall) [kg/s]
clnt.mdot_tot = pemel.N_stk*clnt.mdot_stk;  % Coolant
%h2o.mdot_tot  = pemel.N_stk*h2o.mdot_stk;	% H2O
%h2.mdot_tot   = pemel.N_stk*h2.mdot_stk;    % H2

% Process fluid total mass flow rates [kg/s]
h2o.mdot_reac_tot = pemel.totN_cel*const.M_h2o*pemel.I/(2*const.F);  % Total h2o consumed in reaction [kg/s]
h2o.mdot_in_tot   = h2o.stoich*h2o.mdot_reac_tot;					 % Total h2o inlet mass flow rate [kg/s]
h2o.mdot_out_tot  = h2o.mdot_in_tot - h2o.mdot_reac_tot;			 % Total h2o outlet mass flow rate [kg/s]
h2.mdot_reac_tot  = pemel.totN_cel*const.M_h2*pemel.I/(2*const.F);   % Total h2 mass produced [kg/s]

% Total BP plate mass
bp.m = bp.rho * pemel.totN_cel*(bp.L*bp.W*bp.thk - (prch.N*prch.Vol + clch.N*clch.Vol));


%% Heat exchanger sizing
% Preheater lengths
[ph.L_h2o, ph.L_clnt] = HXsizer_PH(h2o.mdot_in_tot, clnt.mdot_tot, prch.D, clch.D, ...
						amb.T_sea, h2o.T_stk_in, clnt.T_stk_out);

% Heat rejector length & thermal resistance
[HX_rj.L, HX_rj.Rt, HX_rj.As, HX_rj.U] = HXsizer_rjct(clnt.mdot_tot, clch.D, ...
										 amb.T_sea, clnt.T_stk_out, clnt.T_stk_in);
%%%%  DO NOT PUT INPUT PARAMETERS HERE! (put them in first sectiom)  %%%%