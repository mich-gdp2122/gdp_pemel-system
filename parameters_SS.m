%% PLACE INPUT PARAMETERS IN THIS FILE %%

%%% Load constants and gas properties
const_file;
data_gas_file;

%% Ambient conditions.16
amb.T     = 20 + 273.15;       % Ambient temperature  [K]
amb.p     = 1  * 101325;       % Ambient air pressure [Pa]
amb.phi   = 0.60;              % Ambient humidity     [0-1]
amb.T_sea = 10 + 273.15;	   % Ocean background temperature [K]

%% Stack specifications
pemel.N_stk = 8;		   % No. stacks				  []
pemel.N_cel = 30;          % Number of cells in stack []
pemel.A_cel = 0.1*100^2;   % Cell active area         [cm^2]

% Membrane
%mem.A   = pemel.A_cel;   % Membrane area           [m^2]
%mem.thk = 1E-4;          % Membrane thickness      [m]
%mem.M   = 1.1;          % Membrane dry molar mass [kg/mol]
%mem.rho = 1.98E-3;       % Membrane dry density    [kg/cm^3]

% Cooling channels
clch.N   = 60;            % No. cooling tubes per BP plate       []
clch.L   = 0.2;           % Cooling tube length                  [m]
clch.W   = 0.001;		  % Cooling tube width					 [m]
clch.Ac  = clch.W^2;      % Cooling tube cross-section area      [m^2]
clch.Prm = 4*clch.W;      % Cooling tube cross-section perimeter [m]

% Bipolar plate
bp.L   = 0.5;			% BP plate length [m]
bp.W   = 0.2;			% BP plate width [m]
bp.thk = 0.002294;      % BP plate thickness b/tw heat source & channel [m]
bp.cp  = 535.71;        % BP plate Sp. heat (Cp)               [J/(kg*K)]
bp.k   = 20.233;        % BP plate heat conductivity           [W/(m*K)]

% Process channels
prch.N   = 4;			% No. anode/cathode channels each per cell []
prch.W   = 0.8E-3;		% Process channel width
prch.Ac  = prch.W^2;    % Process channel cross-section area       [m^2]
prch.Prm = 4*prch.W;    % Process channel cross-section perimeter  [m]
%prch.L   = 0.01;       % Process channel length                   [m]
%prch.vol = 0.1;        % Process channel volume                   [m^3]


%% Stack operating conditions
% Inputs
pemel.i = 2.2219;	    % Nominal current density [A/cm^2]

% Stack
pemel.p_ca  = 20 * 100000;  % Nominal cathode pressure [bar -> Pa]
pemel.p_an  = amb.p;	    % Nominal anode pressure   [Pa]
pemel.T_stk = 80 + 273.15; %interp1(pemel.i_i, pemel.T_i, pemel.i, 'makima', 'extrap');  % Nominal stack temperature [K]
pemel.V_stk = 60; %pemel.N_cel*interp1(pemel.i_i, pemel.V_i, pemel.i, 'makima', 'extrap');  % Nominal stack voltage [V]
%pemel.q     = 8955.5;  % Stack heat flux density [W/m^2]

% Process water
h2o.stoich    = 1.2;		% Stoichiometric ratio
h2o.mdot_stk  = 12;			% Nominal mass flow rate per stack [kg/s]
h2o.T_stk_in  = 345.7;		% Stack inlet temperature [K]

% Hydrogen
h2.mdot_stk = 6.84E-4;		% Nominal mass flow rate per stack [kg/s]

% Coolant
clnt.mdot_stk  = 0.9;	 % Nominal mass flow rate per stack [kg/s]
clnt.T_stk_in  = 339.15; %interp1(pemel.i_i, pemel.TIn_clnt_i, pemel.i, 'makima', 'extrap');  % Coolant inlet temperature @ stack [K] 
clnt.T_stk_out = 345.15;  % Coolant outlet temperature [K]
clnt.p_stk_in  = 2742.5 + amb.p; %interp1(pemel.i_i, pemel.pIn_clnt_i, pemel.i, 'makima', 'extrap'); % Coolant stack inlet pressure [Pa]
clnt.dp_stk    = 2681; % interp1(pemel.i_i, pemel.dp_clnt_i, pemel.i, 'makima', 'extrap');  % Coolant stack pressure drop [Pa]
clnt.p_stk_out = clnt.p_stk_in - clnt.dp_stk;  % Coolant outlet pressure [Pa]


%% Balance-of-Plant
% Efficiencies
BoP.eff_pmp = 0.90;  % Coolant pump efficiency []
%BoP.eff_fan = 0.90;  % Fan efficiency []


%% Preheater
ph.L = 0.3;	 %% Pipe length

%% Organic Rankine Cycle
%%%Cooling system model
Cooli.T_initial= 345.15; %Input temperature [K]
ORC.mdot       = 1.7;    % Nominal mass flow rate [kg/s]

% Heat exchanger
HX.L           = 4.5;      % Coolant-ORC Heat exchanger length (m)

% Condenser
Con.Fluid_V    = 5;      %Condenser fluid volume [m^3]
Con.PortA_A    = 0.01;   %Condenser port A area [m^2]
Con.PortB_A    = 0.01;   %Condenser port A area [m^2]

% Thermal 
Con.Area       = 1e-4;   %Condenser thermal resistance area [m^2]
Con.Heat_Coeff = 20;     %Condenser thermal resistance heat coefficient [W/(m^2*K)]
Con.Coolant_T  = 283.15;     %Condenser cooling temperature [K]

% Power system
Shaft.speed    = 3600;   %Shaft speed [rpm]


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Derived Parameters
%%%%  DO NOT PUT INPUT PARAMETERS HERE! (put them above this section)  %%%%

%pemel.Q_clt  = pemel.q*(clch.Prm*clch.L);   % Heat flux transfer to single tube [W]
clnt.dT_stk    = 5.98;

% Current density -> current
pemel.I       = pemel.i*pemel.A_cel;	% i -> I [A/cm^2 -> A/m^2 -> A]
%pemel.I_i     = pemel.i_i*pemel.A_cel;  % i -> I [A/cm^2 -> A/m^2 -> A]

% Single cooling channel geometries
clch.Vol = clch.Ac*clch.L;	% Volume [m^3]
clch.As  = clch.Prm*clch.L;	% Surface area [m^2]

% Specs for total no. cells overall 
pemel.totN_cel = pemel.N_cel*pemel.N_stk;  % Total no. cells overall
clch.N_tot     = pemel.totN_cel*clch.N;    % Total no. cool tubes overall
prch.N_tot     = pemel.totN_cel*prch.N;    % Total no. anode/cathode channels (each) overall
clch.Ac_tot    = clch.N_tot*clch.Ac;	   % Cool tube cross-section area [m^2]
clch.Prm_tot   = clch.N_tot*clch.Prm;	   % Cool tube cross-section perimeter [m]
clch.Vol_tot   = clch.N_tot*clch.Vol;	   % Cool tube volume [m^3]
clch.As_tot    = clch.N_tot*clch.As;	   % Cool tube surface area [m^2]
prch.Ac_tot    = prch.N_tot*prch.Ac;	   % Process channel cross-section area [m^2]
%pemel.V_tot    = pemel.N_stk*pemel.V_stk;  % Overall voltage [V]

% Hydraulic diameters [m]
clch.Dh = 4*clch.Ac/clch.Prm;  % Cooling channels
prch.Dh = 4*prch.Ac/prch.Prm;  % Process channels
Con.D   = 2*sqrt(Con.PortA_A/pi);	% ORC condenser diameter (circular pipe assumed)

% Total mass flow rates (for total no. cells overall)
clnt.mdot_tot = pemel.totN_cel*clnt.mdot_stk;

% Process fluid mass flow rates
h2o.mdot_reac_tot = pemel.totN_cel*const.M_h2o*pemel.I/(2*const.F);  % Total h2o consumed in reaction [kg/s]
h2o.mdot_in_tot   = h2o.stoich*h2o.mdot_reac_tot;					 % Total h2o inlet mass flow rate [kg/s]
h2o.mdot_out_tot  = h2o.mdot_in_tot - h2o.mdot_reac_tot;			 % Total h2o outlet mass flow rate [kg/s]
h2.mdot_reac_tot  = pemel.totN_cel*const.M_h2*pemel.I/(2*const.F);   % Total h2 mass produced [kg/s]


%mdot_h2_24h = 473; % kg/24hr
%mdot_h2_s   = round(mdot_h2_24h/86400, 3, 'significant');

%n_h2_tot_calc = 0.83*h2.mdot_reac_tot/const.M_h2;
%n_h2_tot_set = mdot_h2_s/const.M_h2;

% Stack-level mass flow conversion
%pemel.H2out_stk_i = pemel.N_cel*pemel.H2out_i;  % H2 out [kg/s]

%%%%  DO NOT PUT INPUT PARAMETERS HERE! (put them above this section)  %%%%