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

% Membrane
pemel.A_mem   = pemel.A_cel;   % Membrane area           [m^2]
pemel.thk_mem = 1E-4;          % Membrane thickness      [m]
%pemel.M_mem   = 1.1;          % Membrane dry molar mass [kg/mol]
pemel.rho_mem = 1.98E-3;       % Membrane dry density    [kg/cm^3]

% Cooling channels
pemel.N_clt   = 60;            % No. cooling tubes per BP plate       []
pemel.Ac_clt  = (1E-3)^2;      % Cooling tube cross-section area      [m^2]
pemel.Prm_clt = 4*(1E-3);      % Cooling tube cross-section perimeter [m]
pemel.L_clt   = 0.2;           % Cooling tube length                  [m]
pemel.W_clt   = 0.001;		   % Cooling tube width					  [m]

% Bipolar plate
pemel.L_bp    = 0.5;			% BP plate length [m]
pemel.W_bp    = 0.2;			% BP plate width [m]
pemel.thk_bp  = 0.002294;       % BP plate thickness b/tw heat source & channel [m]
pemel.cp_bp   = 535.71;         % BP plate Sp. heat (Cp)               [J/(kg*K)]
pemel.k_bp    = 20.233;         % BP plate heat conductivity           [W/(m*K)]

% Process channels
pemel.N_ch   = 4;				% No. anode/cathode channels each per cell []
pemel.Ac_ch  = (0.8E-3)^2;      % Process channel cross-section area       [m^2]
%pemel.Prm_ch = 0.002;          % Process channel cross-section perimeter  [m]
%pemel.L_ch   = 0.01;           % Process channel length                   [m]
%pemel.vol_ch = 0.1;            % Process channel volume                   [m^3]


%% Operating conditions
% Inputs
input.i = 2.2249;	    % Nominal current density [A/cm^2]

% Stack
pemel.p_ca  = 20 * 100000;  % Nominal cathode pressure [bar -> Pa]
pemel.p_an  = amb.p;	    % Nominal anode pressure   [Pa]
pemel.T_stk = 80 + 273.15; %interp1(pemel.i_i, pemel.T_i, input.i, 'makima', 'extrap');  % Nominal stack temperature [K]
pemel.V_stk = 60; %pemel.N_cel*interp1(pemel.i_i, pemel.V_i, input.i, 'makima', 'extrap');  % Nominal stack voltage [V]
%pemel.q     = 8955.5;  % Stack heat flux density [W/m^2]

% Process water
h2o.stoich    = 1.2;		% Stoichiometric ratio
h2o.mdot_stk  = 12;			% Nominal mass flow rate per stack [kg/s]
h2o.T_stk_in  = 345.7;		% Stack inlet temperature [K]

% Hydrogen
h2.mdot_stk = 6.84E-4;		% Nominal mass flow rate per stack [kg/s]

% Coolant
clnt.mdot_stk  = 0.9;			% Nominal mass flow rate per stack [kg/s]
clnt.T_stk_in  = 341.28; %interp1(pemel.i_i, pemel.TIn_clnt_i, input.i, 'makima', 'extrap');  % Coolant inlet temperature @ stack [K] 
clnt.dT_stk    = 5.98; %interp1(pemel.i_i, pemel.dT_clnt_i, input.i, 'makima', 'extrap');  % Coolant T-diff thru stack [K]
clnt.T_stk_out = clnt.T_stk_in + clnt.dT_stk;  % Coolant outlet temperature [K]
clnt.p_stk_in  = 2742.5 + amb.p; %interp1(pemel.i_i, pemel.pIn_clnt_i, input.i, 'makima', 'extrap'); % Coolant stack inlet pressure [Pa]
clnt.dp_stk    = 2681; % interp1(pemel.i_i, pemel.dp_clnt_i, input.i, 'makima', 'extrap');  % Coolant stack pressure drop [Pa]
clnt.p_stk_out = clnt.p_stk_in - clnt.dp_stk;  % Coolant outlet pressure [Pa]

%% BoP efficiencies
BoP.eff_pmp  = 0.90;  % Coolant pump efficiency []
BoP.eff_fan  = 0.90;  % Fan efficiency []


%% Organic Rankine Cycle
%%%Cooling system model
Cooli.T_initial= 353.15; %Input temperature [K]

%%%Condenser
Con.Fluid_V    = 5;      %Condenser fluid volume [m^3]
Con.PortA_A    = 0.01;   %Condenser port A area [m^2]
Con.PortB_A    = 0.01;   %Condenser port A area [m^2]

%Thermal 
Con.Area       = 1e-4;   %Condenser thermal resistance area [m^2]
Con.Heat_Coeff = 20;     %Condenser thermal resistance heat coefficient [W/(m^2*K)]
Con.Coolant_T  = 283.15;     %Condenser cooling temperature [K]

%%%Power system
Shaft.speed    = 3600;   %Shaft speed [rpm]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Derived Parameters
%%%%  DO NOT PUT INPUT PARAMETERS HERE! (put them above this section)  %%%%

%pemel.Q_clt  = pemel.q*(pemel.Prm_clt*pemel.L_clt);   % Heat flux transfer to single tube [W]
pemel.Dh_clt = 4*pemel.Ac_clt/pemel.Prm_clt;  % Hydraulic diameter [m]

% Current density -> current
input.I       = input.i*pemel.A_cel;	% i -> I [A/cm^2 -> A/m^2 -> A]
%pemel.I_i     = pemel.i_i*pemel.A_cel;  % i -> I [A/cm^2 -> A/m^2 -> A]

% Single cooling channel geometries
pemel.Vol_clt = pemel.Ac_clt*pemel.L_clt;	% Volume [m^3]
pemel.As_clt  = pemel.Prm_clt*pemel.L_clt;	% Surface area [m^2]

% Params for total no. cells per stack
pemel.N_clt_stk      = pemel.N_cel*pemel.N_clt;    % No. cooling tubes in stack
pemel.N_ch_stk       = pemel.N_cel*pemel.N_ch;	   % No. process channels in stack
pemel.totAc_clt_stk  = pemel.N_clt_stk*pemel.Ac_clt;   % [m^2]
pemel.totPrm_clt_stk = pemel.N_clt_stk*pemel.Prm_clt;  % [m]
pemel.totVol_clt_stk = pemel.N_clt_stk*pemel.Vol_clt;  % Volume [m^3]
pemel.totAs_clt_stk  = pemel.N_clt_stk*pemel.As_clt;   % Surface area [m^2]
pemel.totAc_ch_stk   = pemel.N_ch_stk*pemel.Ac_ch; % Total anode/cathode channel area [m^2]

% Params for total no. cells overall 
pemel.totN_cel   = pemel.N_cel*pemel.N_stk;	    % Total no. cells overall
pemel.totN_clt   = pemel.totN_cel*pemel.N_clt;  % Total no. cool tubes overall
pemel.totN_ch    = pemel.totN_cel*pemel.N_ch;   % Total no. proc. channels overall
pemel.totAc_clt  = pemel.totN_clt*pemel.Ac_clt;
pemel.totPrm_clt = pemel.totN_clt*pemel.Prm_clt;
pemel.totVol_clt = pemel.totN_clt*pemel.Vol_clt;
pemel.totAs_clt  = pemel.totN_clt*pemel.As_clt;
pemel.totAc_ch   = pemel.totN_ch*pemel.Ac_ch;

% Process fluid mass flow rates
h2o.mdot_reac_tot = pemel.totN_cel*const.M_h2o*input.I/(2*const.F);  % Total h2o consumed in reaction [kg/s]
h2o.mdot_in_tot   = h2o.stoich*h2o.mdot_reac_tot;					 % Total h2o inlet mass flow rate [kg/s]
h2o.mdot_out_tot  = h2o.mdot_in_tot - h2o.mdot_reac_tot;			 % Total h2o outlet mass flow rate [kg/s]
h2.mdot_reac_tot  = pemel.totN_cel*const.M_h2*input.I/(2*const.F);   % Total h2 mass produced [kg/s]


% Stack-level mass flow conversion
%pemel.H2out_stk_i = pemel.N_cel*pemel.H2out_i;  % H2 out [kg/s]

%%%%  DO NOT PUT INPUT PARAMETERS HERE! (put them above this section)  %%%%