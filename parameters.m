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
pemel.N_stk   = 8;		     % No. stacks		  []
pemel.N_cel   = 30;          % No. cells in stack []
pemel.A_cel   = 0.1*100^2;   % Cell active area   [cm^2]
pemel.eff_stk = 0.79;		 % Stack efficiency   []

% Cooling channels
clch.N   = 60;            % No. cooling tubes per BP plate       []
clch.L   = 0.25;          % Cooling tube length                  [m]
clch.W   = 0.001;		  % Cooling tube width					 [m]
clch.Ac  = clch.W^2;      % Cooling tube cross-section area      [m^2]
clch.Prm = 4*clch.W;      % Cooling tube cross-section perimeter [m]

% Bipolar plate
bp.L   = 0.5;			% BP plate length [m]
bp.W   = 0.2;			% BP plate width [m]
bp.thk = 0.002;			% BP plate thickness [m]
bp.rho = 4494;			% BP density [kg/m^3]
bp.cp  = 534;			% BP plate Sp. heat (Cp)               [J/(kg*K)]
bp.k   = 20.2;			% BP plate heat conductivity           [W/(m*K)]

% Process channels
prch.N   = 4;			% No. anode/cathode channels each per cell []
prch.W   = 0.8E-3;		% Process channel width
prch.Ac  = prch.W^2;    % Process channel cross-section area       [m^2]
prch.Prm = 4*prch.W;    % Process channel cross-section perimeter  [m]
prch.L   = 0.025;       % Process channel length                   [m]


%% Stack operating conditions
% Inputs
pemel.i = 2.2219;	    % Nominal current density [A/cm^2]

% Stack
pemel.p_ca  = 20 * 100000;  % Nominal cathode pressure [bar -> Pa]
pemel.p_an  = amb.p;	    % Nominal anode pressure   [Pa]
pemel.T_stk = 80 + 273.15;  % Nominal stack temperature [K]
pemel.V_stk = 60;			% Nominal stack voltage [V]
%pemel.q     = 8955.5;		% Stack heat flux density [W/m^2]

% Process water
h2o.stoich    = 3.70;		% Stoichiometric ratio
h2o.T_stk_in  = 340;		% Stack inlet temperature [K]
h2o.T_stk_out = 341;		% Stack outlet temperature [K]
h2o.dp_stk    = 5;			% Stack pressure drop [bar]
%h2o.mdot_stk  = 12;		% Nominal mass flow rate per stack [kg/s]

% Coolant
clnt.mdot_stk  = 0.9;						   % Nominal mass flow rate per stack [kg/s]
clnt.T_stk_in  = 339.15;					   % Coolant inlet temperature @ stack [K] 
clnt.T_stk_out = 345.15;					   % Coolant outlet temperature [K]
clnt.p_stk_in  = 2742.5 + amb.p;			   % Coolant stack inlet pressure [Pa]
clnt.dp_stk    = 2681;						   % Coolant stack pressure drop [Pa]
clnt.p_stk_out = clnt.p_stk_in - clnt.dp_stk;  % Coolant outlet pressure [Pa]


%% Balance-of-Plant
% Efficiencies
BoP.eff_pmp  = 0.80;  % Pump efficiency    []
BoP.eff_tbne = 0.80;  % Turbine efficiency []

% Correction factors
cf.L_ph   = 1.4;  % Correction factor for preheater length
cf.L_phCL = 1.4;
cf.L_rjORC = 1;

% Rated powers
htr.PwrRt = 50;		 % Elec. heater rated power [kW]

%% Organic Rankine Cycle
% Operating conditions
%ORC.mdot = 1.7;				% Nominal mass flow rate [kg/s]
ORC.dTh  = 5;					% ORC HX driving temp diff [K]
ORC.dTc  = 10;					% ORC condensor driving temp [K]
%ORC.x1	 = 0.3;					% State 1 quality (vapour mass fraction)
%ORC.x3   = 0.95;				% State 3 quality


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%    DO NOT PUT INPUT PARAMETERS BELOW THIS LINE!    %%%%%%%%%%% %%

% Get derived parameters
parameters_derived;