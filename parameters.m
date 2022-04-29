%% PLACE GENERAL INPUT PARAMETERS IN THIS FILE (used for both channel simulations)

%%% Load constants and gas properties
const_file;
data_gas_file;

%% Ambient conditions
amb.T     = 20 + 273.15;       % Ambient temperature  [K]
amb.p     = 1  * 101325;       % Ambient air pressure [Pa]
amb.phi   = 0.60;              % Ambient humidity     [0-1]
amb.T_sea = 10 + 273.15;	   % Ocean background temperature [K]


%% General stack specifications
pemel.N_stk   = 8;		     % No. stacks		  []
pemel.N_cel   = 30;          % No. cells in stack []
pemel.A_cel   = 0.1*100^2;   % Cell active area   [cm^2]
pemel.eff_stk = 0.79;		 % Stack efficiency   []

% Process channels
prch.N   = 4;			% No. anode/cathode channels each per cell []
prch.W   = 0.8E-3;		% Process channel width
prch.Ac  = prch.W^2;    % Process channel cross-section area       [m^2]
prch.Prm = 4*prch.W;    % Process channel cross-section perimeter  [m]
prch.L   = 0.025;       % Process channel length                   [m]


%% General stack operating conditions
% Stack
pemel.i     = 2.2219;		% Nominal current density [A/cm^2]
pemel.V_stk = 60;			% Nominal stack voltage [V]
pemel.p_ca  = 20 * 100000;  % Nominal cathode pressure [bar -> Pa]
pemel.p_an  = amb.p;	    % Nominal anode pressure   [Pa]
pemel.T_stk = 80 + 273.15;  % Nominal stack temperature [K]

% Process water
h2o.stoich    = 3.70;		% Stoichiometric ratio
h2o.T_stk_in  = 340;		% Stack inlet temperature [K]
h2o.T_stk_out = 341;		% Stack outlet temperature [K]
h2o.dp_stk    = 5;			% Stack pressure drop [bar]


%% Balance-of-Plant
% Efficiencies
BoP.eff_pmp  = 0.85;  % Pump efficiency    []
BoP.eff_tbne = 0.85;  % Turbine efficiency []

% Rated powers
htr.PwrRt = 50;		 % Elec. heater rated power [kW]


%% Organic Rankine Cycle
% Operating conditions
%ORC.mdot = 1.7;				% Nominal mass flow rate [kg/s]
ORC.dTh  = 5;					% ORC HX driving temp diff [K]
ORC.dTc  = 5;					% ORC condensor driving temp [K]
%ORC.x1	 = 0.3;					% State 1 quality (vapour mass fraction)
%ORC.x3   = 0.95;				% State 3 quality
