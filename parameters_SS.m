%% PLACE INPUT PARAMETERS IN THIS FILE %%

%%% Load constants and gas properties
const_file;
data_gas_file;

%% Ambient conditions
amb.T   = 20 + 273.15;       % Ambient temperature  [degC -> K]
amb.p   = 1  * 101325;       % Ambient air pressure [atm -> Pa]
amb.phi = 0.60;              % Ambient humidity     [0-1]

%% Stack specifications
pemel.N_stk = 1;		   % No. stacks				  []
pemel.N_cel = 30;          % Number of cells in stack []
pemel.A_cel = 0.1*100^2;   % Cell active area         [cm^2]

% Membrane
pemel.A_mem   = pemel.A_cel;   % Membrane area           [m^2]
pemel.thk_mem = 1E-4;          % Membrane thickness      [m]
%pemel.M_mem   = 1.1;          % Membrane dry molar mass [kg/mol]
pemel.rho_mem = 1.98E-3;       % Membrane dry density    [kg/cm^3]

% Cooling channels
pemel.N_clt   = 50;            % No. cooling tubes per BP plate       []
pemel.Ac_clt  = 1E-6;          % Cooling tube cross-section area      [m^2]
pemel.Prm_clt = 4E-3;          % Cooling tube cross-section perimeter [m]
pemel.L_clt   = 0.2;           % Cooling tube length                  [m]

% Bipolar plate
pemel.thk_bp  = 0.002;         % BP plate thickness                   [m]
pemel.cp_bp   = 100;           % BP plate Sp. heat (Cp)               [J/(kg*K)]
pemel.k_bp    = 1;             % BP plate heat conductivity           [W/(m*K)]

%pemel.Ac_an  = 2.5E-4;         % Anode channel cross-section area      [m^2]
%pemel.Prm_an = 0.002;          % Anode channel cross-section perimeter [m]
%pemel.L_an   = 0.01;           % Anode channel length                  [m]
%pemel.vol_an = 0.1;            % Cathode channel volume                [m^3]

pemel.Ac_ca  = 2.5E-4;          % Cathode channel cross-section area      [m^2]
%pemel.Prm_ca = 0.002;          % Cathode channel cross-section perimeter [m]
%pemel.L_ca   = 0.01;           % Cathode channel length                  [m]
%pemel.vol_ca = 0.1;            % Cathode channel volume                  [m^3]

% Cell V-i correlation, Temp vect, H2 mass flow rate out
pemel.V_i     =	[1.5,1.6,1.7,1.8,1.9,2,2.1,2.2,2.3,2.4,2.5];  % [V]
pemel.i_i     =	[0.14154,0.39234,0.76815,1.2179,1.7085,2.2219,2.7477,3.2785,3.8078,4.327,4.8181];  % [A/cm^2]
pemel.T_i     =	[340.11,340.37,340.84,341.47,342.2,343.03,343.93,344.9,345.94,347.04,348.21];  % [K]
pemel.H2out_i = [1.44E-06,4.00E-06,7.84E-06,1.25E-05,1.75E-05,2.28E-05,2.82E-05,3.36E-05,3.91E-05,4.45E-05,4.96E-05];  % [kg/s]

%clnt.p_stk_in = 

%% Operating conditions
% Stack
pemel.T_stk = 353.15; %pemel.T_i(6);			% Nominal stack temperature [K]
pemel.p_ca  = 20 * 100000;			% Nominal cathode pressure [bar -> Pa]
pemel.p_an  = amb.p;				% Nominal anode pressure   [Pa]

% Electrical
ctrl.i = pemel.i_i(6);			% Nominal current density [A/cm^2]

% Process water
h2o.mdot_stk  = 12;					% Nominal mass flow rate per stack [kg/s]
h2o.T_stk_in  = 345.7;				% Stack inlet temperature [K]

% Coolant
clnt.T_stk_in = 345;				% Coolant inlet temperature @ stack [K] 
clnt.dT_stk   = 5;					% Coolant T-diff thru stack     [degC or K]  (Default: 5)
clnt.mdot_stk = 12;					% Nominal mass flow rate per stack [kg/s]


%% BoP efficiencies
BoP.eff_pmp  = 0.90;  % Coolant pump efficiency []
BoP.eff_fan  = 0.90;  % Fan efficiency []


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Derived Parameters
%%%%  DO NOT PUT INPUT PARAMETERS HERE! (put them above this section)  %%%%

pemel.N_clt_stk = pemel.N_clt*pemel.N_cel;    % No. cooling tubes in stack

% Total cooling channel cross-section area & perimeter for single pipe
pemel.totAc_clt_stk  = pemel.N_clt_stk*pemel.Ac_clt;   % [m^2]
pemel.totPrm_clt_stk = pemel.N_clt_stk*pemel.Prm_clt;  % [m]
pemel.Dh_clt = 4*pemel.totAc_clt_stk/pemel.totPrm_clt_stk;  % Cooling tube hydraulic diameter [m]

% Electrical conversion
pemel.I_i     = pemel.i_i*pemel.A_cel;  % i -> I [A/cm^2 -> A/m^2 -> A]
pemel.V_stk_i = pemel.V_i*pemel.N_cel;  % Cell -> stack voltage [V]
ctrl.I        = ctrl.i*pemel.A_cel;     % i -> I [A/cm^2 -> A/m^2 -> A]

% Stack-level mass flow conversion
pemel.H2out_stk_i = pemel.N_cel*pemel.H2out_i;  % H2 out [kg/s]

%%%%  DO NOT PUT INPUT PARAMETERS HERE! (put them above this section)  %%%%