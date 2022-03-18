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

pemel.Ac_ch  = 2.5E-4;          % Process channel cross-section area      [m^2]
%pemel.Prm_ch = 0.002;          % Process channel cross-section perimeter [m]
%pemel.L_ch   = 0.01;           % Process channel length                  [m]
%pemel.vol_ch = 0.1;            % Process channel volume                  [m^3]


% Cell V-i correlation, Temp vect, H2 mass flow rate out
pemel.V_i     =	[1.5,1.6,1.7,1.8,1.9,2,2.1,2.2,2.3,2.4,2.5];  % [V]
pemel.i_i     =	[0.14154,0.39234,0.76815,1.2179,1.7085,2.2219,2.7477,3.2785,3.8078,4.327,4.8181];  % [A/cm^2]
pemel.T_i     =	[340.11,340.37,340.84,341.47,342.2,343.03,343.93,344.9,345.94,347.04,348.21];  % [K]

pemel.H2out_i = [1.44E-06,4.00E-06,7.84E-06,1.25E-05,1.75E-05,2.28E-05,2.82E-05,3.36E-05,3.91E-05,4.45E-05,4.96E-05];  % [kg/s]

pemel.pIn_clnt_i  = [101388.049,101387.917,101387.686,101387.382,101387.025,101386.629,101386.204,101385.746,101385.266,101384.758,101384.229];  % [Pa]
pemel.dp_clnt_i   = [62.83903,62.70745,62.47717,62.17412,61.81823,61.42347,60.99981,60.54321,60.06471,59.55826,59.03087];  % [Pa]
pemel.pOut_clnt_i = pemel.pIn_clnt_i - pemel.dp_clnt_i;  % [Pa]

pemel.TIn_clnt_i  = [340.01,340.02,340.05,340.09,340.14,340.19,340.25,340.31,340.38,340.44,340.52];  % [K]
pemel.dT_clnt_i   = [0.04,0.17,0.39,0.67,1,1.38,1.79,2.23,2.7,3.2,3.73];  % [K]
pemel.TOut_clnt_i = pemel.TIn_clnt_i + pemel.dT_clnt_i;  % [K]


%% Operating conditions
% Stack
pemel.p_ca  = 20 * 100000;  % Nominal cathode pressure [bar -> Pa]
pemel.p_an  = amb.p;	    % Nominal anode pressure   [Pa]

% Electrical
input.i = pemel.i_i(6);	    % Nominal current density [A/cm^2]

% Process water
h2o.mdot_stk  = 12;			% Nominal mass flow rate per stack [kg/s]
h2o.T_stk_in  = 345.7;		% Stack inlet temperature [K]

% Coolant
clnt.mdot_stk = 12;			% Nominal mass flow rate per stack [kg/s]


%% BoP efficiencies
BoP.eff_pmp  = 0.90;  % Coolant pump efficiency []
BoP.eff_fan  = 0.90;  % Fan efficiency []


%% Organic Rankine Cycle
%%%Condenser
Con.Fluid_V    = 5;      %Condenser fluid volume [m^3]
Con.PortA_A    = 0.01;   %Condenser port A area [m^2]
Con.PortB_A    = 0.01;   %Condenser port A area [m^2]

%Thermal 
Con.Area       = 1e-4;   %Condenser thermal resistance area [m^2]
Con.Heat_Coeff = 20;     %Condenser thermal resistance heat coefficient [W/(m^2*K)]
Con.Coolant_T  = 10;     %Condenser cooling temperature [K]


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Derived Parameters
%%%%  DO NOT PUT INPUT PARAMETERS HERE! (put them above this section)  %%%%

pemel.N_clt_stk = pemel.N_clt*pemel.N_cel;    % No. cooling tubes in stack

% Total cooling channel cross-section area & perimeter for single pipe
pemel.totAc_clt_stk  = pemel.N_clt_stk*pemel.Ac_clt;   % [m^2]
pemel.totPrm_clt_stk = pemel.N_clt_stk*pemel.Prm_clt;  % [m]
pemel.Dh_clt = 4*pemel.totAc_clt_stk/pemel.totPrm_clt_stk;  % Cooling tube hydraulic diameter [m]

% Total anode/cathode channel area
pemel.totAc_ch_stk = pemel.N_stk*pemel.Ac_c;  % [m^2]

% Electrical conversion
pemel.I_i     = pemel.i_i*pemel.A_cel;  % i -> I [A/cm^2 -> A/m^2 -> A]
pemel.V_stk_i = pemel.V_i*pemel.N_cel;  % Cell -> stack voltage [V]
input.I       = input.i*pemel.A_cel;     % i -> I [A/cm^2 -> A/m^2 -> A]

% Stack-level mass flow conversion
pemel.H2out_stk_i = pemel.N_cel*pemel.H2out_i;  % H2 out [kg/s]

% Nominal stack operating conditions
pemel.T_stk   = interp1(pemel.i_i, pemel.T_i, input.i, 'makima', 'extrap');  % Nominal stack temperature [K]
clnt.T_stk_in = interp1(pemel.i_i, pemel.TIn_clnt_i, input.i, 'makima', 'extrap');  % Coolant inlet temperature @ stack [K] 
clnt.dT_stk   = interp1(pemel.i_i, pemel.dT_clnt_i, input.i, 'makima', 'extrap');  % Coolant T-diff thru stack [K]
clnt.p_stk_in = interp1(pemel.i_i, pemel.pIn_clnt_i, input.i, 'makima', 'extrap'); % Coolant stack inlet pressure [Pa]
clnt.dp_stk   = interp1(pemel.i_i, pemel.dp_clnt_i, input.i, 'makima', 'extrap');  % Coolant stack pressure drop [Pa]

%%%%  DO NOT PUT INPUT PARAMETERS HERE! (put them above this section)  %%%%