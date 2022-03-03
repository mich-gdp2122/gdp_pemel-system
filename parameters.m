%% PLACE INPUT PARAMETERS IN THIS FILE %%

%%% Load constants and gas properties
const_file;
data_gas_file;

%% Ambient conditions
amb.T     = 20 + 273.15;       % Ambient temperature  [degC -> K]
amb.p     = 1  * 101325;       % Ambient air pressure [atm -> Pa]
amb.phi   = 0.60;              % Ambient humidity     [0-1]

%% Control inputs
control.T_stk  = 80 + 273.15;  % Nominal operating temperature [deg C]  (Default: 80)
control.dT_stk = 5;            % Coolant T-diff thru stack     [C or K] (Default: 5)


%% Stack specifications
pemel.N_cel     = 30;          % Number of cells in stack []
pemel.A_cel     = 0.1;         % Cell active area         [m^2]

%pemel.Ac_an  = 2.5E-4;         % Anode channel cross-section area      [m^2]
%pemel.Prm_an = 0.002;          % Anode channel cross-section perimeter [m]
%pemel.L_an   = 0.01;           % Anode channel length                  [m]
pemel.vol_an = 0.1;            % Cathode channel volume                [m^3]

%pemel.Ac_ca  = 2.5E-4;         % Cathode channel cross-section area      [m^2]
%pemel.Prm_ca = 0.002;          % Cathode channel cross-section perimeter [m]
%pemel.L_ca   = 0.01;           % Cathode channel length                  [m]
pemel.vol_ca = 0.1;            % Cathode channel volume                  [m^3]

pemel.A_mem   = pemel.A_cel;   % Membrane area           [m^2]
pemel.thk_mem = 1E-4;          % Membrane thickness      [m]
%pemel.M_mem   = 1.1;           % Membrane dry molar mass [kg/mol]
%pemel.rho     = 1.98E-3;       % Membrane dry density    [kg/cm^3]

pemel.N_clt   = 50;             % No. cooling tubes per BP plate      []
pemel.Ac_clt  = 1E-6;          % Cooling tube cross-section area      [m^2]
pemel.Prm_clt = 4E-3;          % Cooling tube cross-section perimeter [m]
pemel.L_clt   = 0.5;           % Cooling tube length                  [m]

pemel.thk_bp  = 0.002;         % BP plate thickness                   [m]
pemel.cp_bp   = 100;           % BP plate Sp. heat (Cp)               [J/(kg*K)]
pemel.k_bp    = 1;             % BP plate heat conductivity           [W/(m*K)]


%% BoP efficiencies
BoP.eff_pmp  = 0.90;  % Coolant pump efficiency []
BoP.eff_fan  = 0.90;  % Fan efficiency []