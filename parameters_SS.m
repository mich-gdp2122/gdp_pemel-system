%% PLACE INPUT PARAMETERS IN THIS FILE %%

%%% Load constants and gas properties
const_file;
data_gas_file;

%% Ambient conditions
amb.T   = 20 + 273.15;       % Ambient temperature  [degC -> K]
amb.p   = 1  * 101325;       % Ambient air pressure [atm -> Pa]
amb.phi = 0.60;              % Ambient humidity     [0-1]

%% Operating conditions
% Stack
pemel.T_stk = 353;			% Nominal stack temperature [K]
pemel.p_ca  = 20 * 100000;  % Nominal cathode pressure [bar -> Pa]
pemel.p_an  = amb.p;        % Nominal anode pressure   [Pa]

% Process water
h2o.mdot_stk = 12;			% Nominal mass flow rate per stack [kg/s]
h2o.T_stk_in  = 345.7;		% Stack inlet temperature [K]

% Coolant
clnt.T_stk_in = 345;		% Coolant inlet temperature @ stack [K] 
clnt.dT_stk   = 5;          % Coolant T-diff thru stack     [degC or K]  (Default: 5)

% Cell V-i correlation
pemel.V_i = [1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2];        % [V]
pemel.i_i    = [0.010209, 0.062283, 0.24657, 0.58725, ...  % [A/cm^2]
				1.0239,   1.5091,   2.014,   2.5151];


%% Stack specifications
pemel.N_cel     = 30;          % Number of cells in stack []
pemel.A_cel     = 0.1*100^2;   % Cell active area         [cm^2]

pemel.A_mem   = pemel.A_cel;   % Membrane area           [m^2]
pemel.thk_mem = 1E-4;          % Membrane thickness      [m]
%pemel.M_mem   = 1.1;          % Membrane dry molar mass [kg/mol]
pemel.rho_mem = 1.98E-3;       % Membrane dry density    [kg/cm^3]

pemel.N_clt   = 50;            % No. cooling tubes per BP plate       []
pemel.Ac_clt  = 1E-6;          % Cooling tube cross-section area      [m^2]
pemel.Prm_clt = 4E-3;          % Cooling tube cross-section perimeter [m]
pemel.L_clt   = 0.5;           % Cooling tube length                  [m]

pemel.thk_bp  = 0.002;         % BP plate thickness                   [m]
pemel.cp_bp   = 100;           % BP plate Sp. heat (Cp)               [J/(kg*K)]
pemel.k_bp    = 1;             % BP plate heat conductivity           [W/(m*K)]

%pemel.Ac_an  = 2.5E-4;         % Anode channel cross-section area      [m^2]
%pemel.Prm_an = 0.002;          % Anode channel cross-section perimeter [m]
%pemel.L_an   = 0.01;           % Anode channel length                  [m]
%pemel.vol_an = 0.1;            % Cathode channel volume                [m^3]

%pemel.Ac_ca  = 2.5E-4;         % Cathode channel cross-section area      [m^2]
%pemel.Prm_ca = 0.002;          % Cathode channel cross-section perimeter [m]
%pemel.L_ca   = 0.01;           % Cathode channel length                  [m]
%pemel.vol_ca = 0.1;            % Cathode channel volume                  [m^3]

%% TMS block parameters
%%% Cooling water
Pump.CS  = 0.01;  %Pump inlet & outlet cross setion
TL.FlowRate  = 0.5;  %Thermal liquid flow rate Kg/s

%%% Electrolyser heat transfer
Pipe.length  = 1;  %Pipe lenght in m
Heat.load  = 26.875;  %Heat transfer from electrolyser into the cooling water in KW

%% BoP efficiencies
BoP.eff_pmp  = 0.90;  % Coolant pump efficiency []
BoP.eff_fan  = 0.90;  % Fan efficiency []

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Derived Parameters ***DO NOT MODIFY!!***
%%%%  DO NOT PUT INPUT PARAMETERS HERE! (put them above this section)  %%%%

pemel.N_clt_stk = pemel.N_clt*pemel.N_cel;    % No. cooling tubes in stack
pemel.Dh_clt = 4*pemel.Ac_clt/pemel.Prm_clt;  % Cooling tube hydraulic diameter [m]

pemel.I_i = pemel.i_i*pemel.A_cel;      % i -> I [A/cm^2 -> A/m^2 -> A]
pemel.V_stk_i = pemel.V_i*pemel.N_cel;  % Cell -> stack voltage [V]

%%%%  DO NOT PUT INPUT PARAMETERS HERE! (put them above this section)  %%%%