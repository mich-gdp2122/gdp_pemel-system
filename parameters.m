%% PLACE INPUT PARAMETERS IN THIS FILE %%

%%% Constants
const_file;

%%% Gas properties
data_gas_file;


% Ambient conditions
T_amb     = 20;       % Ambient temperature [deg C]
p_amb     = 1;        % Ambient air pressure [atm]
phi_amb   = 0.60;     % Ambient humidity [0-1]

% Stack operating conditions
T_stk     = 80;       % Nominal operating temperature [deg C]			  (Default: 80)
dT_stk    = 5;        % Coolant T-diff between outlet & inlet [C or K]  (Default: 5)

% Stack specifications
A_cel     = 200;      % Cell active area [cm^2]
M_mem     = 1.1;      % Membrane dry molar mass [kg/mol]
rho_mem   = 1.98E-3;  % Membrane dry density [kg/cm^3]
L_mem     = 0.01;     % Membrane thickness [cm]
D_cl_t    = 2;	      % Cooling tube diameter [mm]
N_t_bp    = 50;       % No. cooling tubes per bipolar plate []

% BoP efficiencies
eff_cmp  = 0.90;      % Compressor efficiency []
eff_rcrc = 0.90;	  % Recirculator efficiency []
eff_pmp  = 0.90;      % Coolant pump efficiency []
eff_fan  = 0.90;      % Fan efficiency []