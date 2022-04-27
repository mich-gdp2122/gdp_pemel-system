%% PLACE STRAIGHT CHANNEL INPUT PARAMETERS IN THIS FILE %%

%%% Load general parameters (used for both channel simulations)
parameters;

%% Cooling channel specifications
clch.N   = 2*60;				% No. cooling tubes per BP plate       []
clch.L   = 0.25;				% Cooling tube length                  [m]
clch.W   = 0.001;				% Cooling tube width				   [m]
clch.H   = 0.001;				% Cooling tube height				   [m]
clch.Ac  = clch.W*clch.H;		% Cooling tube cross-section area      [m^2]
clch.Prm = 2*(clch.W+clch.H);   % Cooling tube cross-section perimeter [m]

% % Bipolar plate
% bp.L   = 0.5;			% BP plate length [m]
% bp.W   = 0.2;			% BP plate width [m]
% bp.thk = 0.002;			% BP plate thickness [m]
% bp.rho = 4494;			% BP density [kg/m^3]
% bp.cp  = 534;			% BP plate Sp. heat (Cp)               [J/(kg*K)]
% bp.k   = 20.2;			% BP plate heat conductivity           [W/(m*K)]


%% Coolant operating conditions
clnt.Vldot_ch  = 0.25*1e-6;					   % Volume flow rate per channel [mL/s -> m^3/s]
clnt.T_stk_in  = 339;						   % Coolant inlet temperature @ stack [K] 
clnt.T_stk_out = clnt.T_stk_in + 6;			   % Coolant outlet temperature [K]
%clnt.p_stk_in  = 2742.5 + amb.p;			   % Coolant stack inlet pressure [Pa]
clnt.dp_stk    = 2600;						   % Coolant stack pressure drop [Pa]
%clnt.p_stk_out = clnt.p_stk_in - clnt.dp_stk;  % Coolant outlet pressure [Pa]


%% Balance-of-Plant
% Correction factors
cf.L_ph   = 2;  % Correction factor for preheater length
cf.L_phCL = 1.5;
cf.L_rjORC = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%    DO NOT PUT INPUT PARAMETERS BELOW THIS LINE!    %%%%%%%%%%% %%

% Get derived parameters
parameters_derived;