% FEEG6013 Group Design Project, 2021-2022
% Group 19
%
% Created by Michael and Victor
%
%
% PLACE SERPENTINE INPUT PARAMETERS IN THIS FILE
%
%
%%
% Load general parameters (used for both channel simulations)
parameters;

%% Cooling channel specifications
clch.N   = 2*50;			   % No. cooling tubes per BP plate       []
clch.L   = 0.598;			   % Cooling tube length                  [m]
clch.W   = 0.001;			   % Cooling tube width					  [m]
clch.H   = 0.0005;			   % Cooling tube height				  [m]
clch.Ac  = clch.W*clch.H;	   % Cooling tube cross-section area      [m^2]
clch.Prm = 2*(clch.W+clch.H);  % Cooling tube cross-section perimeter [m]


%% Coolant operating conditions
clnt.Vldot_ch  = 0.33*1e-6;	   % Volume flow rate per channel     [mL/s -> m^3/s]
clnt.T_stk_in  = 341.44;	   % Coolant stack inlet temperature  [K] 
clnt.T_stk_out = 347.32;	   % Coolant stack outlet temperature [K]
clnt.dp_stk    = 14209;		   % Coolant stack pressure drop	  [Pa]


%% Balance-of-Plant
% Correction factors
cf.L_ph   = 1.4;  % Correction factor for preheater length
cf.L_rjORC = 1;	  % Correction factor for heat rejector length


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%    DO NOT PUT INPUT PARAMETERS BELOW THIS LINE!    %%%%%%%%%%% %%

% Get derived parameters
parameters_derived;