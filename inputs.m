%% PLACE INPUT PARAMETERS IN THIS FILE %%
clearvars;                  % Clear all currently defined workspace variables

in  = struct();        % Create workspace structure for input parameters
ref = struct();		   % Reference values (for control)
el	= struct();		   % Electrolyser properties

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Output file properties %%

create_file    = false;     % Set to true to store below input parameters into .mat file
dont_overwrite = false;     % Set to true to create a new input file instead of overwriting existing one

filename = 'params';        % Output .mat file name (if create_file = true)
							% DON'T FORGET TO UPDATE FILENAME IN SIMLNK MODEL MANAGER

%% Parameters  %%

% Ambient conditions
in.T_amb     = 20;       % Ambient temperature [deg C]
in.p_amb     = 1;        % Ambient air pressure [atm]
in.phi_amb   = 0.60;     % Ambient humidity [0-1]

% Stack operating conditions
ref.T_stk     = 80;       % Nominal operating temperature [deg C]			  (Default: 80)
ref.dT_stk    = 5;        % Coolant T-diff between outlet & inlet [C or K]  (Default: 5)

ref.p_ca= 1.2;		% Reactant supply pressure [atm]		   (Default: 1.2)
in.dp_stk    = 0.2;		% Reactant pressure drop in channel [atm]  (Default: 0.2)


% Stack specifications
in.A_fc      = 200;       % Cell active area [cm^2]

in.M_mem     = 1.1;      % Membrane dry molar mass [kg/mol]
in.rho_mem   = 1.98E-3;  % Membrane dry density [kg/cm^3]
in.L_mem     = 0.01;     % Membrane thickness [cm]

in.D_cl_t    = 2;	    % Cooling tube diameter [mm]
in.N_t_bp    = 50;       % No. cooling tubes per bipolar plate []

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% BoP efficiencies
in.eff_cmp  = 0.70;      % Compressor efficiency []
in.eff_rcrc = 0.70;		% Recirculator efficiency []
in.eff_pmp  = 0.60;      % Coolant pump efficiency []
in.eff_fan  = 0.60;      % Fan efficiency []


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Conversions and Derived Parameters %%

in.T_amb     = in.T_amb + 273.15;		% Ambient Temperature [C -> K]
in.T_stk     = in.T_stk + 273.15;		% Stack temperature [C -> K]

in.N_fc		= ceil((in.Prat_stk/50)*360);  % No. cells (based on 50kW stk of 360 cells)
in.m_stk		= 30*(in.N_fc/360);			  % Stack mass [kg] (based on 360 cell stk of 30kg)
in.As_cl_stk = 5*(in.N_fc/360);			  % Stk cooling Q area (based on 360 cell stk w/ 5m^2 area) [m^2]

in.p_stk_out = in.p_stk - in.dp_stk;	% Stack outlet pressure [atm]
in.Af_veh    = in.W_veh*in.H_veh;		% Vehicle frontal area [m^2]


%% Save input file %%

if create_file == true
    if dont_overwrite == false
        % Overwrite any existing input files w/ specified filename
        filename = [filename,'.mat'];
        save(filename, 'input');
        disp([newline,'Input data saved as ''',filename,''' (previous file overwritten).']);
	else
        % Append current date & time to filename, creating new input file
        filename = [filename,'_',datestr(now,'yyyymmdd-HHMMSS'),'.mat'];    
        save(filename, 'input');
        disp([newline,'Input data saved as ''',filename,'''.']);
	end
end

disp([newline,'Input data successfully loaded into workspace.',newline,newline]);

clearvars -except input;   % Take out the trash