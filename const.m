function out = const(param)
% Place all constants in here

%%
param = lower(param);       % Force input value to lowercase

% Const's used in multiple cases
    R     = 8.314;          % Universal gas const [J/(mol*K)]
    M_H2  = 0.0020158;      % H2 molar mass [kg/mol]
    M_O2  = 0.0319988;      % O2 molar mass [kg/mol]
    M_H2O = 0.0180152;      % Water molar mass [kg/mol]

switch param
    case 'r'
        out = R;          % Universal gas const. [J/(mol*K)] 
    case 'p_atm'
        out = 101325;     % 1 atm pressure [Pa]
    case 'f'
        out   = 96485.3;  % Faraday const. [C/mol]
    case 'hhv_h2'
        out = 1.419E8;    % H2 higher heating value [J/kg]
        
	%%%%%%  Density [kg/m^3]  %%%%%%
    case 'rho_h2'
        out  = 0.082;   % H2
    case 'rho_o2'
        out  = 1.31;    % O2
    case 'rho_h2o'
        out = 1000;     % Water
    case 'rho_clnt'
        out = 1000;     % Coolant (liquid water)
    
	%%%%%%  Molar mass [kg/mol]  %%%%%%
    case 'm_h2'
        out  = M_H2;    % H2
    case 'm_o2'
        out  = M_O2;    % O2
    case 'm_h2o'
        out   = M_H2O;  % Water (liquid & vapour)
        
    %%%%%%  Sp. gas const. [J/(kg*K)]  %%%%%
    case 'r_h2'      
        out = R/M_H2;   % H2
    case 'r_o2'
        out = R/M_O2;   % O2
    case 'r_h2o'
        out = R/M_H2O;  % Water (vapour only)

	%%%%%  Membrane  %%%%%
	case 'dw_mem'
		out = 1.28E-10; % Water diffusion coeff. @ Nafion membrane []
	case 'k_darcy'
		out = 1.58E-18; % Nafion membrane permeability to water [m^2]
 

	otherwise
		error(['Unknown const parameter, ', param]);
end
end