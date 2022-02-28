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
        out = R;            % Universal gas const. [J/(mol*K)]
        
    case 'p_atm'
        p_atm = 101325;     % 1 atm pressure [Pa]
        out   = p_atm;

    case 'f'
        F   = 96485.3;      % Faraday const. [C/mol]
        out = F;

    case 'hhv_h2'
        HHV_H2 = 1.419E8;   % H2 higher heating value [J/kg]
        out    = HHV_H2;

        
	%%%%%%  Density [kg/m^3]  %%%%%%
    case 'rho_h2'
        rho_H2  = 0.082;    % H2
        out     = rho_H2;

    case 'rho_o2'
        rho_O2  = 1.31;     % O2
        out     = rho_O2;

    case 'rho_h2o'
        rho_H2O = 1000;     % Water
        out     = rho_H2O;
      
    case 'rho_clnt'
        rho_clnt = 1000;    % Coolant (liquid water)
        out      = rho_clnt;

    
	%%%%%%  Molar mass [kg/mol]  %%%%%%
    case 'm_h2'
        out  = M_H2;        % H2
        
    case 'm_o2'
        out  = M_O2;        % O2

    case 'm_h2o'
        out   = M_H2O;      % Water (liquid & vapour)
        
        
    %%%%%%  Sp. gas const. [J/(kg*K)] %%%%%
    case 'r_h2'      
        out = R/M_H2;      % H2
    
    case 'r_o2'
        out = R/M_O2;      % O2
        
    case 'r_h2o'
        out = R/M_H2O;     % Water (vapour only)

 
	otherwise
		error(['Unknown const parameter, ', param]);
end
end