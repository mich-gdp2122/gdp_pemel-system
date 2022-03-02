%%%% Place all constants here
load('const.mat');

%%
R	   = 8.314;    % Universal gas const [J/(mol*K)]
F	   = 96485.3;  % Faraday  [C/mol]
HHV_h2 = 1.419E8;  % H2 higher heating value [J/kg]

% Density [kg/m^3]
rho_h2  = 0.082;   % H2
rho_o2  = 1.31;    % O2
rho_h2o = 1000;    % Water

% Molar mass [kg/mol]
M_h2  = 0.0020158;  % H2
M_o2  = 0.0319988;  % O2
M_h2o = 0.0180152;  % Water

% Sp. gas  [J/(kg*K)]
R_h2  = round(R/M_h2,  5, 'significant');  % H2
R_o2  = round(R/M_o2,  5, 'significant');  % O2
R_h2o = round(R/M_h2o, 5, 'significant');  % Water

% Membrane params
