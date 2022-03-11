%%%% Place all constants here
%load('const.mat');

%%
const.R	     = 8.314;    % Universal gas const [J/(mol*K)]
const.F	     = 96485.3;  % Faraday  [C/mol]
const.HHV_h2 = 1.419E8;  % H2 higher heating value [J/kg]

% Density [kg/m^3]
const.rho_h2  = 0.082;   % H2
const.rho_o2  = 1.31;    % O2
const.rho_h2o = 1000;    % Water

% Molar mass [kg/mol]
const.M_h2  = 0.0020158;  % H2
const.M_o2  = 0.0319988;  % O2
const.M_h2o = 0.0180152;  % Water

% Sp. gas  [J/(kg*K)]
const.R_h2  = round(const.R/const.M_h2,  5, 'significant');  % H2
const.R_o2  = round(const.R/const.M_o2,  5, 'significant');  % O2
const.R_h2o = round(const.R/const.M_h2o, 5, 'significant');  % Water

% Membrane params
const.Dw_mem  = 1.28E-10;  % Water diffusion coeff. @ Nafion membrane []
const.K_darcy = 1.58E-18;  % Nafion membrane permeability to water [m^2]