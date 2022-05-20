function [L, Rt, As, U] = HXsizer_cond(mdot, Qout, D, Tc, Th, x1, h_c)
% FEEG6013 Group Design Project, 2021-2022
% Group 19
%
% Created by Michael
%
%
% Determines performance and required surface area of ORC condenser
%
% Cold side assumed stationary ocean w/ infinite thermal capacity
% (C_c = 0; dT_c = 0)
%
%
%% 1) Calc thermal properties and heat transf coeff's
% Hot side properties & heat transf coeff [W/(m^2*K)]
data_h = data_r600a_sat(Th, 1, x1);
h_h = calc_h_pipe2P('Ts', data_h.k, data_h.mu, D, mdot, Th, 1, x1);

% Cold-side properties & h_c via nat. convect. correlation (if no h_c value
% provided)
if exist('h_c','var') == 0
	data_c = data_water(Tc, Tc);
	h_c = calc_h_sink(data_c.rho, data_c.k, data_c.mu, data_c.Pr, data_c.beta, ...
		D, Th, Th, Tc);
end

% Overall heat trasnf coeff. [W/(m^2*K)]
U = 1/((1/h_c) + (1/h_h));


%% 2) Calc fluid volume required and thermal resistance
% Surface area [m^2] and pipe length [m], 
% via integrated 1st Law && Newton's Cooling Law relation
As = Qout/( U*(Th - Tc) );
L  = As/(pi*D);

% Thermal resistance [K/W]
Rt = 1/(U*As);
end


