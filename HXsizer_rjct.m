function [L, Rt, As, U] = HXsizer_rjct(mdot_h, D_h, Tc, Th_in, Th_out, h_c)
% Determines required length of coolant heat rejector
%
% Cold side assumed stationary ocean w/ infinite thermal capacity
% (C_c = 0; dT_c = 0)

%% 1) Calc thermal properties and heat transf coeff's
% Hot side properties & heat transf coeff [W/(m^2*K)]
data_h = data_water(Th_in, Th_out);
h_h = calc_h_pipeL('Ts', data_h.k, data_h.mu, data_h.Pr, D_h, mdot_h);  % Hot side

% Cold-side properties & h_c via nat. convect. correlation (if no h_c value
% provided)
if exist('h_c','var') == 0
	data_c = data_water(Tc, Tc);
	h_c = calc_h_sink(data_c.rho, data_c.k, data_c.mu, data_c.Pr, data_c.beta, ...
		D_h, Th_in, Th_out, Tc);
end
% Overall heat trasnf coeff. [W/(m^2*K)]
U = 1/((1/h_c) + (1/h_h));

%% 2) Calc pipe length required and thermal resistance
% Surface area m[^2] and length [m^2]
% via integrated 1st Law && Newton's Cooling Law relation
As = ((mdot_h*data_h.cp)/U)*log((Th_in - Tc)/(Th_out - Tc));
L  = As/(pi*D_h);

% Thermal resistance [K/W]
Rt = 1/(U*As);
end