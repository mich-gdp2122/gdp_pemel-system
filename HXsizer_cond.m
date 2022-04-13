function [L, Rt] = HXsizer_cond(mdot, D, Tc, Th, x1, h_c)
% Determines required length of ORC condenser
%
% Cold side assumed stationary ocean w/ infinite thermal capacity
% (C_c = 0; dT_c = 0)

%% 1) Calc thermal properties and heat transf coeff's
% Hot side properties & heat transf coeff [W/(m^2*K)]
[~,~,~, k_h, mu_h,  ~, dh_41] = data_r600a_sat(Th, 1, x1);
h_h = calc_h_pipe2P('Ts', k_h, mu_h, D, mdot, Th, 1, x1);

% Cold-side properties & h_c via nat. convect. correlation (if no h_c value
% provided)
if exist('h_c','var') == 0
	[rho_c, ~, k_c, mu_c, Pr_c, beta_c] = data_water(Tc, Tc);
	h_c = calc_h_sink(rho_c, k_c, mu_c, Pr_c, beta_c, D, Th, Th, Tc);
end

% Overall heat trasnf coeff. [W/(m^2*K)]
U = 1/((1/h_c) + (1/h_h));


%% 2) Calc pipe length required and thermal resistance
% Length [m] & surface area m[^2], 
% via integrated 1st Law && Newton's Cooling Law relation
L  = ( mdot*abs(dh_41) )/( U*pi*D*(Th - Tc) );
As = pi*D*L;

% Thermal resistance [K/W]
Rt = 1/(U*As);
end


