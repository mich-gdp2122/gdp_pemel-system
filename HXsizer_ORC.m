function [L_c, L_h, Rt] = HXsizer_ORC(mdot_h, mdot_c, D_c, D_h, Th_in, Th_out, Tc, x3)
% Determines required length of hot & cold sides for preheater

%% 1) Calc thermal properties
[~, cp_h, k_h, mu_h, Pr_h] = data_water(Th_in, Th_out);  % Hot-side  (coolant)
[~,~,  ~, k_c, mu_c]	   = data_r600a_sat(Tc, 0, x3);  % Cold-side (r600a)

%% 2) Calc ovrl heat transf coeff, U
% Calc hot & cold heat transf coeff's
h_h = calc_h_pipeL( 'Ts', k_h, mu_h, Pr_h, D_h, mdot_h);
h_c = calc_h_pipe2P('Qs', k_c, mu_c, D_c, mdot_c, Tc, 0, x3);

% Overall heat trasnf coeff.
U = 1/((1/h_c) + (1/h_h));

%% 3) Calc HX surface area and lengths required
As = (mdot_h*cp_h/U)*log( (Th_in - Tc)/(Th_out - Tc) );

% Each-side pipe length required [m]
L_c = As/(pi*D_c);
L_h = As/(pi*D_h);

% Thermal resistance [K/W]
Rt = 1/(U*As);
end