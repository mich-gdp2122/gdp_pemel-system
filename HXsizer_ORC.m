function [L_h, L_c, Rt, As] = HXsizer_ORC(cf_L, mdot_h, mdot_c, D_c, D_h, ...
									Th_in, Th_out, Tc, x3)
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

% Req'd hot-side pipe length [m] & pipe volume occupied [m^3]
L_h   = cf_L*As/(pi*D_h);
%Vol_c = (D_c/4)*As;
L_c   = cf_L*As/(pi*D_c);
%Vol_h = (D_h/4)*As;
%Vol_c = VolRto*Vol_h;

% Thermal resistance [K/W]
Rt = 1/(U*As);
end