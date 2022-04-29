function [L_h, L_c, Rt, As, U] = ...
	HXsizer_ORC(mdot_h, mdot_c, D_c, D_h, Th_in, Th_out, Tpp, Tc2f, x3)
% Determines required length of hot & cold sides for ORC recovery heater
%% 1) Subcooled (preheat) region
% Calc thermal properties
data_h2f3 = data_water(Tpp, Th_out);      % Hot-side  (coolant)
data_c2f3 = data_r600a_sat(Tc2f, 0, x3);  % Cool-side (r600a, state 2)



%% 2) Evaporation region
% Calc thermal properties
data_h2f3 = data_water(Th_in, Tpp);       % Hot-side  (coolant)
data_c2f3 = data_r600a_sat(Tc2f, 0, x3);  % Cool-side (r600a)

% Calc hot & cold heat transf coeff's
h_h = calc_h_pipeL( 'Ts', data_h2f3.k, data_h2f3.mu, data_h2f3.Pr, D_h, mdot_h);
h_c = calc_h_pipe2P('Qs', data_c2f3.k, data_c2f3.mu, D_c, mdot_c, Tc2f, 0, x3);
% Overall heat trasnf coeff.
U = 1/((1/h_c) + (1/h_h));

% Calc HX surface area and lengths required
As = (mdot_h*data_h2f3.cp/U)*log( (Th_in - Tc2f)/(Th_out - Tc2f) );

% Req'd hot-side pipe length [m] & pipe volume occupied [m^3]
L_h   = As/(pi*D_h);
%Vol_c = (D_c/4)*As;
L_c   = As/(pi*D_c);
%Vol_h = (D_h/4)*As;
%Vol_c = VolRto*Vol_h;

% Thermal resistance [K/W]
Rt = 1/(U*As);
end