function [L_h, L_c, Rt, As, U] = HXsizer_ORC(mdot_h, mdot_c, D_c, D_h, Th_in, Th_out, Tpp, Tc2f, Tc2, x3)
% % FEEG6013 Group Design Project, 2021-2022
% Group 19
%
% Created by Michael
%
%
% Determines performance and required surface area of ORC heat recovery unit
%
%
%% 1) Subcooled (boiler) region
% Calc thermal properties
data_hPH = data_water(Tpp, Th_out);  % Hot-side  (coolant)
data_c2f = data_r600a_sat(Tc2f, 0);   % Cool-side (r600a, on satL point)
data_c2  = data_r600a_sat(Tc2, 0);	  % Cool-side (r600a, state 2)

% Get mean subcooled properties
data_cPH.cp = mean([data_c2.cp, data_c2f.cp]);
data_cPH.k  = mean([data_c2.k, data_c2f.k]);
data_cPH.mu = mean([data_c2.mu, data_c2f.mu]);
data_cPH.Pr = mean([data_c2.Pr, data_c2f.Pr]);

% Calc heat capacity rates
C_hPH   = mdot_h*data_hPH.cp;
C_cPH   = mdot_c*data_cPH.cp;
C_minPH = min(C_cPH,C_hPH);
C_maxPH = max(C_cPH,C_hPH);
C_rtoPH = C_minPH/C_maxPH;

% Calc effectiveness and NTU
efct = ( C_hPH*(Tpp - Th_out) )/( C_minPH*(Tpp - Tc2) );		% Effectiveness
NTU  = (1/(1 - C_rtoPH))*log( (1 - efct*C_rtoPH)/(1 - efct) );  % NTU for counter-flow HX

% Calc hot & cold side heat transf coeff's
h_hPH = calc_h_pipeL('Qs', data_hPH.k, data_hPH.mu, data_hPH.Pr, D_h, mdot_h);
h_cPH = calc_h_pipeL('Qs', data_cPH.k, data_cPH.mu, data_cPH.Pr, D_c, mdot_c);
% Overall heat transf coeff.
U_PH  = 1/((1/h_cPH) + (1/h_hPH));

% Heat transf surface area, and hence each pipe length
AsPH = NTU*C_minPH/U_PH;


%% 2) Evaporation region
% Calc thermal properties
data_hEvp = data_water(Th_in, Tpp);       % Hot-side  (coolant)
data_cEvp = data_r600a_sat(Tc2f, 0, x3);  % Cool-side (r600a, states 2f->3)

% Calc hot & cold heat transf coeff's
h_hEvp = calc_h_pipeL( 'Ts', data_hEvp.k, data_hEvp.mu, data_hEvp.Pr, D_h, mdot_h);
h_cEvp = calc_h_pipe2P('Qs', data_cEvp.k, data_cEvp.mu, D_c, mdot_c, Tc2f, 0, x3);
% Overall heat trasnf coeff.
U_Evp = 1/((1/h_cEvp) + (1/h_hEvp));

% Calc HX surface area and lengths required
AsEvp = (mdot_h*data_hEvp.cp/U_Evp)*log( (Th_in - Tc2f)/(Tpp - Tc2f) );


%% 3) Combine sections to one HX
% Total surface area [m^2]
As = AsPH + AsEvp;

% Req'd hot-side pipe length [m] & pipe volume occupied [m^3]
L_h   = As/(pi*D_h);
L_c   = As/(pi*D_c);

% Average heat transf coeff [W/(m^2*K)]
U = (AsPH*U_PH + AsEvp*U_Evp)/As;

% Thermal resistance [K/W]
Rt = 1/(U*As);
end