function [L_c, L_h, Rt, As, U, Th_out] = HXsizer_PH(mdot_c, mdot_h, D_c, D_h, Tc_in, Tc_out, Th_in, Th_out)
% FEEG6013 Group Design Project, 2021-2022
% Group 19
%
% Created by Michael
%
%
% Determines required length of hot & cold sides for preheater
%
%
%% 1) Calc thermal properties @ avg. temperatures
% Cool-side thermal properties, heat capacity rate
data_c = data_water(Tc_in, Tc_out);
C_c = mdot_c*data_c.cp;

% Iterate only if no Th_out value specified
if exist('Th_out','var') == 0
	Th_out = Th_in - 1;  % Initial Th_out guess
	while true
	% Iterate to find Th_out
		Th_out_prev = Th_out; % Hold prev value to determine convergence
	
		% Hot-side thermal properties, heat capacity rate
		data_h = data_water(Th_in, Th_out);
		C_h = mdot_h*data_h.cp;
		Th_out = Th_in - (C_c/C_h)*(Tc_out - Tc_in);  % Calc new Th_out
	
		% Compare new outlet T w/ prev. outlet T, then break loop if converged
		dTh_out = max(Th_out, Th_out_prev) - min(Th_out, Th_out_prev);
		if dTh_out < 0.0001
			break
		end
	end
	clear dTh_out Th_out_prev;
else
	% Hot-side thermal properties, heat capacity rate
	[~, data_h.cp, data_h.k, data_h.mu, data_h.Pr] = data_water(Th_in, Th_out);
	C_h = mdot_h*data_h.cp;
end

%% 3) Calc min & max heat capacity rates
C_min = min(C_c,C_h);
C_max = max(C_c,C_h);
C_rto = C_min/C_max;

%% 4) Calc effectiveness and NTU
efct = ( C_h*(Th_in - Th_out) )/( C_min*(Th_in - Tc_in) );  % Effectiveness
NTU  = (1/(1 - C_rto))*log( (1 - efct*C_rto)/(1 - efct) );  % NTU for counter-flow HX


%% 5) Calc ovrl heat transf coeff, U
% Calc hot & cold heat transf coeff's
h_c = calc_h_pipeL('Qs', data_c.k, data_c.mu, data_c.Pr, D_c, mdot_c);
h_h = calc_h_pipeL('Qs', data_h.k, data_h.mu, data_h.Pr, D_h, mdot_h);

% Overall heat transf coeff.
U = 1/((1/h_c) + (1/h_h));

%% 6) Calc heat transf surface area, and hence each pipe length
As = NTU*C_min/U;

% Each-side pipe length required [m]
L_c = As/(pi*D_c);
L_h = As/(pi*D_h);

% Thermal resistance [K/W]
Rt = 1/(U*As);
end