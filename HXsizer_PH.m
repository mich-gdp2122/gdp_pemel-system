function [L_c, L_h, Rt, As, U, Th_out] = ...
	HXsizer_PH(mdot_c, mdot_h, D_c, D_h, Tc_in, Tc_out, Th_in, Th_out)
% Determines required length of hot & cold sides for preheater

%% 1) Calc thermal properties @ avg. temperatures
% Cool-side thermal properties, heat capacity rate
[~, cp_c, k_c, mu_c, Pr_c] = data_water(Tc_in, Tc_out);
C_c = mdot_c*cp_c;

% Iterate only if no Th_out value specified
if exist('Th_out','var') == 0
	Th_out = Th_in - 1;  % Initial Th_out guess
	while true
	% Iterate to find Th_out
		Th_out_prev = Th_out; % Hold prev value to determine convergence
	
		% Hot-side thermal properties, heat capacity rate
		[~, cp_h, k_h, mu_h, Pr_h] = data_water(Th_in, Th_out);
		C_h = mdot_h*cp_h;
		Th_out = Th_in - (C_c/C_h)*(Tc_out - Tc_in);  % Calc new Th_out
	
		% Compare new outlet T w/ prev. outlet T, then break loop if converged
		dTh_out = max(Th_out, Th_out_prev) - min(Th_out, Th_out_prev);
		if dTh_out < 0.01
			break
		end
	end
	clear dTh_out Th_out_prev;
else
	% Hot-side thermal properties, heat capacity rate
	[~, cp_h, k_h, mu_h, Pr_h] = data_water(Th_in, Th_out);
	C_h = mdot_h*cp_h;
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
h_c = calc_h_pipeL('Qs', k_c, mu_c, Pr_c, D_c, mdot_c);
h_h = calc_h_pipeL('Qs', k_h, mu_h, Pr_h, D_h, mdot_h);

% Overall heat trasnf coeff.
U = 1/((1/h_c) + (1/h_h));

%% 6) Calc heat transf surface area, and hence each pipe length
As = NTU*C_min/U;

% Each-side pipe length required [m]
L_c = As/(pi*D_c);
L_h = As/(pi*D_h);

% Thermal resistance [K/W]
Rt = 1/(U*As);
end