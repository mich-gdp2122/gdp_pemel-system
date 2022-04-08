function [L_c, L_h] = ph_sizer(mdot_c,mdot_h,D_c,D_h,...
											Tc_in,Tc_out,Th_in,Th_out)
% Determines required length of hot & cold sides for preheater

%% 1) Calc thermal properties @ avg. temperatures
% Cool-side thermal properties, heat capacity rate
[rho_c, cp_c, k_c, mu_c, Pr_c] = data_water(Tc_in, Tc_out);
C_c = mdot_c*cp_c;

% Iterate only if no Th_out value specified
if exist('Th_out','var') == 0
	Th_out = Th_in - 1;  % Initial Th_out guess
	while true
	% Iterate to find Th_out
		Th_out_prev = Th_out; % Hold prev value to determine convergence
	
		% Hot-side thermal properties, heat capacity rate
		[rho_h, cp_h, k_h, mu_h, Pr_h] = data_water(Th_in, Th_out);
		C_h = mdot_h*cp_h;
		Th_out = Th_in - (C_c/C_h)*(Tc_out - Tc_in);  % Calc new Th_out
	
		% Compare new outlet T w/ prev. outlet T, then break loop if converged
		dTh_out = max(Th_out, Th_out_prev) - min(Th_out, Th_out_prev);
		if dTh_out < 0.01
			break
		end
	end
	clear dTh_out Th_out_prev;
end

%% 3) Calc min & max heat capacity rates
C_min = min(C_c,C_h);
C_max = max(C_c,C_h);
C_rat = C_min/C_max;

%% 4) Calc effectiveness and NTU
efct = ( C_h*(Th_in - Th_out) )/( C_min*(Th_in - Tc_in) );  % Effectiveness
NTU  = (1/(1 - C_rat))*log( (1 - efct*C_rat)/(1 - efct) );  % NTU for counter-flow HX


%% 5) Calc ovrl heat transf coeff, U
% Calc hot & cold heat transf coeff's
h_c = calc_h(rho_c, k_c, mu_c, Pr_c, D_c, mdot_c);
h_h = calc_h(rho_h, k_h, mu_h, Pr_h, D_h, mdot_h);

% Overall heat trasnf coeff.
U = 1/((1/h_c) + (1/h_h));

%% 6) Calc heat transf surface area, and hence each pipe length
As = NTU*C_min/U;

% Each-side pipe length required [m]
L_c = As/(pi*D_c);
L_h = As/(pi*D_h);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calc heat transfer coefficient
function h = calc_h(rho, k, mu, Pr, Dh, mdot)
	% Avg. velocity [m/s], Reynold's number
	u = mdot/(rho*(pi/4)*Dh);
	Re = (rho*u*Dh)/mu;

	% Reynold's no. transition region
	Re_lam = 2300;
	Re_tur = 3000;
	
	% Darcy friction factor []
	f_lam = 64/Re;							% Laminar case
	f_tur = ( 0.790*log(Re) - 1.64 )^(-2);  % Turbulent case (Petukov)

	if Re <= Re_lam						% Laminar case
		f = f_lam;
	elseif Re >= Re_tur					% Turbulent case
		f = f_tur;
	elseif Re > Re_lam && Re < Re_tur	% Transitional case (via interpolation)
		f = blend(f_lam, f_tur, Re_lam, Re_tur, Re);
	end
	
	% Coolant Nusselt number []
	Nu_lam= 4.36;							% Laminar case (const. heat flux)
	Nu_tur = ( (f/8)*(Re - 1000)*Pr )/...	% Turbulent case (Gnielinski)
			 ( 1 + 12.7*((f/8)^0.5)*( Pr^(2/3) - 1 ) );  

	if Re <= Re_lam						% Laminar case
		Nu = Nu_lam;
	elseif Re >= Re_tur					% Turbulent case
		Nu = Nu_tur;
	elseif Re > Re_lam && Re < Re_tur	% Transitional case (via interpolation)
		Nu = blend(Nu_lam, Nu_tur, Re_lam, Re_tur, Re);
	end
	h = k*Nu/Dh;  % Coolant heat transfer coefficient [J/(m^2 K)]
end

%% Blend function
function y = blend(y1,y2,x1,x2,x)
% Blend between two values y1 and y2 as x varies between x1 and x2. The
% blending function ensures y is continuous and differentiable in the
% range x1 <= x <= x2.

% Copyright 2017-2018 The MathWorks, Inc.

% [Source: simscape.function.blend.scc]
	u = (x-x1)/(x2-x1);
	transition = 3*u^2 - 2*u^3;
	if le(x,x1)
		y = y1;
	elseif ge(x,x2)
    	y = y2;
	else
    	y = (1-transition)*y1 + transition*y2;
	end
end