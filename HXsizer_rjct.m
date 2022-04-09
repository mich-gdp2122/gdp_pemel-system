function [L_h, Rt, As_h, U] = HXsizer_rjct(mdot_h,D_h,Tc,Th_in,Th_out, h_c)
% Determines required length of coolant heat rejector
%
% Cold side assumed stationary ocean w/ infinite thermal capacity
% (C_c = 0; dT_c = 0)

%% 1) Calc thermal properties and heat transf coeff's
% Hot side properties & heat transf coeff [W/(m^2*K)]
[rho_h, cp_h, k_h, mu_h, Pr_h] = data_water(Th_in, Th_out);
h_h = calc_h_h(rho_h, k_h, mu_h, Pr_h, D_h, mdot_h);  % Hot side

% Cold-side properties & h_c via nat. convect. correlation (if no h_c value
% provided)
if exist('h_c','var') == 0
	[rho_c, ~, k_c, mu_c, Pr_c, beta_c] = data_water(Tc, Tc);
	h_c = calc_h_c(rho_c, k_c, mu_c, Pr_c, beta_c, D_h, Th_in, Th_out, Tc);
end

% Overall heat trasnf coeff. [W/(m^2*K)]
U = 1/((1/h_c) + (1/h_h));


%% 2) Calc pipe length required and thermal resistance
% Length [m] & surface area m[^2], 
% via integrated 1st Law && Newton's Cooling Law relation
L_h = ((mdot_h*cp_h)/(U*pi*D_h))*log((Th_in - Tc)/(Th_out - Tc));
As_h = pi*D_h*L_h;

% Thermal resistance [K/W]
Rt = 1/(U*As_h);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calc h_h coefficient
function h_h = calc_h_h(rho, k, mu, Pr, Dh, mdot)
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

	% Heat transfer coefficient [W/(m^2 K)]
	h_h = k*Nu/Dh;  
end


%% Calc h_c coefficient
function h_c = calc_h_c(rho, k, mu, Pr, beta, D_h, Th_in, Th_out, Tc)
	g    = 9.81;				% Gravitational acceleration [m/s^2]
	Ts_h = (Th_in + Th_out)/2;	% Avg. pipe surface temperature [K]
	
	% Rayleigh number
	Ra = ( (g*beta*(Ts_h - Tc)*D_h^3)/((mu/rho)^2) )*Pr;
	
	% Nusselt correl. for nat. convect. of horizontal cylinder [Donald and Chu]
	Nu = (  0.6  +  ( 0.387*Ra^(1/6) ) / ( (1 + (0.559/Pr)^(9/16))^(8/27) )  )^2;
	
	% Heat transfer coefficient [W/(m^2 K)]
	h_c = k*Nu/D_h;  
end


%% Blend function
function y = blend(y1,y2,x1,x2,x)
% Blend between two values y1 and y2 as x varies between x1 and x2. The
% blending function ensures y is continuous and differentiable in the
% range x1 <= x <= x2.

% Copyright 2017-2018 The MathWorks, Inc.

% [Source: simscape.function.blend.scc]
	u = (x-x1)/(x2-x1);
	t = 3*u^2 - 2*u^3;
    y = (1-t)*y1 + t*y2;
end