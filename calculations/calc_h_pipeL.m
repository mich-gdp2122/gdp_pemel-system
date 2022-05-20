function h = calc_h_pipeL(sfCase, k, mu, Pr, Dh, mdot)
% FEEG6013 Group Design Project, 2021-2022
% Group 19
%
% Created by Michael
%
%
% Calculate heat transfer coefficient for liquid flow in pipe
%
%
%%
% Reynold's number
Re = (4*mdot)/(pi*Dh*mu);

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
switch sfCase
	case 'Qs'
		Nu_lam = 4.36;						% Laminar case (const. surface Q)
	case 'Ts'
		Nu_lam = 3.66;						% Laminar case (const. surface T)
	otherwise
		error(['Unknown pipe surface heat profile, ''', sfCase, '''.', newline ...
			  '(Choose ''Qs'' for const. surface heat flux, or ',...
			  '''Ts'' for const. surface temperature)']);
end
	Nu_tur = ( (f/8)*(Re - 1000)*Pr )/...	% Turbulent case (Gnielinski)
		 ( 1 + 12.7*((f/8)^0.5)*( Pr^(2/3) - 1 ) );  

if Re <= Re_lam								% Laminar case
	Nu = Nu_lam;
elseif Re >= Re_tur							% Turbulent case
	Nu = Nu_tur;
elseif Re > Re_lam && Re < Re_tur			% Transitional case (via interpolation)
	Nu = blend(Nu_lam, Nu_tur, Re_lam, Re_tur, Re);
end

% Heat transfer coefficient [W/(m^2 K)]
h = k*Nu/Dh;  
end