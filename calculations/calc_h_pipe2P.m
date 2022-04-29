function h = calc_h_pipe2P(sfCase, k, mu, Dh, mdot, T, xA, xB)
% Calculate heat transfer coefficient for liquid flow in pipe
%  T: Saturation temperaturedata_c.mu
% xA: Inlet quality/vapour mass fraction
% xB: Outlet quality/vapour mass fraction

% Saturated properties
data_sL = data_r600a_sat(T, 0);  % Liquid
data_sV = data_r600a_sat(T, 1);  % Vapour

% Reynold's number
Re    = (4*mdot)/(pi*Dh*mu);	% Avg. 2P fluid
Re_sL = (4*mdot)/(pi*Dh*data_sL.mu);	% Sat. liquid

% Reynold's no. transition region
Re_lam = 2300;
Re_tur = 3000;

%% Coolant Nusselt number []
% Laminar value
switch sfCase
	case 'Qs'
		Nu_lam = 4.36;	% Const. surface Q
	case 'Ts'
		Nu_lam = 3.66;	% Const. surface T
	otherwise
		error(['Unknown pipe surface heat profile, ''', sfCase, '''.', newline ...
			  '(Choose ''Qs'' for const. surface heat flux, or ',...
			  '''Ts'' for const. surface temperature)']);
end

% Turbulent value (Cavallini and Zecchin)
Nu_tur = ...
	(	0.05 * Re_sL^0.8 * data_sL.Pr^0.33 * ...
		(  ( (sqrt(data_sL.rho/data_sV.rho) - 1)*xB + 1 )^1.8 - ...
		   ( (sqrt(data_sL.rho/data_sV.rho) - 1)*xA + 1 )^1.8   ...
		) ...
	) / ( 1.8*(sqrt(data_sL.rho/data_sV.rho) - 1)*(xB - xA) );

if Re <= Re_lam						% Laminar case
	Nu = Nu_lam;
elseif Re >= Re_tur					% Turbulent case
	Nu = Nu_tur;
elseif Re > Re_lam && Re < Re_tur	% Transitional case (via interpolation)
	Nu = blend(Nu_lam, Nu_tur, Re_lam, Re_tur, Re);
end
%%
% Heat transfer coefficient [W/(m^2 K)]
h = k*Nu/Dh;  
end