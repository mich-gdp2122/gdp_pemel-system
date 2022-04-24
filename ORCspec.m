function ORCstruct = ORCspec(dTc, dTh, mdot_h, Th_in, Th_out, Tc, eff_pmp, eff_tbn)
% Calculate mass flow and pipe area requirements for ORC
%% Calculate roperties
% Max & min cycle temperatures, via driving temps [K]
Tmax = Th_out - dTh;
Tmin = Tc + dTc;

% State properties on saturation line
[rho2, ~,~,~,~,~,~, pmax, c2, h2, s2] = data_r600a_sat(Tmax, 0);  % State 2 (sat L)
[rho4, ~,~,~,~,~,~, pmin, c4, h4, s4] = data_r600a_sat(Tmin, 1);  % State 4 (sat v)

% State properties in sat region, via isentropic processes
[rho1, c1, h1, x1] = isentrop(Tmin, s2, 0);  % State 1, via process 1->2
[rho3, c3, h3, x3] = isentrop(Tmax, s4, 1);  % State 1, via process 1->2
%[rho1, ~,~,~,~,~,~,~, c1, h1] = data_r600a_sat(Tmin, x1); % State 1
%[rho3, ~,~,~,~,~,~,~, c3, h3] = data_r600a_sat(Tmax, x3); % State 3

[rho23, ~,~,~,~,~,~,~, c23] = data_r600a_sat(Tmax, 0, x3);  % Process 2->3
[rho41, ~,~,~,~,~,~,~, c41] = data_r600a_sat(Tmin, x1, 1);  % Process 4->1

% Sp. enthalpy changes [J/kg]
dh_12 = h2 - h1;
dh_23 = h3 - h2;
dh_34 = h4 - h3;
dh_41 = h1 - h4;

% Vapour volume fractions
% y1 = y(x1, rho2, rho4);
% y3 = y(x3, rho2, rho4);


%% Mass flow required to transfer coolant heat
% Coolant cp [J/(kg*K)] and heat transfer req'd [W]
[~, cp_h, ~, ~, ~] = data_water(Th_in, Th_out);
Qin = mdot_h*cp_h*abs(Th_in - Th_out);

% ORC mass flow req'd, via energy balance [kg/s]
mdot = Qin/abs(dh_23);

% Heat rejection req'd [W]
Qout = mdot*abs(dh_41);


%% Pipe area required to maintain incompressible flow
Ma_max = 0.27;		% Max incompressible Mach number

% Minimum rho*c
rhoc_min = ...
	min([rho1*c1, ...					% State 1
		 rho2*c2, ...					% State 2
		 rho3*c3, ...					% State 3
		 rho4*c4, ...					% State 4
		 mean([rho1*c1 rho2*c2]), ...	% Process 1->2		 
		 rho23*c23, ...					% Process 2->3
		 mean([rho3*c3 rho4*c4]), ...	% Process 3->4
		 rho41*c41, ...					% Process 4->1
		 ]);

% Required diameter [m] and cross-section area [m^2]
D  = sqrt( (4*mdot)/(pi*rhoc_min*Ma_max) );
Ac = pi*(D/2)^2;

%% Calculate cycle performance
% Pump and turbine power [W]
pwr_pmp  = (1/eff_pmp)*mdot*abs(dh_12);
pwr_tbne = eff_tbn*mdot*abs(dh_34); 
% Thermal efficiency [%]
eff = 100*(pwr_tbne - pwr_pmp)/Qin;


%% Put all into specified ORC structure
ORCstruct.dTh      = dTh;
ORCstruct.dTc      = dTc;
ORCstruct.Tmin     = Tmin;
ORCstruct.Tmax     = Tmax;
ORCstruct.pmin     = pmin;
ORCstruct.pmax     = pmax; 
ORCstruct.x1       = x1;
ORCstruct.x3       = x3; 
ORCstruct.mdot     = mdot; 
ORCstruct.Ac       = Ac; 
ORCstruct.D        = D; 
ORCstruct.pwr_pmp  = pwr_pmp; 
ORCstruct.pwr_tbne = pwr_tbne; 
ORCstruct.Qin      = Qin; 
ORCstruct.Qout     = Qout;
ORCstruct.eff      = eff;
end


%%
function [rhoB, cB, hB, xB] = isentrop(TB, sA, xA)
% Calc. outlet properties of isentropic process
	if xA == 1
		x = flip(0:01:1);	% Search from sat v line
	else
		x = 0:01:1;			% Search from sat L line otherwise
	end
	% Initial guess
	[~,~,~,~,~,~,~,~,~,~,sB_prev] = data_r600a_sat(TB, x(1));
	% Cycle thru vapour qualities until sp. entropy out = sp. entropy in
	for i = 2:length(x)
		[~,~,~,~,~,~,~,~,~,~, sB] = data_r600a_sat(TB, x(i));
		if (sB >= sA && sB_prev <= sA) || (sB <= sA && sB_prev >= sA)
			% Interp for xB between found sB vals of interest
			xB = interp1([sB_prev, sB], [x(i-1), x(i)], sA, 'linear');
			% Get output data
			[rhoB, ~,~,~,~,~,~,~, cB, hB, sB] = data_r600a_sat(TB, xB);
			break
		else
			% Otherwise set current sB val as prev then continue loop
			sB_prev = sB;
		end
	end
end

% function y_out = y(x, rho_sL, rho_sv)
% % Convert vapour quality/mass fraction to volume fraction
% 	y_out = ( (1/rho_sv)*x )/( (1/rho_sv)*x + (1/rho_sL)*(1 - x) );
% end