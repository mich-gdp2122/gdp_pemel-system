function [pmin, pmax, v3, mdot, Ac, D, y1, y3] = ...
	ORCspec(Tmin, Tmax, x1, x3, mdot_h, Th_in, Th_out)
% Calculate mass flow and pipe area requirements for ORC
%%
% Get relevant state properties
[rho1,  ~, ~,~,~,~,~, pmin, c1] = data_r600a_sat(Tmin, x1); % State 1
[rho2,  ~, ~,~,~,~,~,    ~, c2] = data_r600a_sat(Tmax, 0);  % State 2 (sat L)
[rho3, v3, ~,~,~,~,~, pmax, c3] = data_r600a_sat(Tmax, x3); % State 3
[rho4,  ~, ~,~,~,~,~,    ~, c4] = data_r600a_sat(Tmin, 1);  % State 4 (sat v)

[rho23, ~,~,~,~,~, dh_23, ~, c23] = data_r600a_sat(Tmax, 0, x3);  % Process 2->3
[rho41, ~,~,~,~,~,     ~, ~, c41] = data_r600a_sat(Tmin, x1, 1);  % Process 4->1

% Vapour volume fractions
y1 = y(x1, rho2, rho4);
y3 = y(x3, rho2, rho4);


%% Mass flow required to transfer coolant heat
% Coolant cp [J/(kg*K)] and heat transfer req'd [W]
[~, cp_h, ~, ~, ~] = data_water(Th_in, Th_out);
Qdot_h = mdot_h*cp_h*abs(Th_in - Th_out);

% ORC mass flow req'd, via energy balance [kg/s]
mdot = Qdot_h/abs(dh_23);


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
end
%%
function y_out = y(x, rho_sL, rho_sv)
% Convert vapour quality/mass fraction to volume fraction
	y_out = ( (1/rho_sv)*x )/( (1/rho_sv)*x + (1/rho_sL)*(1 - x) );
end