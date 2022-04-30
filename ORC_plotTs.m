clearvars;
close all;
parameters_straight;
% Generate saturation lines
Tsat = [257:2:407, 407.01:0.1:407.81];
for i = 1:length(Tsat)
	s_satL(i) = Ts_sat_r600a(Tsat(i), 'T_L');
	s_satV(i) = Ts_sat_r600a(Tsat(i), 'T_v');
end
% % Blend between gap
% s_crit = linspace(s_satL(end),s_satV(end));
% for i = 1:length(s_crit)
% 	Tcrit(i) = interp1([s_satL(end), s_satV(end)], )
% 	%Tcrit(i) = blend(Tsat(end), Tsat(end), s_satL(end), s_satV(end), s_crit(i));
% end

% Combine into single array
Tsat  = [Tsat, flip(Tsat)];
s_sat = [s_satL, flip(s_satV)];

figure
hold on
grid on
plot(s_sat,Tsat,'k-', 'LineWidth', 1.6);
pl_Ts = plot(ORC.plotTs.s,ORC.plotTs.T, 'o-', ...
	'LineWidth',1.1, 'MarkerSize',6, 'Color','#4DBEEE');
%plot(ORC.plotTs.sh,ORC.plotTs.Th, 'r-', 'LineWidth',0.8)
pl_Th = plot(ORC.plotTs.sh_n,ORC.plotTs.Th_n, 'ro--');
%pl_Tc = plot(ORC.plotTs.sc_n,ORC.plotTs.Tc_n, 'o--', 'Color', '#0072BD');
yline(ORC.plotTs.Tc_n,'--', ...
	'Color','#0072BD', 'Interpreter', 'tex', 'LabelHorizontalAlignment', 'left');
ylim([270, 360])
xlim([800, 2600])
ylabel('T [K]')
xlabel('s [J/(kg*K)]')
lbl_Ts   =				  {'1', '2', '2f', '3', '4', ''};
lbl_Ts_y = ORC.plotTs.T + [-6,   1,  -1.5, 1.5, -6,   0];
lbl_ts_x = ORC.plotTs.s + [50, -15,   -40,  20, 60,   0];
lbl_Th   =                   {'T_{h,in}', 'T_{pp}', 'T_{h,out}'};
lbl_Th_y = ORC.plotTs.Th_n +  2 + zeros(1,3);
lbl_Th_x = ORC.plotTs.sh_n + 40 + zeros(1,3);
lbl_Tc   = 'T_{sea}';
lbl_Tc_y = ORC.plotTs.Tc_n;
lbl_Tc_x = 930; %mean(ORC.plotTs.sc_n);
labelpoints(lbl_ts_x,lbl_Ts_y,lbl_Ts, 'FontSize',11, 'FontWeight','bold');
labelpoints(lbl_Th_x,lbl_Th_y,lbl_Th, 'Interpreter','tex', 'FontSize',10);
labelpoints(lbl_Tc_x,lbl_Tc_y,lbl_Tc, 'Interpreter','tex', 'FontSize',10);
% lbl_Ts   =				  {'1', '2', '2f',   '3', '4', ''};
% lbl_Ts_y = ORC.plotTs.T + [-4,   5,    1,     4,  -4,   0];
% lbl_ts_x = ORC.plotTs.s + [10, -30, -100,   -10,  30,   0];
% lbl_Th   =                   {'T_{h,in}', 'T_{pp}', 'T_{h,out}'};
% lbl_Th_y = ORC.plotTs.Th_n + [  5,          5,        5];
% lbl_Th_x = ORC.plotTs.sh_n + [-75,        -60,        -70];
% text(lbl_ts_x,lbl_Ts_y,lbl_Ts, 'FontSize',11);
% text(lbl_Th_x,lbl_Th_y,lbl_Th, 'Interpreter','tex', 'FontSize',10);

function out = Ts_sat_r600a(X, X_ph)
	% r600a sp entropy on saturation curve
	T_i = [257:5:407, 407.01:0.1:407.81];	% Temperature vector [K]
	
	s_L_i = 1E3*[	0.86368, ...
					0.90626, ...
					0.94849, ...
					0.9904, ...
					1.032, ...
					1.0734, ...
					1.1144, ...
					1.1553, ...
					1.196, ...
					1.2365, ...
					1.2768, ...
					1.317, ...
					1.3571, ...
					1.3972, ...
					1.4372, ...
					1.4771, ...
					1.5171, ...
					1.5572, ...
					1.5973, ...
					1.6376, ...
					1.6781, ...
					1.7189, ...
					1.76, ...
					1.8017, ...
					1.844, ...
					1.8873, ...
					1.9318, ...
					1.9783, ...
					2.0279, ...
					2.0839, ...
					2.1651, ...
					2.1654, ...
					2.1683, ...
					2.1715, ...
					2.1751, ...
					2.179, ...
					2.1837, ...
					2.1894, ...
					2.1976, ...
					2.2259];	% r600a sp. entropy (sat liq) [J/(kg*K)]
	
	s_v_i = 1E3*[	2.2994, ...
					2.2977, ...
					2.297, ...
					2.2971, ...
					2.2979, ...
					2.2995, ...
					2.3017, ...
					2.3044, ...
					2.3077, ...
					2.3114, ...
					2.3155, ...
					2.32, ...
					2.3247, ...
					2.3297, ...
					2.3349, ...
					2.3402, ...
					2.3455, ...
					2.3508, ...
					2.356, ...
					2.361, ...
					2.3657, ...
					2.3699, ...
					2.3735, ...
					2.3764, ...
					2.3782, ...
					2.3787, ...
					2.3772, ...
					2.3729, ...
					2.3641, ...
					2.3462, ...
					2.2937, ...
					2.2934, ...
					2.2906, ...
					2.2875, ...
					2.2839, ...
					2.2798, ...
					2.2748, ...
					2.2684, ...
					2.2589, ...
					2.2259];	% r600a sp. entropy (sat vap) [J/(kg*K)]
	
	% Output s_sat(Tsat) or Tsat(s_sat), depending on spec'd case
	switch X_ph
		case 'T_L'
			out = interp1(T_i, s_L_i, X, 'makima', 'extrap');
		case 'T_v'
			out = interp1(T_i, s_v_i, X, 'makima', 'extrap');
		case 's_L'
			out = interp1(s_L_i, T_i, X, 'makima', 'extrap');
		case 's_v'
			out = interp1(s_v_i, T_i, X, 'makima', 'extrap');
	end
end