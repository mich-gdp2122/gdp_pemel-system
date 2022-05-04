% Main data-gathering script for TMS models
%
% SIMULINK MODEL OUTPUTS:
%	
%	TMS(i)_pwr   [stk, BoP, tot]				Stack, Total BoP, Total [kW]
%	TMS(i)_BoP   [pmpFW, pmpClnt, htr]		H2O pump, Clnt pump, heater [W]
%	TMS(i)_Qdot [stkIn, PHOut, rjOut]		Stack in, Phtr out, RJ out [kW]
%	TMS(i)_eff   [stk, BoP, tot]				Stack, BoP, Overall [%]
%   TMS(i)_ORC   [P_Grs, P_net, Qin, eff]	
%
%% Set program modes
mode_test = 0;	% Test mode (don't run simulations if != 0)

%%
clearvars('-except', 'mode_*')
close all

% Make output directory
outdir = ['output_', datestr(now,'yyyy-mm-dd--HH-MM-SS')];
mkdir(outdir);


%% Run simulations
% Run simulations if not testing
if mode_test == 0
	runSim('straight', outdir);					% 1) STRAIGHT CHANNELS
	clearvars('-except', 'mode_*', 'outdir')	% Clear for next run
	runSim('serpentine', outdir);				% 2) SERPENTINE CHANNELS
	clearvars('-except', 'mode_*', 'outdir')	% Clear for loading data
	disp('Simulations complete')
end
% Load simulation data
strt = load('output_straight.mat');
sptn = load('output_serpentine.mat');


%% Figure properties
b_h  = 360;  % Height [px]
b_ar = 1.8;  % Aspect ratio [w:h]

% Bar chart colours (RGB)
c.red    = [0.6350 0.0780 0.1840];
c.blue   = [0.0000 0.4470 0.7410];
c.yellow = [0.9290 0.6940 0.1250];
c.green  = [0.4660 0.6740 0.1880];
c.orange = [0.8500 0.3250 0.0980];
c.cyan   = [0.3010 0.7450 0.9330];
c.purple = [0.4940 0.1840 0.5560];


%% Total BoP Pwr Consumption data [kW]
% Organise data and round it
y_totBoP  = round([strt.TMS1_pwr(2), strt.TMS2_pwr(2), strt.TMS3_pwr(2), strt.TMS4_pwr(2); ...
				   sptn.TMS1_pwr(2), sptn.TMS2_pwr(2), sptn.TMS3_pwr(2), sptn.TMS4_pwr(2)], ...
			3, 'significant');
% Bar colours and respective legend labels
c_totBoP = [c.orange; c.blue; c.yellow; c.green];
l     = {'Config 1', 'Config 2', 'Config 3', 'Config 4'};
% Plot data
plotBarData(outdir, y_totBoP, c_totBoP, 'BoP power consumption [kW]', [],...
	'grouped', [], [], '1_totBoP', l, b_h, b_ar);


%% BoP Power Consumption data [% tot]
% Organise data and round it
y1_BoP = round(0.1*[strt.TMS1_BoP/strt.TMS1_pwr(2); sptn.TMS1_BoP/sptn.TMS1_pwr(2)], ...
			3, 'significant');
y2_BoP = round(0.1*[strt.TMS2_BoP/strt.TMS2_pwr(2); sptn.TMS2_BoP/sptn.TMS2_pwr(2)], ...
			3, 'significant');
y3_BoP = round(0.1*[strt.TMS3_BoP/strt.TMS3_pwr(2); sptn.TMS3_BoP/sptn.TMS3_pwr(2)], ...
			3, 'significant');
y4_BoP = round(0.1*[strt.TMS4_BoP/strt.TMS4_pwr(2); sptn.TMS4_BoP/sptn.TMS4_pwr(2)], ...
			3, 'significant');
% Set any NaN values to 0
y1_BoP(isnan(y1_BoP)) = 0;
y2_BoP(isnan(y2_BoP)) = 0;
y3_BoP(isnan(y3_BoP)) = 0;
y4_BoP(isnan(y4_BoP)) = 0;
% Figure dimensions
p_h = 900;
p_ar = 0.85;
p_ar2 = 1.6;
% Bar colours and respective legend labels
c_BoP = [c.blue; c.cyan; c.orange; c.green];
l     = {'Water pump', 'Coolant pump', 'Heater', 'ORC pump'};
% Title subtitles
subtitles_totBoP = append('Total BoP pwr: ', cellfun(@num2str, num2cell(y_totBoP), 'UniformOutput', false),' kW');
%subtitles1 = append('Total BoP pwr: ',{num2str(y_totBoP(1,1)); num2str(y_totBoP(1,2)); num2str(y_totBoP(1,3)); num2str(y_totBoP(1,4))},' kW');
% Plot data
plotPieData(outdir, y1_BoP(1,:), y2_BoP(1,:), y3_BoP(1,:), y4_BoP(1,:), c_BoP, 1,...
			'Straight Channels', [], subtitles_totBoP(1,:), '2A_BoP_strt', l, p_h, p_ar);
plotPieData(outdir, y1_BoP(2,:), y2_BoP(2,:), y3_BoP(2,:), y4_BoP(2,:), c_BoP, 1, ...
			'Serpentine Channels', [], subtitles_totBoP(2,:), '2B_BoP_sptn', l, p_h, p_ar);

plotPieData2(outdir, y1_BoP, y2_BoP, y3_BoP, y4_BoP, c_BoP, 1, subtitles_totBoP,...
			 '2_BoP', l, p_h, p_ar2);


%% Heat recovery data [kW]
% Organise data and round it
y1_Qdot(1:2,:) = round([strt.TMS1_Qdot(2:end); ...
						strt.TMS2_Qdot(2:end) ...
					   ], 4, 'significant');
y1_Qdot(3,:)   = round([strt.TMS3_Qdot(2), ...
						strt.TMS3_Qdot(3)-abs(strt.TMS3_Qdot(4)), ...
						strt.TMS3_Qdot(4)...
					   ], 4, 'significant');
y1_Qdot(4,:)   = round([strt.TMS4_Qdot(2), ...
						strt.TMS4_Qdot(3)-abs(strt.TMS4_Qdot(4)), ...
						strt.TMS4_Qdot(4)...
					   ], 4, 'significant');

y2_Qdot(1:2,:) = round([sptn.TMS1_Qdot(2:end); ...
						sptn.TMS2_Qdot(2:end) ...
					   ], 4, 'significant');
y2_Qdot(3,:)   = round([sptn.TMS3_Qdot(2), ...
					    sptn.TMS3_Qdot(3)-abs(sptn.TMS3_Qdot(4)), ...
						sptn.TMS3_Qdot(4)...
					   ], 4, 'significant');
y2_Qdot(4,:)   = round([sptn.TMS4_Qdot(2), ...
						sptn.TMS4_Qdot(3)-abs(sptn.TMS4_Qdot(4)), ...
						sptn.TMS4_Qdot(4)...
					   ], 4, 'significant');

% Figure dimensions
bQ_h = 560;
bQ_ar = 2;
% Bar colours and respective legend labels
c_Qdot = [c.cyan; c.yellow; c.blue];
l       = {'Preheater', 'ORC net pwr', 'Rejected'};
% Plot titles
subtitle1_Qdot = ['PEMEL heat load: ', num2str(round(strt.TMS1_Qdot(1), 3, 'significant')), ' kW'];
subtitle2_Qdot = ['PEMEL heat load: ', num2str(round(sptn.TMS1_Qdot(1), 3, 'significant')), ' kW'];
 % Plot data
plotBarData2(outdir, y1_Qdot, y2_Qdot, c_Qdot, 'Heat recovered [kW]', ...
	subtitle1_Qdot, subtitle2_Qdot,'stacked', [], [], '3_Qdot', l, bQ_h, bQ_ar);


%% Heat recovery data [% stack in]
% Organise data and round it
y1_QdotPct(1:2,:) = round(100*[abs(strt.TMS1_Qdot(2:end))/strt.TMS1_Qdot(1); ...
							   abs(strt.TMS2_Qdot(2:end))/strt.TMS2_Qdot(1)...
							  ], 4, 'significant');
y1_QdotPct(3,:)   = round(100*[strt.TMS3_Qdot(2), ...
							   strt.TMS3_Qdot(3)-abs(strt.TMS3_Qdot(4)), ...
							   abs(strt.TMS3_Qdot(4))...
							  ]/strt.TMS3_Qdot(1), 4, 'significant');
y1_QdotPct(4,:)   = round(100*[strt.TMS4_Qdot(2), ...
							   strt.TMS4_Qdot(3)-abs(strt.TMS4_Qdot(4)), ...
							   abs(strt.TMS4_Qdot(4))...
							  ]/strt.TMS4_Qdot(1), 4, 'significant');

y2_QdotPct(1:2,:) = round(100*[abs(sptn.TMS1_Qdot(2:end))/sptn.TMS1_Qdot(1); ...
							   abs(sptn.TMS2_Qdot(2:end))/sptn.TMS2_Qdot(1)...
							  ], 4, 'significant');
y2_QdotPct(3,:)   = round(100*[sptn.TMS3_Qdot(2), ...
							   sptn.TMS3_Qdot(3)-abs(sptn.TMS3_Qdot(4)), ...
							   abs(sptn.TMS3_Qdot(4))...
							  ]/sptn.TMS3_Qdot(1), 4, 'significant');
y2_QdotPct(4,:)   = round(100*[sptn.TMS4_Qdot(2), ...
							   sptn.TMS4_Qdot(3)-abs(sptn.TMS4_Qdot(4)), ...
							   abs(sptn.TMS4_Qdot(4))...
							  ]/sptn.TMS4_Qdot(1), 4, 'significant');
% Set any NaN values to 0
y1_QdotPct(isnan(y1_QdotPct)) = 0;
y2_QdotPct(isnan(y2_QdotPct)) = 0;

% Figure dimensions
pQ_h = 720;
pQ_ar = 0.9;
pQ_ar2 = 1.8;
% Bar colours and respective legend labels
c_QdotPct = c_Qdot;
l       = {'Preheater', 'ORC net pwr', 'Rejected'};
% Plot data
plotPieData(outdir, y1_QdotPct(1,:), y1_QdotPct(2,:), y1_QdotPct(3,:), y1_QdotPct(4,:), c_QdotPct, 2,...
			'Straight Channels', subtitle1_Qdot, [], '4A_QdotPct_strt', l, pQ_h, pQ_ar);
plotPieData(outdir, y2_QdotPct(1,:), y2_QdotPct(2,:), y2_QdotPct(3,:), y2_QdotPct(4,:), c_QdotPct, 2,...
			'Serpentine Channels', subtitle2_Qdot, [], '4B_QdotPct_sptn', l, pQ_h, pQ_ar);

plotPieData2(outdir, [y1_QdotPct(1,:); y2_QdotPct(1,:)], [y1_QdotPct(2,:); y2_QdotPct(2,:)], ...
			 [y1_QdotPct(3,:); y2_QdotPct(3,:)], [y1_QdotPct(4,:); y2_QdotPct(4,:)], c_QdotPct, 2, [], ...
			 '4_QdotPct', l, pQ_h, pQ_ar2);


%% Efficiency data [%]
% Organise data and round it
y_eff  = round([strt.TMS1_eff(2), strt.TMS2_eff(2), strt.TMS3_eff(2), strt.TMS4_eff(2); ...
				   sptn.TMS1_eff(2), sptn.TMS2_eff(2), sptn.TMS3_eff(2), sptn.TMS4_eff(2)], ...
			4, 'significant');
% Bar colours and respective legend labels
c_eff = c_totBoP;
l     = {'Config 1', 'Config 2', 'Config 3', 'Config 4'};
b_ar2 = 2.3;
% Title
title_eff = ['PEMEL efficiency: ', num2str(round(strt.TMS1_eff(1), 3, 'significant')), '%'];
% Plot data
plotBarData(outdir, y_eff, c_eff, 'Overall efficiency [%]', title_eff, ...
	'grouped', [], [], '5_eff', l, b_h, b_ar2);


%% Total Pwr consumption data [kW]
% Organise data and round it
y_totPwr  = round([strt.TMS1_pwr(3), strt.TMS2_pwr(3), strt.TMS3_pwr(3), strt.TMS4_pwr(3); ...
				   sptn.TMS1_pwr(3), sptn.TMS2_pwr(3), sptn.TMS3_pwr(3), sptn.TMS4_pwr(3)], ...
			3, 'significant');
% Bar colours and respective legend labels
c_totPwr = c_totBoP;
l     = {'Config 1', 'Config 2', 'Config 3', 'Config 4'};
% Title
title_pwr = ['PEMEL power consumption: ', num2str(round(strt.TMS1_pwr(1), 3, 'significant')), ' kW'];
% Plot data
plotBarData(outdir, y_totPwr, c_totPwr, 'Overall power consumption [kW]', title_pwr,...
	'grouped', [], [], '6_totPwr', l, b_h, b_ar2);


%% ORC T-s plot

orc_h = 750;
orc_ar = 4/3;

% Generate saturation lines
Tsat = [257:2:407, 407.01:0.1:407.81];
for i = 1:length(Tsat)
	s_satL(i) = Ts_sat_r600a(Tsat(i), 'T_L');
	s_satV(i) = Ts_sat_r600a(Tsat(i), 'T_v');
end
% Combine into single array
Tsat  = [Tsat, flip(Tsat)];
s_sat = [s_satL, flip(s_satV)];


% Create plots
fig_ORC = figure('Position',[300,100, orc_h*orc_ar, orc_h]);
t_ORC   = tiledlayout(fig_ORC,2,2, 'TileSpacing','compact', 'Padding','tight');
t_ORC.YLabel.String = 'Temperature, T [K]';
t_ORC.XLabel.String = 'Specific entropy, s [J kg^{-1} K^{-1}]';

% Config 3 (straight)
t11_ORC = nexttile;
title(t11_ORC, 'Config 3, Straight Channels')
hold on
grid on
plot(t11_ORC,s_sat,Tsat,'k-', 'LineWidth', 1.6);
pl_Ts = plot(t11_ORC,strt.ORC.plotTs.s,strt.ORC.plotTs.T, 'o-', ...
	'LineWidth',1.1, 'MarkerSize',6, 'Color','#D95319');
pl_Th = plot(t11_ORC,strt.ORC.plotTs.sh_n,strt.ORC.plotTs.Th_n, 'o--', 'Color','#A2142F');
yline(t11_ORC,strt.ORC.plotTs.Tc_n,'--', ...
	'Color','#0072BD', 'Interpreter', 'tex', 'LabelHorizontalAlignment', 'left');
ylim(t11_ORC,[270, 360])
xlim(t11_ORC,[800, 2600])
lbl_Ts   =					   {'1', '2', '2f', '3', '4', ''};
lbl_Ts_y = strt.ORC.plotTs.T + [-6,   1,  -1.5, 1.5, -6,   0];
lbl_ts_x = strt.ORC.plotTs.s + [50, -15,   -40,  20, 60,   0];
lbl_Th   = {'T_{h,in}', 'T_{pp}', 'T_{h,out}'};
lbl_Th_y = strt.ORC.plotTs.Th_n +  2 + zeros(1,3);
lbl_Th_x = strt.ORC.plotTs.sh_n + 40 + zeros(1,3);
lbl_Tc   = 'T_{sea}';
lbl_Tc_y = strt.ORC.plotTs.Tc_n;
lbl_Tc_x = 930;
labelpoints(lbl_ts_x,lbl_Ts_y,lbl_Ts, 'FontSize',11, 'FontWeight','bold');
labelpoints(lbl_Th_x,lbl_Th_y,lbl_Th, 'Interpreter','tex', 'FontSize',10);
labelpoints(lbl_Tc_x,lbl_Tc_y,lbl_Tc, 'Interpreter','tex', 'FontSize',10);

% Config 4 (straight)
t12_ORC = nexttile;
title(t12_ORC, 'Config 4, Straight Channels')
hold on
grid on
plot(t12_ORC,s_sat,Tsat,'k-', 'LineWidth', 1.6);
pl_Ts = plot(t12_ORC,strt.ORCph.plotTs.s,strt.ORCph.plotTs.T, 'o-', ...
	'LineWidth',1.1, 'MarkerSize',6, 'Color','#D95319');
pl_Th = plot(t12_ORC,strt.ORCph.plotTs.sh_n,strt.ORCph.plotTs.Th_n, 'o--', 'Color','#A2142F');
yline(t12_ORC,strt.ORCph.plotTs.Tc_n,'--', ...
	'Color','#0072BD', 'Interpreter', 'tex', 'LabelHorizontalAlignment', 'left');
ylim(t12_ORC,[270, 360])
xlim(t12_ORC,[800, 2600])
lbl_Ts   =					   {'1', '2', '2f', '3', '4', ''};
lbl_Ts_y = strt.ORCph.plotTs.T + [-6,   1,  -1.5, 1.5, -6,   0];
lbl_ts_x = strt.ORCph.plotTs.s + [50, -15,   -40,  20, 60,   0];
lbl_Th   = {'T_{h,in}', 'T_{pp}', 'T_{h,out}'};
lbl_Th_y = strt.ORCph.plotTs.Th_n +  2 + zeros(1,3);
lbl_Th_x = strt.ORCph.plotTs.sh_n + 40 + zeros(1,3);
lbl_Tc   = 'T_{sea}';
lbl_Tc_y = strt.ORCph.plotTs.Tc_n;
lbl_Tc_x = 930;
labelpoints(lbl_ts_x,lbl_Ts_y,lbl_Ts, 'FontSize',11, 'FontWeight','bold');
labelpoints(lbl_Th_x,lbl_Th_y,lbl_Th, 'Interpreter','tex', 'FontSize',10);
labelpoints(lbl_Tc_x,lbl_Tc_y,lbl_Tc, 'Interpreter','tex', 'FontSize',10);

% Config 3 (serpentine)
t21_ORC = nexttile;
title(t21_ORC, 'Config 3, Serpentine Channels')
hold on
grid on
plot(t21_ORC,s_sat,Tsat,'k-', 'LineWidth', 1.6);
pl_Ts = plot(t21_ORC,sptn.ORC.plotTs.s,sptn.ORC.plotTs.T, 'o-', ...
	'LineWidth',1.1, 'MarkerSize',6, 'Color','#D95319');
pl_Th = plot(t21_ORC,sptn.ORC.plotTs.sh_n,sptn.ORC.plotTs.Th_n, 'o--', 'Color','#A2142F');
yline(t21_ORC,sptn.ORC.plotTs.Tc_n,'--', ...
	'Color','#0072BD', 'Interpreter', 'tex', 'LabelHorizontalAlignment', 'left');
ylim(t21_ORC,[270, 360])
xlim(t21_ORC,[800, 2600])
lbl_Ts   =					   {'1', '2', '2f', '3', '4', ''};
lbl_Ts_y = sptn.ORC.plotTs.T + [-6,   1,  -1.5, 1.5, -6,   0];
lbl_ts_x = sptn.ORC.plotTs.s + [50, -15,   -40,  20, 60,   0];
lbl_Th   = {'T_{h,in}', 'T_{pp}', 'T_{h,out}'};
lbl_Th_y = sptn.ORC.plotTs.Th_n +  2 + zeros(1,3);
lbl_Th_x = sptn.ORC.plotTs.sh_n + 40 + zeros(1,3);
lbl_Tc   = 'T_{sea}';
lbl_Tc_y = sptn.ORC.plotTs.Tc_n;
lbl_Tc_x = 930;
labelpoints(lbl_ts_x,lbl_Ts_y,lbl_Ts, 'FontSize',11, 'FontWeight','bold');
labelpoints(lbl_Th_x,lbl_Th_y,lbl_Th, 'Interpreter','tex', 'FontSize',10);
labelpoints(lbl_Tc_x,lbl_Tc_y,lbl_Tc, 'Interpreter','tex', 'FontSize',10);

% Config 4 (serpentine)
t22_ORC = nexttile;
title(t22_ORC,'Config 4, Serpentine Channels')
hold on
grid on
plot(t22_ORC,s_sat,Tsat,'k-', 'LineWidth', 1.6);
pl_Ts = plot(t22_ORC,sptn.ORCph.plotTs.s,sptn.ORCph.plotTs.T, 'o-', ...
	'LineWidth',1.1, 'MarkerSize',6, 'Color','#D95319');
pl_Th = plot(t22_ORC,sptn.ORCph.plotTs.sh_n,sptn.ORCph.plotTs.Th_n, 'o--', 'Color','#A2142F');
yline(t22_ORC,sptn.ORCph.plotTs.Tc_n,'--', ...
	'Color','#0072BD', 'Interpreter', 'tex', 'LabelHorizontalAlignment', 'left');
ylim(t22_ORC,[270, 360])
xlim(t22_ORC,[800, 2600])
lbl_Ts   =					   {'1', '2', '2f', '3', '4', ''};
lbl_Ts_y = sptn.ORCph.plotTs.T + [-6,   1,  -1.5, 1.5, -6,   0];
lbl_ts_x = sptn.ORCph.plotTs.s + [50, -15,   -40,  20, 60,   0];
lbl_Th   = {'T_{h,in}', 'T_{pp}', 'T_{h,out}'};
lbl_Th_y = sptn.ORCph.plotTs.Th_n +  2 + zeros(1,3);
lbl_Th_x = sptn.ORCph.plotTs.sh_n + 40 + zeros(1,3);
lbl_Tc   = 'T_{sea}';
lbl_Tc_y = sptn.ORCph.plotTs.Tc_n;
lbl_Tc_x = 930;
labelpoints(lbl_ts_x,lbl_Ts_y,lbl_Ts, 'FontSize',11, 'FontWeight','bold');
labelpoints(lbl_Th_x,lbl_Th_y,lbl_Th, 'Interpreter','tex', 'FontSize',10);
labelpoints(lbl_Tc_x,lbl_Tc_y,lbl_Tc, 'Interpreter','tex', 'FontSize',10);

% Save figure
saveas(fig_ORC, [outdir,'/fig7_ORC_Ts.fig']);
saveas(fig_ORC, [outdir,'/fig7_ORC_Ts.png'], 'png');
saveas(fig_ORC, [outdir,'/fig7_ORC_Ts.svg'], 'svg');

%% HX sizing data
% Heat exchanger data
HX.cpmnt = {'';'Preheater';  ...
			'Heat rejector'; '';  ...
			'ORC heat exchanger'; '';  ...
			'ORC condenser'; ''};
HX.config = {'';'2, 4';   ...		% Preheater
		 	'1';	  ...		% Heat rejector (no preheat)
		 	'2';      ...		% Heat rejector (preheat)
		 	'3';	  ...		% ORC heat exchanger (no preheat)
		 	'4';      ...		% ORC heat exchanger (preheat)
		 	'3';      ...		% ORC Condenser (no preheat)
		 	'4'};				% ORC Condenser (preheat)
HX.As_strt = {'Straight'};
HX.As_strt(2:8,1) = num2cell(round(...
		[strt.HX_ph.As;        ...
	 	strt.HX_rj.As;        ...
     	strt.HX_rjPh.As;      ...
	 	strt.HX_ORC.As;       ...
	 	strt.HX_ORCph.As;     ...
	 	strt.HX_cond.As;      ...
	 	strt.HX_condPh.As],   ...
	 	4,'significant'));	% Surface areas (straight) [m^2]
HX.As_sptn(1) = "Serpentine";
HX.As_sptn(2:8,1) = num2cell(round(...
		[sptn.HX_ph.As;        ...
	 	sptn.HX_rj.As;        ...
     	sptn.HX_rjPh.As;      ...
	 	sptn.HX_ORC.As;       ...
	 	sptn.HX_ORCph.As;     ...
	 	sptn.HX_cond.As;      ...
	 	sptn.HX_condPh.As],   ...
	 	4,'significant'));	% Surface areas (serpentine) [m^2]
HX.U_strt(1) = "Straight";
HX.U_strt(2:8,1) = num2cell(round(...
		[strt.HX_ph.U;        ...
	 	strt.HX_rj.U;        ...
     	strt.HX_rjPh.U;      ...
	 	strt.HX_ORC.U;       ...
	 	strt.HX_ORCph.U;     ...
	 	strt.HX_cond.U;      ...
	 	strt.HX_condPh.U],    ...
	 	4,'significant'));	% Overall heat transf coeff. (straight) [W/(m^2*K)]
HX.U_sptn(1) = "Serpentine";
HX.U_sptn(2:8,1) = num2cell(round(...
		[sptn.HX_ph.U;        ...
	 	sptn.HX_rj.U;        ...
     	sptn.HX_rjPh.U;      ...
	 	sptn.HX_ORC.U;       ...
	 	sptn.HX_ORCph.U;     ...
	 	sptn.HX_cond.U;      ...
	 	sptn.HX_condPh.U],    ...
	 	4,'significant'));	% Overall heat transf coeff. (serpentine) [W/(m^2*K)]
% Create cell fields
Tb_HX = table(HX.cpmnt, HX.config, HX.As_strt, HX.As_sptn, HX.U_strt, HX.U_sptn); ...
Tb_HX.Properties.VariableNames = {'Component','Configuration', 'Surface area [m2]', '.', 'U [W m2 K-1)]', '..'};
%Tb_HX.Properties.VariableUnits = {'', '', 'm^2', 'W/(m^2*K)'};
% Write data to file
writetable(Tb_HX, [outdir,'/HXdata.txt'], ...
	'WriteVariableNames', true, 'Delimiter','tab');
writetable(Tb_HX, [outdir,'/HXdata.xlsx'], ...
	'WriteVariableNames', true);
type([outdir,'/HXdata.txt'])

%%
% Save workspace variables to file
save([outdir,'/output.mat'], 'strt', 'sptn', 'y_*', 'y1_*', 'y2_*');

% Reload workspace file to clear clutter
clearvars -except outdir;
load([outdir,'/output.mat']);


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function runSim(TMSch,outdir)	
	% Run parameters script & save workspace (except 'c','mode_' & 'outdir' variables)
	run(['parameters_',TMSch]);
	save([outdir,'/output_',TMSch,'.mat'], '-regexp', '^(?!(c|b_x|mode_[\w|\d]*|outdir[\w|\d]*)$).');
	
	% Run simulations
	disp(['Simulating TMS1 (',TMSch,')...']);
	sim(['TMS/',TMSch,'_TMS1_noPH_noORC.slx']);			% 1) No preheat; no ORC
	disp(['Simulating TMS2 (',TMSch,')...']);
	sim(['TMS/',TMSch,'_TMS2_PH_noORC.slx']);			% 2) Preheat; no ORC
	disp(['Simulating TMS3 (',TMSch,')...']);
	sim(['TMS/',TMSch,'_TMS3_noPH_ORC.slx']);			% 3) No preheat; ORC
	disp(['Simulating TMS4 (',TMSch,')...']);
	sim(['TMS/',TMSch,'_TMS4_PH_ORC.slx']);				% 4) Preheat; ORC
	% Save output data for testing (if needed)
	save(['output_',TMSch,'.mat'], '-regexp', '^(?!(c|b_x|mode_[\w|\d]*|outdir[\w|\d]*)$).');
end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
function plotBarData(outdir,y,c,yLabel,plottitle,bartype,yLim,yScale,name,l, b_h, b_ar)
	% x-axis labels
	x = categorical({'Straight channels', 'Serpentine channels'});
	x = reordercats(x, {'Straight channels', 'Serpentine channels'});
	% Create tiled figure
	fig = figure('Position',[300,100, b_h*b_ar, b_h]);
	
	% Plot data
	if strcmp(bartype,'stacked')
		b = bar(x, y, bartype, 'FaceColor', 'flat'); %, 'BarWidth', 0.4);
	else
		b = bar(x, y, bartype, 'FaceColor', 'flat', 'BarWidth', 0.9);
	end
	set(gca, 'Ticklength', [0 0])
	xlabel('Configuration');
	ylabel(yLabel);
	if strcmp(yScale,'log')
		% Set log scale (if requested)
		set('YScale','log')
	end
	if isempty(yLim) == 0
		% Set y-axis limits (if requested)
		ylim(yLim)
	end
	if  isempty(plottitle) == 0
		% Set title (if requested)
		title(plottitle, 'FontWeight','normal')
	end
	for i = 1:length(y(1,:))
		% Show values on charts
		xtips = b(i).XEndPoints;
		ytips = b(i).YEndPoints;
		labels = string(b(i).YData);
		text(xtips,ytips,labels,'HorizontalAlignment','center',...
    		'VerticalAlignment','bottom')
		% Generate colour data
		b(i).CData = c(i,:);
	end
	% Create legend labels
	if isempty(l) == 0
		lg = legend(b, l, 'Location','northeast');
		%lg.Layout.Tile = 'east';
	end

	% Save figure
	saveas(fig, [outdir,'/fig',name,'.fig']);
	saveas(fig, [outdir,'/fig',name,'.png'], 'png');
	saveas(fig, [outdir,'/fig',name,'.svg'], 'svg');
end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
function plotBarData2(outdir,y1,y2,c,yLabel,subtitle1,subtitle2,bartype,yLim, ...
					  yScale,name,l, b_h, b_ar)
	% x-axis labels
	x = categorical({'Config 1','Config 2','Config 3','Config 4'});
	x = reordercats(x, {'Config 1','Config 2','Config 3','Config 4'});

	% Create tiled figure
	fig = figure('Position',[300,100, b_h*b_ar, b_h]);
	t   = tiledlayout(fig,1,2, 'TileSpacing','compact', 'Padding','tight');
	%t.XLabel.String = 'Configuration';
	t.YLabel.String = yLabel;
	
	% Plot Straight channel data
	t1 = nexttile;
	if strcmp(bartype,'stacked')
		b1 = bar(t1, x, y1, bartype, 'FaceColor', 'flat', 'BarWidth', 0.8);
	else
		b1 = bar(t1, x, y1, bartype, 'FaceColor', 'flat', 'BarWidth', 0.9);
	end
	if isempty(subtitle1) == 0
		title(t1,'Straight Channels',subtitle1, 'Interpreter','tex')
	else
		title(t1,'Straight Channels')
	end
	t1.XAxis.TickLength = [0 0];
	if strcmp(yScale,'log')
		set(t1,'YScale','log')
	end
	if isempty(yLim) == 0
		ylim(t1, yLim)
	end
	for i = 1:length(y1(1,:))
		% Show values on charts
		xtips1 = b1(i).XEndPoints;
		ytips1 = b1(i).YEndPoints;
		labels1 = string(b1(i).YData);
		text(xtips1,ytips1,labels1,'HorizontalAlignment','center',...
    		'VerticalAlignment','bottom')
% 		labels1 = string(abs(b1(i).YData));
% 		text(xtips1,ytips1,[labels1, ' %'],'HorizontalAlignment','center',...
%     		'VerticalAlignment','bottom')
		% Generate colour data
		b1(i).CData = c(i,:);
	end
	
	% Plot Serpentine channel data
	t2 = nexttile;
	if strcmp(bartype,'stacked')
		b2 = bar(t2, x, y2, bartype, 'FaceColor', 'flat', 'BarWidth', 0.8);
	else
		b2 = bar(t2, x, y2, bartype, 'FaceColor', 'flat', 'BarWidth', 0.9);
	end
	if isempty(subtitle2) == 0
		title(t2,'Serpentine Channels',subtitle2, 'Interpreter','tex')
	else
		title(t2,'Serpentine Channels')
	end
	t2.XAxis.TickLength = [0 0];
	if strcmp(yScale,'log')
		set(t2,'YScale','log')
	end
	if isempty(yLim) == 0
		ylim(t2, yLim)
	end
	for i = 1:length(y2(1,:))
		% Show values on charts
		xtips2 = b2(i).XEndPoints;
		ytips2 = b2(i).YEndPoints;
		labels2 = string(b2(i).YData);
		text(xtips2,ytips2,labels2,'HorizontalAlignment','center',...
    		'VerticalAlignment','bottom')
% 		labels2 = string(abs(b2(i).YData));
% 		text(xtips2,ytips2,[labels2, ' %'],'HorizontalAlignment','center',...
%     		'VerticalAlignment','bottom')
		% Generate colour data
		b2(i).CData = c(i,:);
	end
	% Create legend labels
	if isempty(l) == 0
		lg = legend(t1, l, 'Location','northwest');
% 		lg = legend(l);
% 		lg.Layout.Tile = 'north';
	end

	% Save figure
	saveas(fig, [outdir,'/fig',name,'.fig']);
	saveas(fig, [outdir,'/fig',name,'.png'], 'png');
	saveas(fig, [outdir,'/fig',name,'.svg'], 'svg');
end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
function plotPieData(outdir,y1,y2,y3,y4,c,dp,TMSch,subtitle,tilesubtitle,name,l, p_h, p_ar)
	fig = figure('Position',[300,100, p_h*p_ar, p_h]);
	t   = tiledlayout(fig,2,2, 'TileSpacing','compact', 'Padding','tight');
	if isempty(subtitle) == 0
		ttl = title(t,['\bf',TMSch],['\rm',subtitle], 'Interpreter','tex');
		%ttl.TitleFontWeight = 'bold';
		%ttl.SubtitleFontWeight = 'normal';
	else
		title(t,TMSch, 'FontWeight', 'bold')
	end

	% Plot charts
	% Straight
	t1 = nexttile;
	p1 = pie(t1,y1, ['%.',num2str(dp),'f%%']);
	ax = gca();
	ax.Colormap = c;

	t2 = nexttile;
	p2 = pie(t2,y2, ['%.',num2str(dp),'f%%']);
	ax = gca();
	ax.Colormap = c;

	t3 = nexttile;
	p3 = pie(t3,y3, ['%.',num2str(dp),'f%%']);
	ax = gca();
	ax.Colormap = c;

	t4 = nexttile;
	p4 = pie(t4,y4, ['%.',num2str(dp),'f%%']);
	ax = gca();
	ax.Colormap = c;

	% Create title subtitles (if requested)
	if isempty(tilesubtitle) == 0
		title(t1,'Config 1',tilesubtitle(1));
		title(t2,'Config 2',tilesubtitle(2));
		title(t3,'Config 3',tilesubtitle(3));
		title(t4,'Config 4',tilesubtitle(4));
	else
		title(t1,'Config 1');
		title(t2,'Config 2');
		title(t3,'Config 3');
		title(t4,'Config 4');
	end

	% Create legend labels
	if isempty(l) == 0
		lg = legend(l);
		lg.Layout.Tile = 'north';
	end

	% rm 0% labels
	txt1 = findobj(p1,'Type','Text');
	txt2 = findobj(p2,'Type','Text');
	txt3 = findobj(p3,'Type','Text');
	txt4 = findobj(p4,'Type','Text');
	Opc1 = startsWith({txt1.String}, ["0","<"]);
	Opc2 = startsWith({txt2.String}, ["0","<"]);
	Opc3 = startsWith({txt3.String}, ["0","<"]);
	Opc4 = startsWith({txt4.String}, ["0","<"]);
	delete(txt1(Opc1))
	delete(txt2(Opc2))
	delete(txt3(Opc3))
	delete(txt4(Opc4))

	% Save figure
	saveas(fig, [outdir,'/fig',name,'.fig']);
	saveas(fig, [outdir,'/fig',name,'.png'], 'png');
	saveas(fig, [outdir,'/fig',name,'.svg'], 'svg');
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
function plotPieData2(outdir,y1,y2,y3,y4,c,dp,tilesubtitle,name,l, p_h, p_ar)
	fig = figure('Position',[300,100, p_h*p_ar, p_h]);
	t   = tiledlayout(fig,2,4, 'TileSpacing','compact', 'Padding','tight');
	
	% Plot charts
	% Straight
	t11 = nexttile;
	p11 = pie(t11,y1(1,:), ['%.',num2str(dp),'f%%']);
	ax = gca();
	ax.Colormap = c;

	t12 = nexttile;
	p12 = pie(t12,y2(1,:), ['%.',num2str(dp),'f%%']);
	ax = gca();
	ax.Colormap = c;

	t13 = nexttile;
	p13 = pie(t13,y3(1,:), ['%.',num2str(dp),'f%%']);
	ax = gca();
	ax.Colormap = c;

	t14 = nexttile;
	p14 = pie(t14,y4(1,:), ['%.',num2str(dp),'f%%']);
	ax = gca();
	ax.Colormap = c;

	% Serpentine
	t21 = nexttile;
	p21 = pie(t21,y1(2,:), ['%.',num2str(dp),'f%%']);
	ax = gca();
	ax.Colormap = c;

	t22 = nexttile;
	p22 = pie(t22,y2(2,:), ['%.',num2str(dp),'f%%']);
	ax = gca();
	ax.Colormap = c;

	t23 = nexttile;
	p23 = pie(t23,y3(2,:), ['%.',num2str(dp),'f%%']);
	ax = gca();
	ax.Colormap = c;

	t24 = nexttile;
	p24 = pie(t24,y4(2,:), ['%.',num2str(dp),'f%%']);
	ax = gca();
	ax.Colormap = c;

	% Create legend labels
	if isempty(l) == 0
		lg = legend(l);
		lg.Layout.Tile = 'north';
	end

	% Create title subtitles (if requested)
	if isempty(tilesubtitle) == 0
		title(t11,'Config 1, Straight',tilesubtitle(1,1));
		title(t12,'Config 2, Straight',tilesubtitle(1,2));
		title(t13,'Config 3, Straight',tilesubtitle(1,3));
		title(t14,'Config 4, Straight',tilesubtitle(1,4));
		title(t21,'Config 1, Serpentine',tilesubtitle(2,1));
		title(t22,'Config 2, Serpentine',tilesubtitle(2,2));
		title(t23,'Config 3, Serpentine',tilesubtitle(2,3));
		title(t24,'Config 4, Serpentine',tilesubtitle(2,4));
	else
		title(t11,'Config 1, Straight');
		title(t12,'Config 2, Straight');
		title(t13,'Config 3, Straight');
		title(t14,'Config 4, Straight');
		title(t21,'Config 1, Serpentine');
		title(t22,'Config 2, Serpentine');
		title(t23,'Config 3, Serpentine');
		title(t24,'Config 4, Serpentine');
	end

	% rm 0% labels
	txt11 = findobj(p11,'Type','Text');
	txt12 = findobj(p12,'Type','Text');
	txt13 = findobj(p13,'Type','Text');
	txt14 = findobj(p14,'Type','Text');
	txt21 = findobj(p21,'Type','Text');
	txt22 = findobj(p22,'Type','Text');
	txt23 = findobj(p23,'Type','Text');
	txt24 = findobj(p24,'Type','Text');
	Opc11 = startsWith({txt11.String}, ["0","<"]);
	Opc12 = startsWith({txt12.String}, ["0","<"]);
	Opc13 = startsWith({txt13.String}, ["0","<"]);
	Opc14 = startsWith({txt14.String}, ["0","<"]);
	Opc21 = startsWith({txt21.String}, ["0","<"]);
	Opc22 = startsWith({txt22.String}, ["0","<"]);
	Opc23 = startsWith({txt23.String}, ["0","<"]);
	Opc24 = startsWith({txt24.String}, ["0","<"]);
	delete(txt11(Opc11))
	delete(txt12(Opc12))
	delete(txt13(Opc13))
	delete(txt14(Opc14))
	delete(txt21(Opc21))
	delete(txt22(Opc22))
	delete(txt23(Opc23))
	delete(txt24(Opc24))

% 	set(txt11(Opc11),'String', '')
% 	set(txt12(Opc12),'String', '')
% 	set(txt14(Opc14),'String', '')
% 	set(txt21(Opc21),'String', '')
% 	set(txt22(Opc22),'String', '')
% 	set(txt24(Opc24),'String', '')

	% Save figure
	saveas(fig, [outdir,'/fig',name,'.fig']);
	saveas(fig, [outdir,'/fig',name,'.png'], 'png');
	saveas(fig, [outdir,'/fig',name,'.svg'], 'svg');
end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
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