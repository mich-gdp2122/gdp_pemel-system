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
mode_test = 1;	% Test mode (don't run simulations if != 0)

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
% Round data to 3 sf
y_totBoP  = round([strt.TMS1_pwr(2), strt.TMS2_pwr(2), strt.TMS3_pwr(2), strt.TMS4_pwr(2); ...
				   sptn.TMS1_pwr(2), sptn.TMS2_pwr(2), sptn.TMS3_pwr(2), sptn.TMS4_pwr(2)], ...
			3, 'significant');
% Bar colours and respective legend labels
c_totBoP = [c.orange; c.blue; c.yellow; c.cyan];
l     = {'Config 1', 'Config 2', 'Config 3', 'Config 4'};
% Plot data
plotBarData(outdir, y_totBoP, c_totBoP, 'Power consumption [kW]', ...
	'grouped', [], [], '1_totBoP', l, b_h, b_ar);


%% BoP Power Consumption data [% tot]
% Round data to 3 sf
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
% Bar colours and respective legend labels
c_BoP = [c.blue; c.cyan; c.orange; c.purple];
l     = {'Water pump', 'Coolant pump', 'Heater', 'ORC pump'};
% Plot data
plotPieData(outdir, y1_BoP, y2_BoP, y3_BoP, y4_BoP, c_BoP, '2_BoP', l, b_h, b_ar);


%% Coolant Heat Transfer data [kW]
% Round data to 3 sf
y1_Qclnt = round([TMS1_Qclnt; TMS2_Qclnt; TMS3_Qclnt; TMS4_Qclnt], ...
	3, 'significant');
y2_Qclnt = round([TMS1A_Qclnt; TMS2A_Qclnt; TMS3A_Qclnt; TMS4A_Qclnt], ...
	3, 'significant');
% Bar colours and respective legend labels
c_Qclnt = [c.orange; c.green; c.cyan];
l       = {'Stack', 'Preheater', 'Rejected'};
% Plot data
plotBarData(outdir, y1_Qclnt, y2_Qclnt, c_Qclnt, 'Heat transferred [kW]', ...
	'stacked', [], [], '3_Qclnt', l, b_h, b_ar);


%% Efficiency data [%]
% Round data to 3 sf
y1_eff = round([TMS1_eff; TMS2_eff; TMS3_eff; TMS4_eff], ...
	3, 'significant');
y2_eff = round([TMS1A_eff; TMS2A_eff; TMS3A_eff; TMS4A_eff], ...
	3, 'significant');
% Bar colours and respective legend labels
c_eff = [c.orange; c.green];
l     = {'Stack', 'Overall'};
% Plot data
plotBarData(outdir, y1_eff, y2_eff, c_eff, 'Efficiency [%]', ...
	'grouped',[], [], '4_eff', l, b_h, b_ar);


%% ORC data
% Round data to 3 sf
y1_ORC = round([sort(TMS3_ORC(1:3)); sort(TMS3A_ORC(1:3)); ...
	sort(TMS4_ORC(1:3)); sort(TMS4A_ORC(1:3))], ...
	3, 'significant');
y2_ORC = round([TMS3_ORC(4); TMS3A_ORC(4); TMS4_ORC(4); TMS4A_ORC(4)], ...
	3, 'significant');

% x labels
x_ORC = categorical({'3','3A','4','4A'});
x_ORC = reordercats(x_ORC, {'3','3A','4','4A'});

% Bar colours
c1_ORC = [c.yellow; c.cyan; c.orange];
c2_ORC = c.green;

% Create tiled figure
fig = figure('Position',[300,250, b_h*b_ar, b_h]);
t   = tiledlayout(fig,2,1, 'TileSpacing','compact', 'Padding','tight');
t.XLabel.String = 'Configuration';

% Plot ORC power data
t1 = nexttile;
b1 = bar(t1, x_ORC, y1_ORC, 'stacked', 'FaceColor', 'flat', 'BarWidth', 0.4);
ylabel('Power [kW]')
%title(t1,'Power Balance')
t1.XAxis.TickLength = [0 0];
for i = 1:length(y1_ORC(1,:))
	% Show values on charts
	xtips1 = b1(i).XEndPoints;
	ytips1 = b1(i).YEndPoints;
	labels1 = string(b1(i).YData);
	text(xtips1,ytips1,labels1,'HorizontalAlignment','center',...
		'VerticalAlignment','bottom')
	% Generate colour data
	b1(i).CData = c1_ORC(i,:);
end

% Plot ORC efficiency
t2 = nexttile;
b2 = bar(t2, x_ORC, y2_ORC, 'FaceColor', 'flat', 'BarWidth', 0.4);
ylabel('Efficiency [%]')
ylim(t2,[0, 50])
%title(t2,'Thermal Efficiency')
t2.XAxis.TickLength = [0 0];
for i = 1:length(y2_ORC(1,:))
	% Show values on charts
	xtips2 = b2(i).XEndPoints;
	ytips2 = b2(i).YEndPoints;
	labels2 = string(b2(i).YData);
	text(xtips2,ytips2,labels2,'HorizontalAlignment','center',...
		'VerticalAlignment','bottom')
	% Generate colour data
	b2(i).CData = c2_ORC;
end
% Create legend labels
l  = {'Net pwr', 'Gross pwr', 'Heat in'};
lg = legend(t1, l, 'Location','northeast');
%lg.Layout.Tile = 'east';

% Save figure
saveas(fig, [outdir,'/fig5_ORC.fig']);
saveas(fig, [outdir,'/fig5_ORC.png'], 'png');
saveas(fig, [outdir,'/fig5_ORC.svg'], 'svg');


%% HX sizing data
% Heat exchanger data
HX.cpmnt = {'Preheater'; ''; ...
			'Heat rejector'; ''; ''; ...
			'ORC heat exchanger'; ''; ''; ...
			'ORC condenser'; ''; ''};
HX.config = {'2, 4';   ...		% Preheater (open FW loop)
		 	'2A, 4A'; ...		% Preheater (closed FW loop)
		 	'1, 1A';  ...		% Heat rejector (no preheat)
		 	'2';      ...		% Heat rejector (preheat, open loop)
		 	'2A';     ...		% Heat rejector (preheat, closed loop)
		 	'3, 3A';  ...		% ORC heat exchanger (no preheat)
		 	'4';      ...		% ORC heat exchanger (preheat, open loop)
		 	'4A';     ...		% ORC heat exchanger (preheat, closed loop)
		 	'3, 3A';  ...		% ORC Condenser (no preheat)
		 	'4';      ...		% ORC Condenser (preheat, open loop)
		 	'4A'};				% ORC Condenser (preheat, closed loop)
HX.As = round(...
		[HX_ph.As;        ...
	 	HX_phCL.As;      ...
	 	HX_rj.As;        ...
     	HX_rjPh.As;      ...
	 	HX_rjPhCL.As;    ...
	 	HX_ORC.As;       ...
	 	HX_ORCph.As;     ...
	 	HX_ORCphCL.As;   ...
	 	HX_cond.As;      ...
	 	HX_condPh.As;    ...
	 	HX_condPhCL.As], ...
	 	4,'significant');	% Surface areas [m^2]
HX.U = round(...
		[HX_ph.U;        ...
	 	HX_phCL.U;      ...
	 	HX_rj.U;        ...
     	HX_rjPh.U;      ...
	 	HX_rjPhCL.U;    ...
	 	HX_ORC.U;       ...
	 	HX_ORCph.U;     ...
	 	HX_ORCphCL.U;   ...
	 	HX_cond.U;      ...
	 	HX_condPh.U;    ...
	 	HX_condPhCL.U], ...
	 	4,'significant');	% Overall heat transf coeff. [W/(m^2*K)]
% Create cell fields
Tb_HX = table(HX.cpmnt, HX.config, HX.As, HX.U); ...
Tb_HX.Properties.VariableNames = {'Component','Configuration', 'Surface area [m^2]', 'U [W/(m^2*K)]'};
Tb_HX.Properties.VariableUnits = {'', '', 'm^2', 'W/(m^2*K)'};
% Write data to file
writetable(Tb_HX, [outdir,'/HXdata.txt'], ...
	'WriteVariableNames', true, 'Delimiter','tab');
writetable(Tb_HX, [outdir,'/HXdata.xlsx'], ...
	'WriteVariableNames', true);
type([outdir,'/HXdata.txt'])

%%
% Save workspace variables to file
save([outdir,'/output.mat'], 'TMS*', 'y1_*', 'y2_*', '-append');

% Reload workspace file to clear clutter
clearvars -except outdir outdir_ch;
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
function plotBarData(outdir,y,c,yLabel,bartype,yLim,yScale,name,l, b_h, b_ar)
	% x-axis labels
	x = categorical({'Straight channels', 'Serpentine channels'});
	x = reordercats(x, {'Straight channels', 'Serpentine channels'});
	% Create tiled figure
	fig = figure('Position',[300,250, b_h*b_ar, b_h]);
	
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
function plotPieData(outdir,y1,y2,y3,y4,c,name,l, p_h, p_ar)
	fig = figure('Position',[300,250, p_h*p_ar, p_h]);
	t   = tiledlayout(fig,2,4); %, 'TileSpacing','compact', 'Padding','tight');

	% Plot charts
	% Straight
	t11 = nexttile;
	p11 = pie(t11,y1(1,:), 'explode');
	title(t11,'Config 1');
	pText11 = findobj(p11,'Type','text');	% rm 0% labels
	pText11(4).String = '';

	t12 = nexttile;
	p12 = pie(t12,y2(1,:), 'explode');
	title(t12,'Config 2');
	pText12 = findobj(p12,'Type','text');	% rm 0% labels
	pText12(3).String = '';
	pText12(4).String = '';

	t13 = nexttile;
	p13 = pie(t13,y3(1,:), 'explode');
	title(t13,'Config 3');

	t14 = nexttile;
	p14 = pie(t14,y4(1,:), 'explode');
	title(t14,'Config 4');
	pText14 = findobj(p14,'Type','text');	% rm 0% labels
	pText14(3).String = '';

	% Serpentine
	t21 = nexttile;
	p21 = pie(t21,y1(2,:), 'explode');
	title(t21,'Config 1');
	pText21 = findobj(p21,'Type','text');	% rm 0% labels
	pText21(4).String = '';

	t22 = nexttile;
	p22 = pie(t22,y2(2,:), 'explode');
	title(t22,'Config 2');
	pText22 = findobj(p22,'Type','text');	% rm 0% labels
	pText22(3).String = '';
	pText22(4).String = '';

	t23 = nexttile;
	p23 = pie(t23,y3(2,:), 'explode');
	title(t23,'Config 3');

	t24 = nexttile;
	p24 = pie(t24,y4(2,:), 'explode');
	title(t24,'Config 4');
	pText24 = findobj(p24,'Type','text');	% rm 0% labels
	pText24(3).String = '';


	% Create legend labels
	if isempty(l) == 0
		lg = legend(l);
		lg.Layout.Tile = 'north';
	end

	% Save figure
	saveas(fig, [outdir,'/fig',name,'.fig']);
	saveas(fig, [outdir,'/fig',name,'.png'], 'png');
	saveas(fig, [outdir,'/fig',name,'.svg'], 'svg');
end