% Main data-gathering script for TMS models
%
% SIMULINK MODEL OUTPUTS:
%	
%	outTMS(i)_pwr   [stk, BoP, tot]				Stack, Total BoP, Total [kW]
%	outTMS(i)_BoP   [pmpFW, pmpClnt, htr]		H2O pump, Clnt pump, heater [W]
%	outTMS(i)_Qclnt [stkIn, PHOut, rjOut]		Stack in, Phtr out, RJ out [kW]
%	outTMS(i)_eff   [stk, BoP, tot]				Stack, BoP, Overall [%]
%   outTMS(i)_ORC   [P_Grs, P_net, Qin, eff]	
%% Set program modes
mode_test = 1;	% Test mode (don't run simulations if != 0)

%%
clearvars('-except', 'mode_*')
close all

% Make output directory
outdir = ['output_', datestr(now,'yyyy-mm-dd--HH-MM-SS')];
mkdir(outdir);


%% Run simulations
% 1) STRAIGHT CHANNELS
outputData(mode_test, 'straight', outdir);

% Clear for next run
clearvars('-except', 'mode_*', 'outdir')

% 2) SERPENTINE CHANNELS
outputData(mode_test, 'serpentine', outdir);


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
function outputData(mode_test, TMSch, outdir)
	% Figure dimensions
	b_h  = 720;  % Height [px]
	b_ar = 1.6; % Aspect ratio [w:h]
	
	% Bar chart colours (RGB)
	c.red    = [0.6350 0.0780 0.1840];
	c.blue   = [0.0000 0.4470 0.7410];
	c.yellow = [0.9290 0.6940 0.1250];
	c.green  = [0.4660 0.6740 0.1880];
	c.orange = [0.8500 0.3250 0.0980];
	c.cyan   = [0.3010 0.7450 0.9330];
	c.purple = [0.4940 0.1840 0.5560];

	% Make output directory for current channels data
	outdir_ch = [outdir, '/',TMSch];
	mkdir(outdir_ch);

	% Load existing output data if testing (skip running simulations)
	if mode_test ~= 0 && exist(['output_',TMSch,'.mat'], 'file')
		load(['output_',TMSch,'.mat']);
	end

	% Run parameters script & save workspace (except 'c','mode_' & 'outdir' variables)
	run(['parameters_',TMSch]);
	save([outdir_ch,'/output.mat'], '-regexp', '^(?!(c|b_x|mode_[\w|\d]*|outdir[\w|\d]*)$).');

	if mode_test == 0
		% Run simulations
		disp(['Simulating TMS1 (',TMSch,')...']);
		sim(['TMS/',TMSch,'_TMS1_noPH_noORC.slx']);			% 1) No preheat; no ORC
		disp(['Simulating TMS2 (',TMSch,')...']);
		sim(['TMS/',TMSch,'_TMS2_PH_noORC.slx']);			% 2) Preheat; no ORC
		disp(['Simulating TMS3 (',TMSch,')...']);
		sim(['TMS/',TMSch,'_TMS3_noPH_ORC.slx']);			% 3) No preheat; ORC
		disp(['Simulating TMS4 (',TMSch,')...']);
		sim(['TMS/',TMSch,'_TMS4_PH_ORC.slx']);				% 4) Preheat; ORC
		disp(['Simulating TMS1A (',TMSch,')...']);
		sim(['TMS/',TMSch,'_TMS1A_noPH_noORC_CL.slx']);		% 1A) No preheat; no ORC (closed FW loop)
		disp(['Simulating TMS2A (',TMSch,')...']);
		sim(['TMS/',TMSch,'_TMS2A_PH_noORC_CL.slx']);		% 2A) Preheat; no ORC (closed FW loop)
		disp(['Simulating TMS3A (',TMSch,')...']);
		sim(['TMS/',TMSch,'_TMS3A_noPH_ORC_CL.slx']);		% 3A) No preheat; ORC (closed FW loop)
		disp(['Simulating TMS4A (',TMSch,')...']);
		sim(['TMS/',TMSch,'_TMS4A_PH_ORC_CL.slx']);			% 4A) Preheat; ORC (closed FW loop)
		
		% Save output data for testing (if needed)
		save(['output_',TMSch,'.mat'], '-regexp', '^(?!(c|b_x|mode_[\w|\d]*|outdir[\w|\d]*)$).');
		
		disp('Simulations complete')
	end
	
	
	%% Power Consumption data [kW]
	% Round data to 3 sf
	y1_pwr  = round([outTMS1_pwr; outTMS2_pwr; outTMS3_pwr; outTMS4_pwr], ...
		3, 'significant');
	y2_pwr =  round([outTMS1A_pwr; outTMS2A_pwr; outTMS3A_pwr; outTMS4A_pwr], ...
		3, 'significant');
	% Bar colours and respective legend labels
	c_pwr = [c.green; c.yellow; c.orange];
	l     = {'Stack', 'BoP', 'Total'};
	% Plot data
	plotData(outdir_ch, y1_pwr, y2_pwr, c_pwr, 'Power consumption [kW]', ...
		'grouped', [], [], '1_pwr', l, b_h, b_ar);
	
	
	%% BoP Power Consumption data [W]
	% Round data to 3 sf
	y1_BoP = round([outTMS1_BoP; outTMS2_BoP; outTMS3_BoP; outTMS4_BoP], ...
		3, 'significant');
	y2_BoP = round([outTMS1A_BoP; outTMS2A_BoP; outTMS3A_BoP; outTMS4A_BoP], ...
		3, 'significant');
	% Bar colours and respective legend labels
	c_BoP = [c.blue; c.cyan; c.orange; c.purple];
	l     = {'Water pump', 'Coolant pump', 'Heater', 'ORC pump'};
	% Plot data
	plotData(outdir_ch, y1_BoP, y2_BoP, c_BoP, 'Power consumption [W]', ...
		'grouped', [0 1E5], 'log', '2_BoP', l, b_h, b_ar);
	
	
	%% Coolant Heat Transfer data [kW]
	% Round data to 3 sf
	y1_Qclnt = round([outTMS1_Qclnt; outTMS2_Qclnt; outTMS3_Qclnt; outTMS4_Qclnt], ...
		3, 'significant');
	y2_Qclnt = round([outTMS1A_Qclnt; outTMS2A_Qclnt; outTMS3A_Qclnt; outTMS4A_Qclnt], ...
		3, 'significant');
	% Bar colours and respective legend labels
	c_Qclnt = [c.orange; c.green; c.cyan];
	l       = {'Stack', 'Preheater', 'Rejected'};
	% Plot data
	plotData(outdir_ch, y1_Qclnt, y2_Qclnt, c_Qclnt, 'Heat transferred [kW]', ...
		'stacked', [], [], '3_Qclnt', l, b_h, b_ar);
	
	
	%% Efficiency data [%]
	% Round data to 3 sf
	y1_eff = round([outTMS1_eff; outTMS2_eff; outTMS3_eff; outTMS4_eff], ...
		3, 'significant');
	y2_eff = round([outTMS1A_eff; outTMS2A_eff; outTMS3A_eff; outTMS4A_eff], ...
		3, 'significant');
	% Bar colours and respective legend labels
	c_eff = [c.orange; c.green];
	l     = {'Stack', 'Overall'};
	% Plot data
	plotData(outdir_ch, y1_eff, y2_eff, c_eff, 'Efficiency [%]', ...
		'grouped',[], [], '4_eff', l, b_h, b_ar);
	
	
	%% ORC data
	% Round data to 3 sf
	y1_ORC = round([sort(outTMS3_ORC(1:3)); sort(outTMS3A_ORC(1:3)); ...
		sort(outTMS4_ORC(1:3)); sort(outTMS4A_ORC(1:3))], ...
		3, 'significant');
	y2_ORC = round([outTMS3_ORC(4); outTMS3A_ORC(4); outTMS4_ORC(4); outTMS4A_ORC(4)], ...
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
	saveas(fig, [outdir_ch,'/fig5_ORC.fig']);
	saveas(fig, [outdir_ch,'/fig5_ORC.png'], 'png');
	saveas(fig, [outdir_ch,'/fig5_ORC.svg'], 'svg');
	
	
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
	writetable(Tb_HX, [outdir_ch,'/HXdata.txt'], ...
		'WriteVariableNames', true, 'Delimiter','tab');
	writetable(Tb_HX, [outdir_ch,'/HXdata.xlsx'], ...
		'WriteVariableNames', true);
	type([outdir_ch,'/HXdata.txt'])
	
	%%
	% Save workspace variables to file
	save([outdir_ch,'/output.mat'], 'outTMS*', 'y1_*', 'y2_*', '-append');
	
	% Reload workspace file to clear clutter
	clearvars -except outdir outdir_ch;
	load([outdir_ch,'/output.mat']);
end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plotData(outdir_ch,y1,y2,c,yLabel,bartype,yLim,yScale,name,l, b_h, b_ar)
	% x-axis labels
	x1 = categorical({'1','2','3','4'});
	x1 = reordercats(x1, {'1','2','3','4'});
	x2 = categorical({'1A','2A','3A','4A'});
	x2 = reordercats(x2, {'1A','2A','3A','4A'});

	% Create tiled figure
	fig = figure('Position',[300,250, b_h*b_ar, b_h]);
	t   = tiledlayout(fig,2,1, 'TileSpacing','compact', 'Padding','tight');
	t.XLabel.String = 'Configuration';
	t.YLabel.String = yLabel;
	
	% Plot open FW loop data
	t1 = nexttile;
	if strcmp(bartype,'stacked')
		b1 = bar(t1, x1, y1, bartype, 'FaceColor', 'flat', 'BarWidth', 0.4);
	else
		b1 = bar(t1, x1, y1, bartype, 'FaceColor', 'flat', 'BarWidth', 0.9);
	end
	title(t1,'Open Feedwater Loop')
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
		% Generate colour data
		b1(i).CData = c(i,:);
	end
	
	% Plot closed FW loop data
	t2 = nexttile;
	if strcmp(bartype,'stacked')
		b2 = bar(t2, x2, y2, bartype, 'FaceColor', 'flat', 'BarWidth', 0.4);
	else
		b2 = bar(t2, x2, y2, bartype, 'FaceColor', 'flat', 'BarWidth', 0.9);
	end
	title(t2,'Closed Feedwater Loop')
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
		% Generate colour data
		b2(i).CData = c(i,:);
	end
	% Create legend labels
	lg = legend(t1, l, 'Location','northeast');
	%lg.Layout.Tile = 'east';

	% Save figure
	saveas(fig, [outdir_ch,'/fig',name,'.fig']);
	saveas(fig, [outdir_ch,'/fig',name,'.png'], 'png');
	saveas(fig, [outdir_ch,'/fig',name,'.svg'], 'svg');
end