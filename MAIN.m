% Main data-gathering script for TMS models
%
% SIMULINK MODEL OUTPUTS:
%	
%	outTMS(i)_pwr   [stk, BoP, tot]		   Stack, Total BoP, Total [kW]
%	outTMS(i)_BoP   [pmpFW, pmpClnt, htr]  H2O pump, Clnt pump, heater [W]
%	outTMS(i)_Qclnt [stkIn, PHOut, rjOut]  Stack in, Phtr out, RJ out [kW]
%	outTMS(i)_eff   [stk, BoP, tot]		   Stack, BoP, Overall [%]
%   outTMS(i)_ORC   [
clear all
close all

% Bar chart properties
b_h  = 480;  % Chart figure height [px]
b_ar = 2.5;  % Chart figure aspect ratio [w:h]
%b_x = categorical({'A','B','C','D'});
%b_x = reordercats(b_x, {'A','B','C','D'});
b_x = categorical({'A','B'});
b_x = reordercats(b_x, {'A','B'});
% Bar chart colours (RGB)
c.red    = [0.6350 0.0780 0.1840];
c.blue   = [0.0000 0.4470 0.7410];
c.yellow = [0.9290 0.6940 0.1250];
c.green  = [0.4660 0.6740 0.1880];
c.orange = [0.8500 0.3250 0.0980];
c.cyan   = [0.3010 0.7450 0.9330];
c.purple = [0.4940 0.1840 0.5560];

% Make output directory
outdir = ['output_', datestr(now,'yyyy-mm-dd--HH-MM-SS')];
mkdir(outdir);

% Make output directory for straight channels data
outdir_ch = [outdir, '/straight'];
mkdir(outdir_ch);

%% 1) Plot data
parameters_SS;  % Run parameters script

% Run simulations
sim("TMS1_noPH_noORC.slx");	% No preheat; no ORC
sim("TMS2_PH_noORC.slx");		% Preheat; no ORC
%sim("TMS3_noPH_ORC.slx");		% No preheat; ORC
%sim("TMS4_PH_ORC.slx");		% Preheat; ORC

% Plot Power Consumption data [kW]
y1_pwr = round([outTMS1_pwr; outTMS2_pwr], 3, 'significant');  % round data to 3 sf
fig11 = figure('Position',[300,300, b_h*b_ar, b_h]);
b1_pwr = bar(b_x, y1_pwr, 'grouped', 'FaceColor', 'flat');
hold on
ylabel('Power consumption [kW]')
%set(gca,'YScale','log')
c11 = [c.green; c.yellow; c.orange];  % [Stack; BoP; Total] colours
	
for i = 1:length(y1_pwr(1,:))
	% Show values on chart
	xtips1 = b1_pwr(i).XEndPoints;
	ytips1 = b1_pwr(i).YEndPoints;
	labels1 = string(b1_pwr(i).YData);
	text(xtips1,ytips1,labels1,'HorizontalAlignment','center',...
    	'VerticalAlignment','bottom')
	% Generate colour data
	b1_pwr(i).CData = c11(i,:);
end
% Create legend labels
L11 = cell(1,length(y1_pwr(1,:)));
L11{1} = 'Stack';
L11{2} = 'BoP';
L11{3} = 'Total';
legend(b1_pwr, L11, 'Location','northeast');
% Save figure
saveas(fig11, [outdir_ch,'/fig1_pwr.fig']);
saveas(fig11, [outdir_ch,'/fig1_pwr.png'], 'png');
saveas(fig11, [outdir_ch,'/fig1_pwr.svg'], 'svg');


% Plot BoP Power Consumption data [W]
y1_BoP = round([outTMS1_BoP;   outTMS2_BoP], 3, 'significant');  % round data to 3 sf
fig12 = figure('Position',[300,300, b_h*b_ar, b_h]);
b1_BoP = bar(b_x, y1_BoP, 'grouped', 'FaceColor', 'flat');
hold on
ylabel('Power consumption [W]')
set(gca,'YScale','log')
c12 = [c.blue; c.cyan;c.orange];  % [H2O pump; Clnt pump; Heater] colours
for i = 1:length(y1_BoP(1,:))
	% Show values on chart
	xtips1 = b1_BoP(i).XEndPoints;
	ytips1 = b1_BoP(i).YEndPoints;
	labels1 = string(b1_BoP(i).YData);
	text(xtips1,ytips1,labels1,'HorizontalAlignment','center',...
    	'VerticalAlignment','bottom')
	% Generate colour data
	b1_BoP(i).CData = c12(i,:);
end
% Create legend labels
L12 = cell(1,length(y1_BoP(1,:)));
L12{1} = 'Water pump';
L12{2} = 'Coolant pump';
L12{3} = 'Heater';
legend(b1_BoP, L12, 'Location','northeast');
% Save figure
saveas(fig12, [outdir_ch,'/fig2_BoP.fig']);
saveas(fig12, [outdir_ch,'/fig2_BoP.png'], 'png');
saveas(fig12, [outdir_ch,'/fig2_BoP.svg'], 'svg');


% Plot Coolant Heat Transfer data [kW]
y1_Qclnt = round([outTMS1_Qclnt; outTMS2_Qclnt], 3, 'significant');  % round data to 3 sf
fig13 = figure('Position',[300,300, b_h*b_ar, b_h]);
b1_Qclnt = bar(b_x, y1_Qclnt, 'stacked', 'FaceColor', 'flat');
hold on
ylabel('Heat transferred [kW]')
%set(gca,'YScale','log')
c13 = [c.orange; c.blue; c.cyan];  % [Stack; Preheater; Rejector] colours
for i = 1:length(y1_Qclnt(1,:))
	% Show values on chart
	xtips1 = b1_Qclnt(i).XEndPoints;
	ytips1 = b1_Qclnt(i).YEndPoints;
	labels1 = string(b1_Qclnt(i).YData);
	text(xtips1,ytips1,labels1,'HorizontalAlignment','center',...
    	'VerticalAlignment','bottom')
	% Generate colour data
	b1_Qclnt(i).CData = c13(i,:);
end
% Create legend labels
L12 = cell(1,length(y1_Qclnt(1,:)));
L12{1} = 'Stack';
L12{2} = 'Prehater';
L12{3} = 'Rejected';
%L12{1} = 'PHeater out';
%L12{2} = 'Rejector out';
legend(b1_Qclnt, L12, 'Location','northeast');
% Save figure
saveas(fig13, [outdir_ch,'/fig3_Qclnt.fig']);
saveas(fig13, [outdir_ch,'/fig3_Qclnt.png'], 'png');
saveas(fig13, [outdir_ch,'/fig3_Qclnt.svg'], 'svg');

% % Plot Coolant Heat Transfer data [kW]
% y1_Qclnt = round(abs([outTMS1_Qclnt(2:end); outTMS2_Qclnt(2:end)]), 3, 'significant');  % round data to 3 sf
% fig13 = figure('Position',[300,300, b_h*b_ar, b_h]);
% b1_Qclnt = bar(b_x, y1_Qclnt, 'stacked', 'FaceColor', 'flat');  % Plot
% hold on
% ylabel('Heat transferred [kW]')
% %set(gca,'YScale','log')
% c13 = [c.blue; c.cyan];			% [Preheater; Rejector] colours
% for i = 1:length(y1_Qclnt(1,:))
% 	% Show values on chart
% 	xtips1 = b1_Qclnt(i).XEndPoints;
% 	ytips1 = b1_Qclnt(i).YEndPoints;
% 	labels1 = string(b1_Qclnt(i).YData);
% 	text(xtips1,ytips1,labels1,'HorizontalAlignment','center',...
%     	'VerticalAlignment','bottom')
% 	% Generate colour data
% 	b1_Qclnt(i).CData = c13(i,:);
% end
% % Plot line of total heat in
% y1_QclntS = round(outTMS1_Qclnt(1));
% b1_QclntS = yline(y1_QclntS,'--k', [num2str(y1_QclntS), ' kW']);
% b1_QclntS.LabelHorizontalAlignment = 'center';
% % Create legend labels
% L12 = cell(1,length(y1_Qclnt(1,:))-1);
% L12{1} = 'PHeater';
% L12{2} = 'Rejector';
% legend([b1_QclntS b1_Qclnt], ['Stack in' L12], 'Location','southeast');
% 
% 
% % Plot Coolant Heat Transfer data [kW]
% y1_Qclnt = round(abs([outTMS1_Qclnt(2:end); outTMS2_Qclnt(2:end)]), 3, 'significant');  % round data to 3 sf
% fig13 = figure('Position',[300,300, b_h*b_ar, b_h]);
% b1_Qclnt = bar(b_x, y1_Qclnt, 'stacked', 'FaceColor', 'flat');  % Plot
% hold on
% ylabel('Heat transferred [kW]')
% %set(gca,'YScale','log')
% c13 = [c.blue; c.cyan];			% [Preheater; Rejector] colours
% for i = 1:length(y1_Qclnt(1,:))
% 	% Show values on chart
% 	xtips1 = b1_Qclnt(i).XEndPoints;
% 	ytips1 = b1_Qclnt(i).YEndPoints;
% 	labels1 = string(b1_Qclnt(i).YData);
% 	text(xtips1,ytips1,labels1,'HorizontalAlignment','center',...
%     	'VerticalAlignment','bottom')
% 	% Generate colour data
% 	b1_Qclnt(i).CData = c13(i,:);
% end
% % Plot line of total heat in
% y1_QclntS = round(outTMS1_Qclnt(1));
% b1_QclntS = yline(y1_QclntS,'--k', ['$\mathrm{\dot{Q}_{clnt,stk}}$ = ', num2str(y1_QclntS), ' kW'],...
% 	'Interpreter', 'latex');
% b1_QclntS.LabelHorizontalAlignment = 'center';
% % Create legend labels
% L12 = cell(1,length(y1_Qclnt(1,:))-1);
% L12{1} = 'PHeater';
% L12{2} = 'Rejector';
% legend(b1_Qclnt,  L12, 'Location','southeast');


% Plot efficiency data [%]
y1_eff = round([outTMS1_eff; outTMS2_eff], 3, 'significant');  % round data to 3 sf
fig14 = figure('Position',[300,300, b_h*b_ar, b_h]);
b1_eff = bar(b_x, y1_eff, 'grouped', 'FaceColor', 'flat');
hold on
ylabel('Efficiency [%]')
%set(gca,'YScale','log')
c14 = [c.orange; c.yellow; c.green];  % [Stack; BoP; Overall] colours
for i = 1:length(y1_eff(1,:))
	% Show values on chart
	xtips1 = b1_eff(i).XEndPoints;
	ytips1 = b1_eff(i).YEndPoints;
	labels1 = string(b1_eff(i).YData);
	text(xtips1,ytips1,labels1,'HorizontalAlignment','center',...
    	'VerticalAlignment','bottom')
	% Generate colour data
	b1_eff(i).CData = c14(i,:);
end
% Create legend labels
L11 = cell(1,length(y1_eff(1,:)));
L11{1} = 'Stack';
L11{2} = 'BoP';
L11{3} = 'Overall';
legend(b1_eff, L11, 'Location','northeast');
% Save figure
saveas(fig14, [outdir_ch,'/fig4_eff.fig']);
saveas(fig14, [outdir_ch,'/fig4_eff.png'], 'png');
saveas(fig14, [outdir_ch,'/fig4_eff.svg'], 'svg');

%% 2) Record HX sizing data
HX_nvars = 2;	% No. heat exchanger variables
HX_n     = 6;	% No. heat exchangers
% Heat exchanger data
HX.As = round(...
		[HX_ph.As;   ...
		HX_rj.As;    ...
        HX_rjPH.As;  ...
		HX_ORC.As;   ...
		HX_ORCph.As; ...
		HX_cond.As], ...
		4,'significant');	% Surface areas [m^2]
HX.U = round(...
	   [HX_ph.U;    ...
		HX_rj.U;    ...
		HX_rjPH.U;  ...
		HX_ORC.U;   ...
		HX_ORCph.U; ...
		HX_cond.U], ...
		4,'significant');	% Overall heat transf coeff. [W/(m^2*K)]
% Create cell fields
Tb_HX = table(HX.As, HX.U); ...
Tb_HX.Properties.VariableNames = {'Surface area [m^2]', 'U [W/(m^2*K)]'};
Tb_HX.Properties.VariableUnits = {'m^2', 'W/(m^2*K)'};
Tb_HX.Properties.RowNames = {'Preheater'; ...
							 'Heat exchanger (config A)'; ...
							 'Heat exchanger (config B)'; ...
							 'Heat exchanger (config C)'; ...
							 'Heat exchanger (config D)'; ...
							 'ORC condenser'};
% Write data to file
writetable(Tb_HX, [outdir_ch,'/HXdata.txt'], ...
	'WriteRowNames', true, 'WriteVariableNames', true, 'Delimiter','tab');
type([outdir_ch,'/HXdata.txt'])

% Save workspace variables to file
%save([outdir_ch,'/output.mat']);