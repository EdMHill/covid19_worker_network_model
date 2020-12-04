
%% make final size plots
make_violinplots('adherence','final_size')
make_violinplots('workpercent','final_size')
%make_violinplots('workpercent_and_backwardsCT','final_size')
make_violinplots('worker_patterns','final_size')

%% make peak infections plots
make_violinplots('adherence','peak_inf')
make_violinplots('workpercent','peak_inf')
%make_violinplots('workpercent_and_backwardsCT','peak_inf')
make_violinplots('worker_patterns','peak_inf')

%% total self-isolating days
make_violinplots('adherence','total_isolation')
make_violinplots('workpercent','total_isolation')
%make_violinplots('workpercent_and_backwardsCT','total_isolation')
make_violinplots('worker_patterns','total_isolation')

%% peak fraction self-isolating
make_violinplots('adherence','peak_isolation')
make_violinplots('workpercent','peak_isolation')
%make_violinplots('workpercent_and_backwardsCT','peak_isolation')
make_violinplots('worker_patterns','peak_isolation')

%% outbreak duration
make_violinplots('adherence','duration')
make_violinplots('workpercent','duration')
%make_violinplots('workpercent_and_backwardsCT','duration')
make_violinplots('worker_patterns','duration')

%% Threshold event plots: Adherence
variable_names = {'final_size','peak_inf','avg_isolation','duration','peak_isolation'};
criteria_thresholds = [0.5,0.01,0.01,150,0.05];
display_legend = true;
make_threshold_event_plots('adherence',variable_names,criteria_thresholds,display_legend)

%% Threshold event plots: Workpercent
variable_names = {'final_size','peak_inf','avg_isolation','duration','peak_isolation'};
criteria_thresholds = [0.5,0.01,0.01,150,0.05];
display_legend = true;
make_threshold_event_plots('workpercent',variable_names,criteria_thresholds,display_legend)

%% Threshold event plots: Synch & asynch worker patterns
variable_names = {'final_size','peak_inf','avg_isolation','duration','peak_isolation'};
criteria_thresholds = [0.5,0.01,0.01,150,0.05];
display_legend = true;
make_threshold_event_plots('worker_patterns',variable_names,criteria_thresholds,display_legend)

%% COVID-secure temporal outputs
make_temporal_grid('COVID_secure','inf_prevalence')
make_temporal_grid('COVID_secure','isol_prevalence')
make_temporal_grid('COVID_secure_no_isol','inf_prevalence')
%% COVID-secure heatmap outputs
variable_names = {'final_size','peak_inf','avg_isolation','peak_isolation','duration'};
make_heatmaps('COVID_secure',variable_names)
make_heatmaps('COVID_secure_no_isol',variable_names)

%% temporal plots (workpercent & worker patterns)
clear variables
% Set up plot
if ispc==1
    scale = 0.9/0.75;
else
    scale = 1;
end
position = [10, 10, 3.5*550*scale, 4.5*450*scale];
set(0, 'DefaultFigurePosition', position);
%set(0, 'DefaultFigurePosition', position,'Units','points');
fig = figure();
clf;
set(fig,'Color', [1 1 1])
fontsize=26;
max_y = 0.2;
subplot(3,3,1); hold on; maxy(1) = make_temporal_plots({'workpercent'},'inf_prevalence',1); ylabel('Propn infectious'); title('Working from home'); set(gca,'FontSize',fontsize); xlim([0,200]);
yticks([0 0.05 0.10 0.15 0.2])
ytickformat('%.2f')
subplot(3,3,2); hold on; maxy(2) = make_temporal_plots({'synch'},'inf_prevalence',1); title('Synchronous'); set(gca,'FontSize',fontsize); xlim([0,200]);
yticks([0 0.05 0.10 0.15 0.2])
ytickformat('%.2f')
subplot(3,3,3); hold on; maxy(3) = make_temporal_plots({'asynch'},'inf_prevalence',1); title('Asynchronous'); set(gca,'FontSize',fontsize); xlim([0,200]);
yticks([0 0.05 0.10 0.15 0.2])
ytickformat('%.2f')
subplot(3,3,1); ylim([0,max_y]); subplot(3,3,2); ylim([0,max_y]); subplot(3,3,3); ylim([0,max_y]);

subplot(3,3,4); hold on; maxy(1) = make_temporal_plots({'workpercent'},'isol_prevalence',1); ylabel('Propn isolating'); set(gca,'FontSize',fontsize); xlim([0,200]);
yticks([0 0.05 0.10 0.15 0.20 0.25])
ytickformat('%.2f')
subplot(3,3,5); hold on; maxy(2) = make_temporal_plots({'synch'},'isol_prevalence',1); set(gca,'FontSize',fontsize); xlim([0,200]);
yticks([0 0.05 0.10 0.15 0.20 0.25])
ytickformat('%.2f')
subplot(3,3,6); hold on; maxy(3) = make_temporal_plots({'asynch'},'isol_prevalence',1); set(gca,'FontSize',fontsize); xlim([0,200]);
yticks([0 0.05 0.10 0.15 0.20 0.25])
ytickformat('%.2f')
subplot(3,3,4); ylim([0,max(maxy)]); subplot(3,3,5); ylim([0,max(maxy)]); subplot(3,3,6); ylim([0,max(maxy)]);

subplot(3,3,7); hold on; maxy(1) = make_temporal_plots({'workpercent'},'Rt',1); xlabel('Time (days)'); ylabel('R_t'); set(gca,'FontSize',fontsize); xlim([0,100]); plot([0,100],[1,1],'--','Linewidth',1.5,'Color',[0.5 0.5 0.5]);
yticks([0 1 2 3 4])
subplot(3,3,8); hold on; maxy(2) = make_temporal_plots({'synch'},'Rt',1); xlabel('Time (days)'); set(gca,'FontSize',fontsize); xlim([0,100]); plot([0,100],[1,1],'--','Linewidth',1.5,'Color',[0.5 0.5 0.5]);
yticks([0 1 2 3 4])
subplot(3,3,9); hold on; maxy(3) = make_temporal_plots({'asynch'},'Rt',1); xlabel('Time (days)'); set(gca,'FontSize',fontsize); xlim([0,100]); plot([0,100],[1,1],'--','Linewidth',1.5,'Color',[0.5 0.5 0.5]);
yticks([0 1 2 3 4])
subplot(3,3,7); ylim([0,4]); subplot(3,3,8); ylim([0,4]); subplot(3,3,9); ylim([0,4]);
export_fig('worker_pattern_plots/worker_pattern_temporal_plots','-pdf','-r600')

% num_lines = 6;
% basecolour = [0 0 0];
% light = [0.8 0.8 0.8];
% colour_vec = [linspace(light(1),basecolour(1),num_lines)', linspace(light(2),basecolour(2),num_lines)', linspace(light(3),basecolour(3),num_lines)'];
% legend_label = {'0 days','1 day','2 days','3 days','4 days','5 days'};
% H = gobjects(num_lines,1);
% for line_itr = 1:num_lines
%     H(line_itr) = fill(NaN,NaN,colour_vec(line_itr,:),...
%                                 'DisplayName',legend_label{line_itr});
%                                 %'FaceAlpha',0.25);
% end
% legend(H,...
%     'LineWidth',1.5,...
%     'FontSize',fontsize,...
%     'Position',[0.839220779220779 0.77 0.05 0.10]);

%% Temporal plots (Adherence)
clear variables
% Set up plot
if ispc==1
    scale = 0.9/0.75;
else
    scale = 1;
end
position = [10, 10, 3.5*550*scale, 1.5*450*scale];
set(0, 'DefaultFigurePosition', position);
fig = figure();
clf;
set(fig,'Color', [1 1 1])
fontsize=26;
subplot(1,3,1); hold on; maxy(1) = make_temporal_plots({'adherence'},'inf_prevalence',1); xlabel('Time (days)'); ylabel('Propn infectious'); set(gca,'FontSize',fontsize); xlim([0,200]);
yticks([0 0.05 0.10 0.15 0.20 0.25 0.30 0.35 0.40])
ytickformat('%.2f')
subplot(1,3,2); hold on; maxy(2) = make_temporal_plots({'adherence'},'isol_prevalence',1); xlabel('Time (days)'); ylabel('Propn isolating'); set(gca,'FontSize',fontsize); xlim([0,200]);
yticks([0 0.05 0.10 0.15 0.20 0.25 0.30])
ytickformat('%.2f')
subplot(1,3,3); hold on; maxy(3) = make_temporal_plots({'adherence'},'Rt',1); xlabel('Time (days)'); ylabel('R_t'); set(gca,'FontSize',fontsize); xlim([0,100]); plot([0,100],[1,1],'--','Linewidth',1.5,'Color',[0.5 0.5 0.5])
yticks([0 1 2 3 4])
subplot(1,3,1); ylim([0,maxy(1)]); subplot(1,3,2); ylim([0,maxy(2)]); subplot(1,3,3); ylim([0,maxy(3)]);
export_fig('adherence_sens_plots/adherence_temporal_plots','-pdf','-r600')

%% Stash of plots no longer used

% %% Threshold event plots: Backwards CT & Workpercent
% variable_names = {'final_size','peak_inf','avg_isolation','duration','peak_isolation'};
% criteria_thresholds = [0.5,0.01,0.01,150,0.05];
% display_legend = true;
% make_threshold_event_plots('workpercent_and_backwardsCT',variable_names,criteria_thresholds,display_legend)


% %% temporal plots (adherence, backwards CT, work percent)
% clear variables
% % Set up plot
% if ispc==1
%     scale = 0.9/0.75;
% else
%     scale = 1;
% end
% position = [10, 10, 3.5*550*scale, 4.5*450*scale];
% set(0, 'DefaultFigurePosition', position);
% %set(0, 'DefaultFigurePosition', position,'Units','points');
% fig = figure();
% clf;
% set(fig,'Color', [1 1 1])
% fontsize=26;
% subplot(3,3,1); hold on; maxy(1) = make_temporal_plots({'adherence'},'inf_prevalence',1); ylabel('Propn infectious'); title('Adherence'); set(gca,'FontSize',fontsize); xlim([0,200]);
% yticks([0 0.05 0.10 0.15 0.20 0.25 0.30 0.35 0.40])
% subplot(3,3,2); hold on; maxy(2) = make_temporal_plots({'backwards_CT'},'inf_prevalence',1); title('Identify infector'); set(gca,'FontSize',fontsize); xlim([0,200]);
% yticks([0 0.05 0.10 0.15 0.20 0.25 0.30 0.35 0.40])
% subplot(3,3,3); hold on; maxy(3) = make_temporal_plots({'workpercent'},'inf_prevalence',1); title('Work percentage'); set(gca,'FontSize',fontsize); xlim([0,200]);
% yticks([0 0.05 0.10 0.15 0.20 0.25 0.30 0.35 0.40])
% subplot(3,3,1); ylim([0,max(maxy)]); subplot(3,3,2); ylim([0,max(maxy)]); subplot(3,3,3); ylim([0,max(maxy)]);
%
% subplot(3,3,4); hold on; maxy(1) = make_temporal_plots({'adherence'},'isol_prevalence',1); ylabel('Propn isolating'); set(gca,'FontSize',fontsize); xlim([0,200]);
% yticks([0 0.05 0.10 0.15 0.20 0.25 0.30])
% subplot(3,3,5); hold on; maxy(2) = make_temporal_plots({'backwards_CT'},'isol_prevalence',1); set(gca,'FontSize',fontsize); xlim([0,200]);
% yticks([0 0.05 0.10 0.15 0.20 0.25 0.30])
% subplot(3,3,6); hold on; maxy(3) = make_temporal_plots({'workpercent'},'isol_prevalence',1); set(gca,'FontSize',fontsize); xlim([0,200]);
% yticks([0 0.05 0.10 0.15 0.20 0.25 0.30])
% subplot(3,3,4); ylim([0,max(maxy)]); subplot(3,3,5); ylim([0,max(maxy)]); subplot(3,3,6); ylim([0,max(maxy)]);
%
% subplot(3,3,7); hold on; maxy(1) = make_temporal_plots({'adherence'},'Rt',1); ylabel('R_t'); set(gca,'FontSize',fontsize); xlim([0,100]); plot([0,100],[1,1],'--','Linewidth',1.5,'Color',[0.5 0.5 0.5])
% yticks([0 1 2 3 4])
% subplot(3,3,8); hold on; maxy(2) = make_temporal_plots({'backwards_CT'},'Rt',1); xlabel('Time(days)'); set(gca,'FontSize',fontsize); xlim([0,100]); plot([0,100],[1,1],'--','Linewidth',1.5,'Color',[0.5 0.5 0.5])
% yticks([0 1 2 3 4])
% subplot(3,3,9); hold on; maxy(3) = make_temporal_plots({'workpercent'},'Rt',1); set(gca,'FontSize',fontsize); xlim([0,100]); plot([0,100],[1,1],'--','Linewidth',1.5,'Color',[0.5 0.5 0.5])
% yticks([0 1 2 3 4])
% subplot(3,3,7); ylim([0,4]); subplot(3,3,8); ylim([0,4]); subplot(3,3,9); ylim([0,4]);
% export_fig('misc_plots/temporal_plots','-pdf','-r600')

% %% temporal plots (worker patterns)
% clear variables
% % Set up plot
% if ispc==1
%     scale = 0.9/0.75;
% else
%     scale = 1;
% end
% position = [10, 10, 3.5*550*scale, 4.5*450*scale];
% set(0, 'DefaultFigurePosition', position);
% %set(0, 'DefaultFigurePosition', position,'Units','points');
% fig = figure();
% clf;
% set(fig,'Color', [1 1 1])
% fontsize=26;
% subplot(3,2,1); hold on; maxy(1) = make_temporal_plots({'synch'},'inf_prevalence',1); ylabel('Propn infectious'); title('Synchronous'); set(gca,'FontSize',fontsize); xlim([0,200]);
% yticks([0 0.05 0.10 0.15 0.20])
% subplot(3,2,2); hold on; maxy(2) = make_temporal_plots({'asynch'},'inf_prevalence',1); title('Asynchronous'); set(gca,'FontSize',fontsize); xlim([0,200]);
% yticks([0 0.05 0.10 0.15 0.20])
% subplot(3,2,1); ylim([0,max(maxy)]); subplot(3,2,2); ylim([0,max(maxy)]);
%
% subplot(3,2,3); hold on; maxy(1) = make_temporal_plots({'synch'},'isol_prevalence',1); ylabel('Propn isolating'); set(gca,'FontSize',fontsize); xlim([0,200]);
% yticks([0 0.05 0.10 0.15 0.20])
% subplot(3,2,4); hold on; maxy(2) = make_temporal_plots({'asynch'},'isol_prevalence',1); set(gca,'FontSize',fontsize); xlim([0,200]);
% yticks([0 0.05 0.10 0.15 0.20])
% subplot(3,2,3); ylim([0,max(maxy)]); subplot(3,2,4); ylim([0,max(maxy)]);
%
% subplot(3,2,5); hold on; maxy(1) = make_temporal_plots({'synch'},'Rt',1); ylabel('R_t'); set(gca,'FontSize',fontsize); xlim([0,100]); plot([0,100],[1,1],'--','Linewidth',1.5,'Color',[0.5 0.5 0.5])
% yticks([0 1 2 3 4])
% subplot(3,2,6); hold on; maxy(2) = make_temporal_plots({'asynch'},'Rt',1); xlabel('Time(days)'); set(gca,'FontSize',fontsize); xlim([0,100]); plot([0,100],[1,1],'--','Linewidth',1.5,'Color',[0.5 0.5 0.5])
% yticks([0 1 2 3 4])
% subplot(3,2,5); ylim([0,4]); subplot(3,2,6); ylim([0,4]);
%
% num_lines = 6;
% basecolour = [0 0 0];
% light = [0.8 0.8 0.8];
% colour_vec = [linspace(light(1),basecolour(1),num_lines)', linspace(light(2),basecolour(2),num_lines)', linspace(light(3),basecolour(3),num_lines)'];
%
% legend_label = {'0 days','1 day','2 days','3 days','4 days','5 days'};
% H = gobjects(num_lines,1);
% for line_itr = 1:num_lines
%     H(line_itr) = fill(NaN,NaN,colour_vec(line_itr,:),...
%                                 'DisplayName',legend_label{line_itr});
%                                 %'FaceAlpha',0.25);
% end
% legend(H,...
%     'LineWidth',1.5,...
%     'FontSize',fontsize,...
%     'Position',[0.839220779220779 0.77 0.05 0.10]);
% export_fig('worker_pattern_plots/worker_pattern_temporal_plots','-pdf','-r600')
