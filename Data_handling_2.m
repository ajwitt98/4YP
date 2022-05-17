% Write a script to statistically handle and plot data
close all
clear all
clc
% Set partitioning parameters
interval = 30;                          %Interval over which stability is required
max_drift = 100;                          %Set drift limit
gradient_cap = 0.0038;                  %Set gradient cap 0.0015
grad_smoothing = 1000;                  %Set gradient vector smoothing value 2050
band_smoothing = 1000;                   %Set the smoothing of the banding function


% Import data sets
rawData = combine;
partitionData = partition(interval,max_drift,gradient_cap,grad_smoothing,band_smoothing);

%%
% Useful constants 
n_runs = length(partitionData(:,1))-1;
rpm_col = rawData(:,2);
rpm_col = cell2mat(rpm_col(2:end));
n_rpms = length(unique(rpm_col));
rpm_vector = unique(rpm_col);
handledData = cell(n_runs+1,25);

% Define headers in cell array
handledData{1,1} = 'RPM';
handledData{1,2} = 'Inlet Flow Rate';
handledData{1,3} = 'PPM';
handledData{1,4} = 'Num SS sections';
handledData{1,5} = 'Median SS Flow Rate';
handledData{1,6} = 'Avg SS Flow Rate';
handledData{1,7} = 'SS Pres Diff';
handledData{1,8} = 'Avg SS Pres Diff';
handledData{1,9} = 'SS Inlet Pres';
handledData{1,10} = 'Avg SS Inlet Pres';
handledData{1,11} = 'SS Run Time';
handledData{1,12} = 'Avg SS Run Time';
handledData{1,13} = 'Avg Retentate Flow %';
handledData{1,14} = 'Avg Permeate Flow %';
handledData{1,15} = 'Retentate Flow %';
handledData{1,16} = 'Permeate Flow %';
handledData{1,17} = 'Mode SS Flow Rate';
handledData{1,18} = 'Mean SS Flow Rate';
handledData{1,19} = 'Fully Open FR';
handledData{1,20} = 'Fully Open PD';
handledData{1,21} = '70% Open FR';
handledData{1,22} = '70% Open PD';
handledData{1,23} = '50% Open FR';
handledData{1,24} = '50% Open PD';
handledData{1,25} = 'FO IP';
% Evaluate Items run by run
for x = 1:n_runs
    
    % Allocate the data from partitionData
    handledData{x+1,1} = partitionData{x+1,2};
    handledData{x+1,2} = partitionData{x+1,11};
    handledData{x+1,3} = partitionData{x+1,12};
    
    % Calculate number of SS sections
    n_SS_sections = length(partitionData{x+1,3});
    handledData{x+1,4} = n_SS_sections;

    % Initialise data handling vectors
    SS_run_time = zeros(n_SS_sections,1);
    SS_inlet_pressure = zeros(n_SS_sections,1);
    SS_pressure_difference = zeros(n_SS_sections,1);
    SS_flow_rate_median = zeros(n_SS_sections,1);
    Retentate_flow = zeros(n_SS_sections,1);
    Permeate_flow = zeros(n_SS_sections,1);
    SS_flow_rate_mode = zeros(n_SS_sections,1);
    SS_flow_rate_mean = zeros(n_SS_sections,1);
    FO_FR = zeros(n_SS_sections,1);
    FO_PD = zeros(n_SS_sections,1);
    FO_IP = zeros(n_SS_sections,1);
    FO_FR_70 = zeros(n_SS_sections,1);
    FO_PD_70 = zeros(n_SS_sections,1);
    FO_FR_50 = zeros(n_SS_sections,1);
    FO_PD_50 = zeros(n_SS_sections,1);

          
    % Calculate the values for each steady state section
    for j = 1:n_SS_sections
        SS_run_time(j) = max(partitionData{x+1,3}{j} - min(partitionData{x+1,3}{j}));
        SS_inlet_pressure(j) = mean(partitionData{x+1,6}{j});
        SS_pressure_difference(j) = mean(partitionData{x+1,10}{j});
        SS_flow_rate_median(j) = median(partitionData{x+1,4}{j});
        SS_flow_rate_mode(j) = mode(partitionData{x+1,4}{j});
        SS_flow_rate_mean(j) = mean(partitionData{x+1,4}{j});
        Retentate_flow(j) = SS_flow_rate_mean(j)/partitionData{x+1,11};
        Permeate_flow(j) = 1 - Retentate_flow(j);       
               
    end 
    max_SS_FR = max(SS_flow_rate_mean);
    SS_FR_diffp = (max_SS_FR - SS_flow_rate_mean)*100/max_SS_FR;
    SS_FR_pom = SS_flow_rate_mean*100/max_SS_FR;
    
    % Find the open valve pressure differences
    for j = 1:n_SS_sections
        if SS_FR_diffp(j) <= 1 && Retentate_flow(j) > 0.9
            FO_FR(j) = SS_flow_rate_mean(j);
            FO_PD(j) = SS_pressure_difference(j);
            FO_IP(j) = SS_inlet_pressure(j);
        end     
    end

    for j = 1:n_SS_sections
        if (60 <= SS_FR_pom(j)) && (SS_FR_pom(j) <= 90)
            FO_FR_70(j) = SS_flow_rate_mean(j);
            FO_PD_70(j) = SS_pressure_difference(j);
        end
        if (SS_FR_pom(j) <= 60)
            FO_FR_50(j) = SS_flow_rate_mean(j);
            FO_PD_50(j) = SS_pressure_difference(j);
        end
    end

    FO_PD = nonzeros(FO_PD);
    FO_FR = nonzeros(FO_FR);
    FO_IP = nonzeros(FO_IP);
    FO_PD_70 = nonzeros(FO_PD_70);
    FO_FR_70 = nonzeros(FO_FR_70);
    FO_PD_50 = nonzeros(FO_PD_50);
    FO_FR_50 = nonzeros(FO_FR_50);

    handledData{x+1,19} = FO_FR;
    handledData{x+1,20} = FO_PD;
    handledData{x+1,21} = FO_FR_70;
    handledData{x+1,22} = FO_PD_70;
    handledData{x+1,23} = FO_FR_50;
    handledData{x+1,24} = FO_PD_50;
    handledData{x+1,25} = FO_IP;

    
    % Calculate the average SS values for each bioreactor/rpm
    avg_SS_run_time = mean(SS_run_time);
    avg_inlet_pressure = mean(SS_inlet_pressure);
    avg_pressure_difference = mean(SS_pressure_difference);
    avg_flow_rate = mean(SS_flow_rate_median);

    % Calculate Flow Percentages
    handledData{x+1,13} = avg_flow_rate/partitionData{x+1,11};
    handledData{x+1,14} = 1 - handledData{x+1,13};
    
    % Write Data into cell array
    handledData{x+1,5} = SS_flow_rate_median;
    handledData{x+1,6} = avg_flow_rate;
    handledData{x+1,7} = SS_pressure_difference;
    handledData{x+1,8} = avg_pressure_difference;
    handledData{x+1,9} = SS_inlet_pressure;
    handledData{x+1,10} = avg_inlet_pressure;
    handledData{x+1,11} = SS_run_time;
    handledData{x+1,12} = avg_SS_run_time;
    handledData{x+1,15} = Retentate_flow;
    handledData{x+1,16} = Permeate_flow;
    handledData{x+1,17} = SS_flow_rate_mode;
    handledData{x+1,18} = SS_flow_rate_mean;
end

%% Plotting Section

% Create a set of vectors with aligned pressure difference and inlet flux
% Find number of total SS values and initialise matrices
p = sum(cell2mat(handledData((2:end),4)));
p_diff_string = zeros(p);
IF_value = zeros(p);
p_inlet_string = zeros(p);
rpm_value = zeros(p);
OF_value = zeros(p);
OF_ratio_value = zeros(p);
FO_pdiff = zeros(p);
FO_inlet = zeros(p);
FO_IP = zeros(p);
perc_pdiff = zeros(p);
perc_flow = zeros(p);
FO_rpm = zeros(p);
perc_rpm = zeros(p);
perc_s_pdiff = zeros(p);
perc_s_flow = zeros(p);
perc_s_rpm = zeros(p);
% Write values into the corresponding spots in the matrix
for j = 1:n_runs
    for k = 1:length(handledData{j+1,7})        
        p_diff_string(j,k) = handledData{j+1,7}(k);
        p_inlet_string(j,k) = handledData{j+1,9}(k);
        IF_value(j,k) = handledData{j+1,2};
        rpm_value(j,k) = handledData{j+1,1};
        OF_value(j,k) = handledData{j+1,18}(k);
        OF_ratio_value(j,k) = handledData{j+1,15}(k);        
    end
    for l = 1:length(handledData{j+1,20})
        FO_pdiff(j,l) = handledData{j+1,20}(l);
        FO_inlet(j,l) = handledData{j+1,2};
        FO_rpm(j,l) = handledData{j+1,1};
        FO_IP(j,l) = handledData{j+1,25}(l);
    end
    for m = 1:length(handledData{j+1,21})
        perc_pdiff(j,m) = handledData{j+1,22}(m);
        perc_flow(j,m) = handledData{j+1,2};
        perc_rpm(j,m) = handledData{j+1,1};
    end
        for n = 1:length(handledData{j+1,23})
        perc_s_pdiff(j,n) = handledData{j+1,24}(n);
        perc_s_flow(j,n) = handledData{j+1,2};
        perc_s_rpm(j,n) = handledData{j+1,1};
    end
end

% Remove Zeros
IF_value = nonzeros(IF_value);
p_diff_string = nonzeros(p_diff_string);
p_inlet_string = nonzeros(p_inlet_string);
rpm_value = nonzeros(rpm_value);
OF_value = nonzeros(OF_value);
OF_ratio_value = nonzeros(OF_ratio_value);
FO_pdiff = nonzeros(FO_pdiff);
FO_inlet = nonzeros(FO_inlet);
FO_rpm = nonzeros(FO_rpm);
FO_IP = nonzeros(FO_IP);
perc_pdiff = nonzeros(perc_pdiff);
perc_flow = nonzeros(perc_flow);
perc_rpm = nonzeros(perc_rpm);
perc_s_pdiff = nonzeros(perc_s_pdiff);
perc_s_flow = nonzeros(perc_s_flow);
perc_s_rpm = nonzeros(perc_s_rpm);
n_FO = length(FO_pdiff);
n_perc = length(perc_rpm);
n_perc_s = length(perc_s_rpm);

% Make inlet flow rates into a string of fixed decimal places
IF_x = cell(length(IF_value),1);
for i = 1:length(IF_value)
    IF_x{i} = num2str(IF_value(i),'%.1f');
end

%% Split Fully Open pressure difference values
FO_PD_split = zeros(n_FO,n_rpms);
FO_IF_split = zeros(n_FO,n_rpms);
FO_IP_split = zeros(n_FO,n_rpms);
FO_matched = [FO_rpm,FO_inlet,FO_pdiff,FO_IP];
FO_cell = cell(2,n_rpms);

perc_PD_split = zeros(n_perc,n_rpms);
perc_IF_split = zeros(n_perc,n_rpms);
perc_matched = [perc_rpm,perc_flow,perc_pdiff];
perc_cell = cell(2,n_rpms);

perc_s_PD_split = zeros(n_perc_s,n_rpms);
perc_s_IF_split = zeros(n_perc_s,n_rpms);
perc_s_matched = [perc_s_rpm,perc_s_flow,perc_s_pdiff];
perc_s_cell = cell(2,n_rpms);

for i = 1:n_FO
    for j = 1:n_rpms
        FO_cell{1,j} = rpm_vector(j);
        if rpm_vector(j) == FO_matched(i,1)
            FO_PD_split(i,j) = FO_matched(i,3);
            FO_IF_split(i,j) = FO_matched(i,2);
            FO_IP_split(i,j) = FO_matched(i,4);
        end
        FO_cell{2,j} = [nonzeros(FO_IF_split(:,j)),nonzeros(FO_PD_split(:,j)),nonzeros(FO_IP_split(:,j))];
    end
end

for i = 1:n_perc
    for j = 1:n_rpms
        perc_cell{1,j} = rpm_vector(j);
        if rpm_vector(j) == perc_matched(i,1)
            perc_PD_split(i,j) = perc_matched(i,3);
            perc_IF_split(i,j) = perc_matched(i,2);
        end
        perc_cell{2,j} = [nonzeros(perc_IF_split(:,j)),nonzeros(perc_PD_split(:,j))];
    end
end

for i = 1:n_perc_s
    for j = 1:n_rpms
        perc_s_cell{1,j} = rpm_vector(j);
        if rpm_vector(j) == perc_s_matched(i,1)
            perc_s_PD_split(i,j) = perc_s_matched(i,3);
            perc_s_IF_split(i,j) = perc_s_matched(i,2);
        end
        perc_s_cell{2,j} = [nonzeros(perc_s_IF_split(:,j)),nonzeros(perc_s_PD_split(:,j))];
    end
end

cm_magma=magma(n_rpms);
cm_inferno=inferno(n_rpms);
cm_plasma=plasma(n_rpms);
cm_viridis=viridis(n_rpms);

%% Outlet Flow Rate vs Pressure Drop i.e. fixed Q1
% Initialise vectors
clist = colormap(cm_plasma);
rsq = zeros(n_rpms,1);
lobf = zeros(n_rpms,2);
PD_split = zeros(p,n_rpms);
IP_split = zeros(p,n_rpms);
IF_split = zeros(p,n_rpms);
OF_split = zeros(p,n_rpms);
OF_ratio_split = zeros(p,n_rpms);
plotting_cell = cell(2,n_rpms);
str = cell(n_rpms,1);

% Make a matrix of all matched values
matched_values = [rpm_value,IF_value,p_inlet_string,p_diff_string,OF_value,OF_ratio_value];

% Split the data into one matrix for each variable that needs to be
% separated
for i = 1:p
    for j = 1:n_rpms
        plotting_cell{1,j} = rpm_vector(j);
        if rpm_vector(j) == matched_values(i,1)
            PD_split(i,j) = matched_values(i,4);
            IP_split(i,j) = matched_values(i,3);
            IF_split(i,j) = matched_values(i,2);
            OF_split(i,j) = matched_values(i,5);
            OF_ratio_split(i,j) = matched_values(i,6);
        end

        plotting_cell{2,j} = [nonzeros(IF_split(:,j)),nonzeros(PD_split(:,j)),nonzeros(IP_split(:,j)),nonzeros(OF_split(:,j)),nonzeros(OF_ratio_split(:,j))];

    end
end

hist_plot = cell(n_rpms,2);
figure(1)
hold on
for i = 1:n_rpms
    x = plotting_cell{2,i}(:,4);
    y = plotting_cell{2,i}(:,2);
    plot(x,y,'x','color',clist(i,:))
    c = polyfit(x,y,1);
    lobf(i,1) = c(1);
    lobf(i,2) = c(2);
    lobf(i,3) = c(2)/inlet(rpm_vector(i));

    % Evaluated equation y = m*x + b
    disp(['Equation is y = ' num2str(c(1)) '*x + ' num2str(c(2))])

    % Evaluate fit equation using polyval
    y_est = polyval(c,x);
    hist_plot{i,1} = x;
    hist_plot{i,2} = y - y_est;   
    
    % Add trend line to plot
    plot(x,y_est,'color',clist(i,:),'LineWidth',2)

    % Find the inlet flow rate
    flow_rate = inlet(plotting_cell{1,i});
    formatSpec_2 = '%.3g';
    str{i} = num2str(flow_rate,formatSpec_2);


   
    % Evaluate R^2 Value
    yresid = y - y_est;
    SSresid = sum(yresid.^2);
    SStotal = (length(y)-1) * var(y);
    rsq(i) = 1 - SSresid/SStotal;

end

xlabel('Outlet Flow Rate (ml/min)','FontSize',16,'Interpreter','latex')
ylabel('Axial Pressure Difference (mbar)','FontSize',16,'Interpreter','latex')
ylim([20 130])
xlim([0 24])
yticks(0:10:130)
xticks(0:2:24)
axis square
legend('',str{1},'',str{2},'',str{3},'',str{4},'',str{5},'',str{6},'Location','best');

lgd = legend;
lgd.FontSize = 12;
lgd.Title.String = 'Inlet Flow Rates (ml/min)';
lgd.Interpreter = 'latex';

figure(2)
subplot(3,2,1)
h1 = histfit(hist_plot{1,2},7,'kernel');
h1(1).FaceColor = clist(1,:);
xlabel('(Measured - Predicted) Pressure Difference (mbar)','Interpreter','latex')
ylabel('Number of Data Points','Interpreter','latex')
title(['\textbf{Inlet Flow Rate = ' str{1} ' ml/min}'],'Interpreter','latex')
xlim([-15 15])
legend(['$R^2 = ' num2str(rsq(1),'%.3f') '$'],'Interpreter','latex','FontSize',11)

subplot(3,2,2)
h2 = histfit(hist_plot{2,2},7,'kernel');
h2(1).FaceColor = clist(2,:);
xlabel('(Measured - Predicted) Pressure Difference (mbar)','Interpreter','latex')
ylabel('Number of Data Points','Interpreter','latex')
title(['\textbf{Inlet Flow Rate = ' str{2} ' ml/min}'],'Interpreter','latex')
xlim([-15 15])
legend(['$R^2 = ' num2str(rsq(2),'%.3f') '$'],'Interpreter','latex','FontSize',11)

subplot(3,2,3)
h3 = histfit(hist_plot{3,2},7,'kernel');
h3(1).FaceColor = clist(3,:);
xlabel('(Measured - Predicted) Pressure Difference (mbar)','Interpreter','latex')
ylabel('Number of Data Points','Interpreter','latex')
title(['\textbf{Inlet Flow Rate = ' str{3} ' ml/min}'],'Interpreter','latex')
xlim([-15 15])
legend(['$R^2 = ' num2str(rsq(3),'%.3f') '$'],'Interpreter','latex','FontSize',11)

subplot(3,2,4)
h4 = histfit(hist_plot{4,2},7,'kernel');
h4(1).FaceColor = clist(4,:);
xlabel('(Measured - Predicted) Pressure Difference (mbar)','Interpreter','latex')
ylabel('Number of Data Points','Interpreter','latex')
title(['\textbf{Inlet Flow Rate = ' str{4} ' ml/min}'],'Interpreter','latex')
xlim([-15 15])
legend(['$R^2 = ' num2str(rsq(4),'%.3f') '$'],'Interpreter','latex','FontSize',11)

subplot(3,2,5)
h5 = histfit(hist_plot{5,2},7,'kernel');
h5(1).FaceColor = clist(5,:);
xlabel('(Measured - Predicted) Pressure Difference (mbar)','Interpreter','latex')
ylabel('Number of Data Points','Interpreter','latex')
title(['\textbf{Inlet Flow Rate = ' str{5} ' ml/min}'],'Interpreter','latex')
xlim([-15 15])
legend(['$R^2 = ' num2str(rsq(5),'%.3f') '$'],'Interpreter','latex','FontSize',11)

subplot(3,2,6)
h6 = histfit(hist_plot{6,2},7,'kernel');
h6(1).FaceColor = clist(6,:);
xlabel('(Measured - Predicted) Pressure Difference (mbar)','Interpreter','latex')
ylabel('Number of Data Points','Interpreter','latex')
title(['\textbf{Inlet Flow Rate = ' str{6} ' ml/min}'],'Interpreter','latex')
xlim([-15 15])
legend(['$R^2 = ' num2str(rsq(6),'%.3f') '$'],'Interpreter','latex','FontSize',11)


%% Fixed Q2
s_dev_FO = zeros(n_rpms,1);
s_dev_perc = zeros(n_rpms,1);
s_dev_perc_s = zeros(n_rpms,1);
mean_perc = zeros(n_rpms,1);
mean_perc_s = zeros(n_rpms,1);
mean_x_perc = zeros(n_rpms,1);
mean_x_perc_s = zeros(n_rpms,1);
FO_mean = zeros(n_rpms,1);
FO_mean_x = zeros(n_rpms,1);
for i = 1:n_rpms
    s_dev_FO(i) = std(FO_cell{2,i}(:,2));
    s_dev_perc(i) = std(perc_cell{2,i}(:,2));
    s_dev_perc_s(i) = std(perc_s_cell{2,i}(:,2));
    FO_mean(i) = mean(FO_cell{2,i}(:,2));
    mean_perc(i) = mean(perc_cell{2,i}(:,2));
    mean_perc_s(i) = mean(perc_s_cell{2,i}(:,2));
    FO_mean_x(i) = FO_cell{2,i}(1,1);
    mean_x_perc(i) = perc_cell{2,i}(1,1);
    mean_x_perc_s(i) = perc_s_cell{2,i}(1,1);
end



cm_plasma_2=plasma(4);
clist_2 = colormap(cm_plasma_2);


figure(3)
hold on
x_FO = FO_inlet;
y_FO = FO_pdiff;
x_perc = perc_flow;
y_perc = perc_pdiff;
x_perc_s = perc_s_flow;
y_perc_s = perc_s_pdiff;
c_FO = polyfit(x_FO,y_FO,1);
c_perc = polyfit(x_perc,y_perc,1);
c_perc_s = polyfit(x_perc_s,y_perc_s,1);
lobf_FO(i,1) = c_FO(1);
lobf_FO(i,2) = c_FO(2);
lobf_perc(i,1) = c_perc(1);
lobf_perc(i,2) = c_perc(2);
lobf_perc_s(i,1) = c_perc_s(1);
lobf_perc_s(i,2) = c_perc_s(2);
errorbar(FO_mean_x,FO_mean,s_dev_FO,'o','Color',clist(1,:))
errorbar(mean_x_perc,mean_perc,s_dev_perc,'o','Color',clist(3,:))
errorbar(mean_x_perc_s,mean_perc_s,s_dev_perc_s,'o','Color',clist(5,:))
% Evaluated equation y = m*x + b
disp(['Equation is y = ' num2str(c_FO(1)) '*x + ' num2str(c_FO(2))])
disp(['Equation is y = ' num2str(c_perc(1)) '*x + ' num2str(c_perc(2))])
disp(['Equation is y = ' num2str(c_perc_s(1)) '*x + ' num2str(c_perc_s(2))])

% Evaluate fit equation using polyval
y_est_FO = polyval(c_FO,x_FO);
y_est_perc = polyval(c_perc,x_perc);
y_est_perc_s = polyval(c_perc_s,x_perc_s);

% Add trend line to plot
plot(x_FO,y_est_FO,'Color',clist(1,:),'LineWidth',2)
plot(x_perc,y_est_perc,'Color',clist(3,:),'LineWidth',2)
plot(x_perc_s,y_est_perc_s,'Color',clist(5,:),'LineWidth',2)

% Evaluate R^2 Value
yresid_FO = y_FO - y_est_FO;
SSresid_FO = sum(yresid_FO.^2);
SStotal_FO = (length(y_FO)-1) * var(y_FO);
rsq_FO = 1 - SSresid_FO/SStotal_FO;
yresid_perc = y_perc - y_est_perc;
SSresid_perc = sum(yresid_perc.^2);
SStotal_perc = (length(y_perc)-1) * var(y_perc);
rsq_perc = 1 - SSresid_perc/SStotal_perc;
yresid_perc_s = y_perc_s - y_est_perc_s;
SSresid_perc_s = sum(yresid_perc_s.^2);
SStotal_perc_s = (length(y_perc_s)-1) * var(y_perc_s);
rsq_perc_s = 1 - SSresid_perc_s/SStotal_perc_s;

str_FO = cell(2,1);
str_perc = cell(2,1);
formatSpec_2 = '%#.3g';
for i = 1:length(c_FO)
    str_FO{i,1} = num2str(c_FO(i),formatSpec_2);
    str_perc{i,1} = num2str(c_perc(i),formatSpec_2);
end

lgd = legend('','','',['$99 \% <$ of Max Flow Rate'],['$60 \% - 90 \%$ of Max Flow Rate'],['$60 \% >$ of Max Flow Rate'],'Location','best');
lgd.FontSize = 12;
lgd.Interpreter = 'latex';

t1 = text(1.02*max(FO_mean_x),max(FO_mean),['$R^2 = ' num2str(rsq_FO,'%#.3f') '$'],'Interpreter','latex','FontSize',14);
t2 = text(1.02*max(FO_mean_x),max(mean_perc),['$R^2 = ' num2str(rsq_perc,'%#.3f') '$'],'Interpreter','latex','FontSize',14);
t3 = text(1.02*max(FO_mean_x),max(mean_perc_s),['$R^2 = ' num2str(rsq_perc_s,'%#.3f') '$'],'Interpreter','latex','FontSize',14);


xlabel('Inlet Flow Rate (ml/min)','FontSize',16,'Interpreter','latex')
ylabel('Axial Pressure Difference, $P_i - P_2$ (mbar)','FontSize',16,'Interpreter','latex')
ylim([20 130])
xlim([8 26])
axis square

%% Normalised Flow Rate 
str_norm = cell(n_rpms,1);
rsq_norm = zeros(n_rpms,1);
lobf_norm = zeros(n_rpms,2);
figure(4)
hold on
for i = 1:n_rpms
    x_norm = plotting_cell{2,i}(:,5);
    y_norm = plotting_cell{2,i}(:,2);
    plot(x_norm,y_norm,'x','color',clist(i,:))
    c_norm = polyfit(x_norm,y_norm,1);
    lobf_norm(i,1) = c_norm(1);
    lobf_norm(i,2) = c_norm(2);

    % Evaluated equation y = m*x + b
    disp(['Equation is y = ' num2str(c_norm(1)) '*x + ' num2str(c_norm(2))])

    % Evaluate fit equation using polyval
    y_est_norm = polyval(c_norm,x_norm);
       
    % Add trend line to plot
    plot(x_norm,y_est_norm,'color',clist(i,:),'LineWidth',2)

    % Find the inlet flow rate
    flow_rate = inlet(plotting_cell{1,i});
    formatSpec_2 = '%.3g';
    str_norm{i} = num2str(flow_rate,formatSpec_2);

    % Evaluate R^2 Value
    yresid_norm = y_norm - y_est_norm;
    SSresid_norm = sum(yresid_norm.^2);
    SStotal_norm = (length(y_norm)-1) * var(y_norm);
    rsq_norm(i) = 1 - SSresid_norm/SStotal_norm;

    % Label Data

    text(1.08,max(y_est_norm),['$\frac{Q_2}{Q_i}= ' num2str(max(x_norm),'%#.3f') '$'],'Interpreter','latex','FontSize',14);

end
text(1.08,130,['Max $(\frac{Q_2}{Q_i})$:'],'Interpreter','latex','FontSize',14);
xlabel('Outlet/Inlet Flow Rate Ratio','FontSize',16,'Interpreter','latex')
ylabel('Axial Pressure Difference, $P_i - P_2$ (mbar)','FontSize',16,'Interpreter','latex')
ylim([0 140])
xlim([0 1.3])
legend('',str_norm{1},'',str_norm{2},'',str_norm{3},'',str_norm{4},'',str_norm{5},'',str_norm{6},'Location','best');

lgd = legend;
lgd.FontSize = 12;
lgd.Title.String = 'Inlet Flow Rates (ml/min)';
lgd.Interpreter = 'latex';
% ylim([0 130])
yticks(0:10:140)
xticks(0:0.1:1.3)
axis square

%% P in
s_dev_PI = zeros(n_rpms,1);
PI_mean = zeros(n_rpms,1);
PI_mean_x = zeros(n_rpms,1);
for i = 1:n_rpms
    s_dev_PI(i) = std(FO_cell{2,i}(:,3));
    PI_mean(i) = mean(FO_cell{2,i}(:,3));
    PI_mean_x(i) = FO_cell{2,i}(1,1);
end

figure(5)
hold on
x_PI = FO_inlet;
y_PI = FO_IP;
c_PI = polyfit(x_PI,y_PI,1);
lobf_PI(i,1) = c_PI(1);
lobf_PI(i,2) = c_PI(2);
errorbar(PI_mean_x,PI_mean,s_dev_PI,'o','Color',clist(1,:))
% Evaluated equation y = m*x + b
disp(['Equation is y = ' num2str(c_PI(1)) '*x + ' num2str(c_PI(2))])

% Evaluate fit equation using polyval
y_est_PI = polyval(c_PI,x_PI);

% Add trend line to plot
plot(x_PI,y_est_PI,'Color',clist(1,:),'LineWidth',2)
% Evaluate R^2 Value
yresid_PI = y_PI - y_est_PI;
SSresid_PI = sum(yresid_PI.^2);
SStotal_PI = (length(y_PI)-1) * var(y_PI);
rsq_PI = 1 - SSresid_PI/SStotal_PI;

str_PI = cell(2,1);
formatSpec_2 = '%#.3g';
% for i = 1:length(c_FO)
%     str_FO{i,1} = num2str(c_FO(i),formatSpec_2);
%     str_perc{i,1} = num2str(c_perc(i),formatSpec_2);
% end
% 
% lgd = legend('','','','99 % < of Max Flow Rate','60 % - 90 % of Max Flow Rate','60 % > of Max Flow Rate','Location','best');
% lgd.FontSize = 12;

xlabel('Inlet Flow Rate (ml/min)','FontSize',16,'Interpreter','latex')
ylabel(['Inlet Pressure, $P_i$ (mbar)'],'Interpreter','latex','FontSize',16,'Interpreter','latex')
% ylim([20 130])
% xlim([8 26])
% axis square

%% Outlet Flow Rate vs Pressure Drop i.e. fixed Q1
% Initialise vectors
rsq_pi = zeros(n_rpms,1);
lobf_pi = zeros(n_rpms,2);
hist_plot_pi = cell(n_rpms,2);
str_pi = cell(n_rpms,1);
figure(6)
hold on
for i = 1:n_rpms
    x_pi = plotting_cell{2,i}(:,4);
    y_pi = plotting_cell{2,i}(:,3);
    plot(x_pi,y_pi,'x','color',clist(i,:))
    c_pi = polyfit(x_pi,y_pi,1);
    lobf_pi(i,1) = c_pi(1);
    lobf_pi(i,2) = c_pi(2);
    str_pi{i} = max(y_pi);
    
    % Evaluated equation y = m*x + b
    disp(['Equation is y = ' num2str(c_pi(1)) '*x + ' num2str(c_pi(2))])

    % Evaluate fit equation using polyval
    y_est_pi = polyval(c_pi,x_pi);
    hist_plot_pi{i,1} = x_pi;
    hist_plot_pi{i,2} = y_pi - y_est_pi;   
    
    % Add trend line to plot
    plot(x_pi,y_est_pi,'color',clist(i,:),'LineWidth',2)
  
    % Evaluate R^2 Value
    yresid_pi = y_pi - y_est_pi;
    SSresid_pi = sum(yresid_pi.^2);
    SStotal_pi = (length(y_pi)-1) * var(y_pi);
    rsq_pi(i) = 1 - SSresid_pi/SStotal_pi;

end

xlabel('Outlet Flow Rate (ml/min)','FontSize',16,'Interpreter','latex')
ylabel('Lumen Inlet Pressure, $P_i$ (mbar)','FontSize',16,'Interpreter','latex')
%ylim([20 130])
xlim([0 24])
%yticks(0:10:130)
xticks(0:2:24)
axis square
legend('',[str{1} ', $P_{i,max} = ' num2str(str_pi{1},'%.3g') '$'],'',[str{2} ', $P_{i,max} = ' num2str(str_pi{2},'%.3g') '$'],'', ...
    [str{3} ', $P_{i,max} = ' num2str(str_pi{3},'%.3g') '$'],'',[str{4} ', $P_{i,max} = ' num2str(str_pi{4},'%.3g') '$'],'', ...
    [str{5} ', $P_{i,max} = ' num2str(str_pi{5},'%.3g') '$'],'', ...
    [str{6} ', $P_{i,max} = ' num2str(str_pi{6},'%.3g') '$'],'Location','best');

lgd = legend;
lgd.FontSize = 12;
lgd.Title.String = 'Inlet Flow Rates (ml/min), and corresponding $P_{i,max}$ (mbar)';
lgd.Interpreter = 'latex';

%%
figure(7)
subplot(3,2,1)
h1 = histfit(hist_plot_pi{1,2},7,'kernel');
h1(1).FaceColor = clist(1,:);
xlabel('(Measured - Predicted) Inlet Pressure, $P_i$ (mbar)','Interpreter','latex')
ylabel('Number of Data Points','Interpreter','latex')
title(['\textbf{Inlet Flow Rate = ' str{1} ' ml/min}'],'Interpreter','latex')
%xlim([-15 15])
legend(['$R^2 = ' num2str(rsq_pi(1),'%.3f') '$'],'Interpreter','latex','FontSize',11)

subplot(3,2,2)
h2 = histfit(hist_plot_pi{2,2},7,'kernel');
h2(1).FaceColor = clist(2,:);
xlabel('(Measured - Predicted) Inlet Pressure, $P_i$ (mbar)','Interpreter','latex')
ylabel('Number of Data Points','Interpreter','latex')
title(['\textbf{Inlet Flow Rate = ' str{2} ' ml/min}'],'Interpreter','latex')
%xlim([-15 15])
legend(['$R^2 = ' num2str(rsq_pi(2),'%.3f') '$'],'Interpreter','latex','FontSize',11)

subplot(3,2,3)
h3 = histfit(hist_plot_pi{3,2},7,'kernel');
h3(1).FaceColor = clist(3,:);
xlabel('(Measured - Predicted) Inlet Pressure, $P_i$ (mbar)','Interpreter','latex')
ylabel('Number of Data Points','Interpreter','latex')
title(['\textbf{Inlet Flow Rate = ' str{3} ' ml/min}'],'Interpreter','latex')
%xlim([-15 15])
legend(['$R^2 = ' num2str(rsq_pi(3),'%.3f') '$'],'Interpreter','latex','FontSize',11)

subplot(3,2,4)
h4 = histfit(hist_plot_pi{4,2},7,'kernel');
h4(1).FaceColor = clist(4,:);
xlabel('(Measured - Predicted) Inlet Pressure, $P_i$ (mbar)','Interpreter','latex')
ylabel('Number of Data Points','Interpreter','latex')
title(['\textbf{Inlet Flow Rate = ' str{4} ' ml/min}'],'Interpreter','latex')
%xlim([-15 15])
legend(['$R^2 = ' num2str(rsq_pi(4),'%.3f') '$'],'Interpreter','latex','FontSize',11)

subplot(3,2,5)
h5 = histfit(hist_plot_pi{5,2},7,'kernel');
h5(1).FaceColor = clist(5,:);
xlabel('(Measured - Predicted) Inlet Pressure, $P_i$ (mbar)','Interpreter','latex')
ylabel('Number of Data Points','Interpreter','latex')
title(['\textbf{Inlet Flow Rate = ' str{5} ' ml/min}'],'Interpreter','latex')
%xlim([-15 15])
legend(['$R^2 = ' num2str(rsq_pi(5),'%.3f') '$'],'Interpreter','latex','FontSize',11)

subplot(3,2,6)
h6 = histfit(hist_plot_pi{6,2},7,'kernel');
h6(1).FaceColor = clist(6,:);
xlabel('(Measured - Predicted) Inlet Pressure, $P_i$ (mbar)','Interpreter','latex')
ylabel('Number of Data Points','Interpreter','latex')
title(['\textbf{Inlet Flow Rate = ' str{6} ' ml/min}'],'Interpreter','latex')
%xlim([-15 15])
legend(['$R^2 = ' num2str(rsq_pi(6),'%.3f') '$'],'Interpreter','latex','FontSize',11)

%%
figure(8)
hold on
a = cell2mat(str);
a = str2num(a);
plot(a,rsq_pi,'x','color',clist(1,:),'LineWidth',1.2)
c_rsq = polyfit(a,rsq_pi,1);
disp(['Equation is y = ' num2str(c_rsq(1)) '*x + ' num2str(c_rsq(2))])
y_est_rsq = polyval(c_rsq,a);
plot(a,y_est_rsq,'color',clist(4,:),'LineWidth',2)
% Evaluate R^2 Value
yresid = rsq_pi - y_est_rsq;
SSresid = sum(yresid.^2);
SStotal = (length(rsq_pi)-1) * var(rsq_pi);
meta_rsq = 1 - SSresid/SStotal;
xlim([8 24])
xticks(8:2:24)
ylim([0.2 1])
legend('',['$R^2 = ' num2str(meta_rsq,'%.3f') '$'],'Interpreter','latex','FontSize',16)
xlabel('Inlet Flow Rate (ml/min)','FontSize',16,'Interpreter','latex')
ylabel('$R^2$ Values for each $P_i$ Line of Best Fit','FontSize',16,'Interpreter','latex')
%%
fr_low = cell(n_rpms,2);
fr_mid = cell(n_rpms,2);
fr_mid_u = cell(n_rpms,2);
fr_high = cell(n_rpms,2);

for i = 1:n_rpms
    flowdata = plotting_cell{2,i}(:,4);
    pressuredata = plotting_cell{2,i}(:,2);
    inletdata = plotting_cell{2,i}(:,1);
    find_fr_low = find(flowdata>0 & flowdata<=8);
    find_fr_mid = find(flowdata>8 & flowdata<=12);
    find_fr_mid_u = find(flowdata>12 & flowdata<=16);
    find_fr_high = find(flowdata>16 & flowdata<25);
    fr_low{i,1} = flowdata(find_fr_low);
    fr_low{i,2} = pressuredata(find_fr_low);    
    fr_mid{i,1} = flowdata(find_fr_mid);
    fr_mid{i,2} = pressuredata(find_fr_mid);
    fr_mid_u{i,1} = flowdata(find_fr_mid_u);
    fr_mid_u{i,2} = pressuredata(find_fr_mid_u);    
    fr_high{i,1} = flowdata(find_fr_high);
    fr_high{i,2} = pressuredata(find_fr_high);
end

y_low = zeros(n_rpms,1);
y_mid = zeros(n_rpms,1);
y_mid_u = zeros(n_rpms,1);
y_high = zeros(n_rpms,1);
std_low = zeros(n_rpms,1);
std_mid = zeros(n_rpms,1);
std_mid_u = zeros(n_rpms,1);
std_high = zeros(n_rpms,1);
of_low = mean(cell2mat(fr_low(:,1)));
of_mid = mean(cell2mat(fr_mid(:,1)));
of_mid_u = mean(cell2mat(fr_mid_u(:,1)));
of_high = mean(cell2mat(fr_high(:,1)));

for i = 1:n_rpms
    y_low(i) = mean(fr_low{i,2});
    std_low(i) = std(fr_low{i,2});
    y_mid(i) = mean(fr_mid{i,2});
    std_mid(i) = std(fr_mid{i,2});
    y_mid_u(i) = mean(fr_mid_u{i,2});
    std_mid_u(i) = std(fr_mid_u{i,2});
    y_high(i) = mean(fr_high{i,2});
    std_high(i) = std(fr_high{i,2});
end

x_low = a(~isnan(y_low));
std_low = std_low(~isnan(y_low));
y_low = y_low(~isnan(y_low));

x_mid = a(~isnan(y_mid));
std_mid = std_mid(~isnan(y_mid));
y_mid = y_mid(~isnan(y_mid));

x_mid_u = a(~isnan(y_mid_u));
std_mid_u = std_mid_u(~isnan(y_mid_u));
y_mid_u = y_mid_u(~isnan(y_mid_u));

x_high = a(~isnan(y_high));
std_high = std_high(~isnan(y_high));
y_high = y_high(~isnan(y_high));

figure(9)
hold on
errorbar(x_low,y_low,std_low,'x','color',clist_2(4,:))
errorbar(x_mid,y_mid,std_mid,'x','color',clist_2(3,:))
errorbar(x_mid_u,y_mid_u,std_mid_u,'x','color',clist_2(2,:))
errorbar(x_high,y_high,std_high,'x','color',clist_2(1,:))
c_low = polyfit(x_low,y_low,1);
c_mid = polyfit(x_mid,y_mid,1);
c_mid_u = polyfit(x_mid_u,y_mid_u,1);
c_high = polyfit(x_high,y_high,1);
y_est_low = polyval(c_low,x_low);
y_est_mid = polyval(c_mid,x_mid);
y_est_mid_u = polyval(c_mid_u,x_mid_u);
y_est_high = polyval(c_high,x_high);
plot(x_high,y_est_high,'color',clist_2(1,:),'LineWidth',2)
plot(x_mid_u,y_est_mid_u,'color',clist_2(2,:),'LineWidth',2)
plot(x_mid,y_est_mid,'color',clist_2(3,:),'LineWidth',2)
plot(x_low,y_est_low,'color',clist_2(4,:),'LineWidth',2)

yresid = y_low - y_est_low;
SSresid = sum(yresid.^2);
SStotal = (length(y_low)-1) * var(y_low);
rsq_low = 1 - SSresid/SStotal;

yresid = y_mid - y_est_mid;
SSresid = sum(yresid.^2);
SStotal = (length(y_mid)-1) * var(y_mid);
rsq_mid = 1 - SSresid/SStotal;

yresid = y_mid_u - y_est_mid_u;
SSresid = sum(yresid.^2);
SStotal = (length(y_mid_u)-1) * var(y_mid_u);
rsq_mid_u = 1 - SSresid/SStotal;

yresid = y_high - y_est_high;
SSresid = sum(yresid.^2);
SStotal = (length(y_high)-1) * var(y_high);
rsq_high = 1 - SSresid/SStotal;


lgd = legend('','','','',['$OFR \le 8 $'], ...
    ['$8 < OFR \le 12 $'], ...
    ['$12<OFR\le16 $'], ...
    ['$16<OFR$'],'Location','northwest');
lgd.FontSize = 12;
lgd.Interpreter = 'latex';
lgd.Title.String = ['Outlet Flow Rate (OFR) Range (ml/min)'];


t1 = text(1.02*max(a),max(y_est_high),['$R^2 = ' num2str(rsq_high,'%#.3f') '$' 10 '$OFR_{mean} = ' num2str(of_high,'%#.3g') '$'],'Interpreter','latex','FontSize',11);
t2 = text(1.02*max(a),max(y_est_mid_u),['$R^2 = ' num2str(rsq_mid_u,'%#.3f') '$' 10 '$OFR_{mean} = ' num2str(of_mid_u,'%#.3g') '$'],'Interpreter','latex','FontSize',11);
t3 = text(1.02*max(a),max(y_mid),['$R^2 = ' num2str(rsq_mid,'%#.3f') '$' 10 '$OFR_{mean} = ' num2str(of_mid,'%#.3g') '$'],'Interpreter','latex','FontSize',11);
t4 = text(1.02*max(a),max(y_est_low),['$R^2 = ' num2str(rsq_low,'%#.3f') '$' 10 '$OFR_{mean} = ' num2str(of_low,'%#.3g') '$'],'Interpreter','latex','FontSize',11);


xlabel('Inlet Flow Rate (ml/min)','FontSize',16,'Interpreter','latex')
ylabel('Axial Pressure Difference, $P_i - P_2$ (mbar)','FontSize',16,'Interpreter','latex')
ylim([20 130])
xlim([8 26])
axis square


