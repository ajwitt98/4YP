%% Data Loading and String Creation
close all
clear all
clc
% Load in Data

[Time, Pressure1_bit, Pressure2_bit] = importfile("C:\Users\Alex Witt\Documents\Engineering work\Fourth Year\4YP\Data_w_Mass\" + ...
    "05_10_t_2", [1, Inf]);

% Required Constants

stamp_no = length(Time);
min_bit = 102.4;
max_bit = 921.6;
range_bit = max_bit - min_bit;
p_max = 500;

% Create arrays

T = Time - datetime(Time(1));
outfmt = 'hh:mm:ss.SSSSS';
timetostring = duration(T,'Format',outfmt);
time_raw_1 = seconds(timetostring);
time_raw_2 = seconds(timetostring);
time_raw_PD = seconds(timetostring);
Pressure1_raw = zeros(stamp_no,1);
Pressure2_raw = zeros(stamp_no,1);
PressurePD_raw = zeros(stamp_no,1);

%% Data Processing 
% Adjust Data for range of measurable values

for i = 1:stamp_no
    if Pressure1_bit(i) <= min_bit
        Pressure1_raw(i) = 0;
    elseif Pressure1_bit(i) >= max_bit
        Pressure1_raw(i) = p_max;
    else Pressure1_raw(i) = p_max*((Pressure1_bit(i) - min_bit)/range_bit);
    end

    if Pressure2_bit(i) <= min_bit
        Pressure2_raw(i) = 0;
    elseif Pressure2_bit(i) >= max_bit
        Pressure2_raw(i) = p_max;
    else Pressure2_raw(i) = p_max*((Pressure2_bit(i) - min_bit)/range_bit);
    end

    PressurePD_raw(i) = Pressure1_raw(i) - Pressure2_raw(i);

    if Pressure1_raw(i) == 0
        PressurePD_raw(i) = 0;
    elseif Pressure2_raw(i) == 0
        PressurePD_raw(i) = 0;
    end
end

%% Data Cleaning and Presentation
% Take lengths of vectors
[a] = find(Pressure1_raw);
[b] = find(Pressure2_raw);
[c] = find(PressurePD_raw);
% Create New Arrays
Pressure1 = nonzeros(Pressure1_raw);
Pressure2 = nonzeros(Pressure2_raw);
PressurePD = nonzeros(PressurePD_raw);
% Take time data readings
time1 = time_raw_1(a);
time2 = time_raw_1(b);
timePD = time_raw_1(c);
% Calculate averages for last 30 seconds
interval = 30;
% Pressure 1
max_time1 = max(time1);
interval_time1 = max_time1 - interval;
[~,interval_stamp_no1] = (min(abs(time1 - interval_time1)));
P1_30_avg = mean(Pressure1(interval_stamp_no1:length(a)));
P1_30 = Pressure1(interval_stamp_no1:length(a));
a_30 = find(Pressure1, 1, 'last' );
t1_30_raw = time1(interval_stamp_no1:a_30);
t1_30 = t1_30_raw - min(t1_30_raw);
P1_30_ravg = movmean(P1_30,100);
%Pressure 2
max_time2 = max(time2);
interval_time2 = max_time2 - interval;
[~,interval_stamp_no2] = (min(abs(time2 - interval_time2)));
P2_30_avg = mean(Pressure2(interval_stamp_no2:length(b)));
P2_30 = Pressure2(interval_stamp_no2:length(b));
b_30 = find(Pressure2, 1, 'last' );
t2_30_raw = time2(interval_stamp_no2:b_30);
t2_30 = t2_30_raw - min(t2_30_raw);
P2_30_ravg = movmean(P2_30,100);
% Pressure Difference
max_timePD = max(timePD);
interval_timePD = max_timePD - interval;
[~,interval_stamp_noPD] = (min(abs(timePD - interval_timePD)));
PPD_30_avg = mean(PressurePD(interval_stamp_noPD:length(c)));
PPD_30 = PressurePD(interval_stamp_noPD:length(c));
c_30 = find(PressurePD, 1, 'last' );
tPD_30_raw = timePD(interval_stamp_noPD:c_30);
tPD_30 = tPD_30_raw - min(tPD_30_raw);
PPD_30_ravg = movmean(PPD_30,100);
% Take Rolling Averages
P_1_ravg = movmean(Pressure1,500);
P_2_ravg = movmean(Pressure2,500);
P_PD_ravg = movmean(PressurePD,500);



%% Plotting Section
% Plot Graphs
cm_magma=magma(10);
cm_inferno=inferno(10);
cm_plasma=plasma(10);
cm_viridis=viridis(10);
clist = colormap(cm_plasma);
figure(1)
hold on

P1 = plot(time1,P_1_ravg,'color',clist(6,:),'LineWidth',1);
P2 = plot(time2,P_2_ravg,'color',clist(9,:),'LineWidth',1);
P3 = plot(timePD,P_PD_ravg,'color',clist(1,:),'LineWidth',1);
legend('Inlet Pressure','Outlet Pressure','Pressure Difference')
xlabel('Time (s)','FontSize',18,'Interpreter','latex')
ylabel('Pressure (mbar)','FontSize',18,'Interpreter','latex')
legend([P1 P2 P3],{'Inlet Pressure, $P_i$','Outlet Pressure, $P_2$','Pressure Difference, $P_i - P_2$'},'Interpreter','latex','FontSize',18)
ylim([20 150])
yticks(20:10:150)



% figure(2)
% plot(timePD,P_PD_ravg,'b')
% legend('Pressure Difference')
% xlabel('Time (Seconds)')
% ylabel('Pressure (mbar)')
% figure(3)
% hold on
% plot(t1_30,P1_30_ravg,'b')
% plot(t2_30,P2_30_ravg,'r')
% plot(tPD_30,PPD_30_ravg,'g')
% yline(P1_30_avg,'k')
% yline(P2_30_avg,'k')
% yline(PPD_30_avg,'k',sprintf('%0.5f',PPD_30_avg))
% xlabel('Time (Seconds)')
% ylabel('Pressure (mbar)')
% xlim([0 interval])
% legend('Inlet Pressure','Outlet Pressure','Pressure Difference')


%% Write results to another txt file for manipulation
%writematrix(P1_30,'P1_30.txt');

% %% Gaussian filtering test
% w = gausswin(1000,8);
% w_prime = w./sum(w);
% Pressure1_f = filter(w_prime,1,Pressure1);
%%
% empties = find(cellfun(@isempty,P_chunks));
% P_chunks(empties) = []; 










