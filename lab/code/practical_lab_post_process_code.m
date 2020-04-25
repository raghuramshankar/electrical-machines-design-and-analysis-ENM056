
% MATLAB Code - IM Parameter Post Process Code

% ENM056 Electrical Machines - Design and Analysis
% Division of Electric Power Engineering
% Chalmers University of Technology 

%   Verion      Date   
%   1.0         2018-10-08
%   1.1         2018-10-11
%   1.2         2018-10-12

% Contact:  Junfei Tang <junfei.tang@chalmers.se>

% Notes:
%
% (1) Please place this code in the same folder with the recorded data from
%     the measurement.
%
% (2) In the 'Load Data' process, please replace the channel number 
%     with the ones that you used in the measurement. 
%
% (3) In the 'Parameter Calculation' Section, please fill in the formulas 
%     to calculate the IM parameters.

%%

close all
clear
clc

%%

figure_configuration_code

%% Scaling

K_voltage        =  650 / 5;
K_current        =   35 / 5;
K_speed          = 1500 / 150;

%% Task 2 - Locked Rotor Test to Get R_series and L_series

% -- Load Data -- %

data = load('Task_2_IM_Locked_Rotor');

u1 = data(:,1) * K_voltage;
u2 = data(:,2) * K_voltage;
u3 = data(:,3) * K_voltage;
i1 = data(:,4) * K_current;
i2 = data(:,5) * K_current;
i3 = data(:,6) * K_current;
n  = data(:,8) * K_speed;

time = (0:1:(length(n)-1))/5e3;

% -- Process Data -- %

% current rms
u_lockedrotor = [u1,u2,u3];
i_lockedrotor = [i1,i2,i3];

U_rms_lockedrotor = mean(rms(u_lockedrotor));
I_rms_lockedrotor = mean(rms(i_lockedrotor));

% active and reactive power
p_lockedrotor = i1 .* u1 + i2 .* u2 + i3 .* u3;
q_lockedrotor = (i1.*(u2-u3) + i2.*(u3-u1) + i3.*(u1-u2)) / sqrt(3);

P_lockedrotor = mean(p_lockedrotor);
Q_lockedrotor = mean(q_lockedrotor);

% -- Plot -- %

f_grid = 50;            % [Hz]
T_grid = 1 / f_grid;    % [s]
temp = time > time(end) - 2.5 * T_grid;

figure('Name','Locked Rotor')

subplot(2,2,1)
hold on
plot(time(temp),u1(temp),'Color',color_2014b_blue)
plot(time(temp),u2(temp),'Color',color_2014b_orange)
plot(time(temp),u3(temp),'Color',color_2014b_green)
hold off
xlabel('Time [s]')
ylabel('Stator Voltage (V)')
legend('{\itu}_{s.1}','{\itu}_{s.2}','{\itu}_{s.3}')
grid on

subplot(2,2,2)
hold on
plot(time(temp),i1(temp),'Color',color_2014b_blue)
plot(time(temp),i2(temp),'Color',color_2014b_orange)
plot(time(temp),i3(temp),'Color',color_2014b_green)
hold off
xlabel('Time [s]')
ylabel('Stator Current [A]')
legend('{\iti}_{s.1}','{\iti}_{s.2}','{\iti}_{s.3}')
grid on

subplot(2,2,3)
plot(time(temp),n(temp),'Color',color_2014b_blue)
xlabel('Time [s]')
ylabel('Rotor Speed [rpm]')
legend('{\it\Omega}_{r}')
grid on

subplot(2,2,4)
hold on
plot(time(temp),p_lockedrotor(temp)/1e3,'Color',color_2014b_blue)
plot(time(temp),q_lockedrotor(temp)/1e3,'Color',color_2014b_orange)
plot(time(temp),P_lockedrotor*ones(1,length(time(temp)))/1e3,'Color',color_2014b_blue,'LineStyle','--')
plot(time(temp),Q_lockedrotor*ones(1,length(time(temp)))/1e3,'Color',color_2014b_orange,'LineStyle','--')
hold off
xlabel('Time [s]')
ylabel('Power [kW] [kVar]')
legend('{\itp}_{s}','{\itq}_{s}','{\itP}_{s}','{\itQ}_{s}')
grid on

result_lockedrotor = {
                    'U_rms_lockedrotor' , U_rms_lockedrotor , 'V';   
                    'I_rms_lockedrotor' , I_rms_lockedrotor , 'A'; 
                    'P_lockedrotor'     , P_lockedrotor/1e3 , 'kW';    
                    'Q_lockedrotor'     , Q_lockedrotor/1e3 , 'kVar';     
                };

disp(result_lockedrotor);

%% Task 3 - No Load Test to Get L_Mag

% -- Load Data -- %

data = load('Task_3_IM_No_Load');

u1 = data(:,1) * K_voltage;
u2 = data(:,2) * K_voltage;
u3 = data(:,3) * K_voltage;
i1 = data(:,4) * K_current;
i2 = data(:,5) * K_current;
i3 = data(:,6) * K_current;
n  = data(:,8) * K_speed;

time = (0:1:(length(n)-1))/5e3;

% For t > 0.4s the mashine is in steady-state no-load operation
temp = time > 0.4;
Speedfakt = 1499/mean(n(temp));
n = n * Speedfakt;

% -- Process Data -- %

% current rms
u_noload = [u1,u2,u3];
i_noload = [i1,i2,i3];

U_rms_noload = mean(rms(u_noload));
I_rms_noload = mean(rms(i_noload));

% active and reactive power
p_noload = i1 .* u1 + i2 .* u2 + i3 .* u3;
q_noload = (i1.*(u2-u3) + i2.*(u3-u1) + i3.*(u1-u2)) / sqrt(3);

P_noload = mean(p_noload);
Q_noload = mean(q_noload);

% -- Plot -- %

f_grid = 50;            % [Hz]
T_grid = 1 / f_grid;    % [s]
temp = time > time(end) - 2.5 * T_grid;

figure('Name','No Load')

clf

subplot(2,2,1)
hold on
plot(time(temp),u1(temp),'Color',color_2014b_blue)
plot(time(temp),u2(temp),'Color',color_2014b_orange)
plot(time(temp),u3(temp),'Color',color_2014b_green)
hold off
xlabel('Time [s]')
ylabel('Stator Voltage (V)')
legend('{\itu}_{s.1}','{\itu}_{s.2}','{\itu}_{s.3}')
grid on

subplot(2,2,2)
hold on
plot(time(temp),i1(temp),'Color',color_2014b_blue)
plot(time(temp),i2(temp),'Color',color_2014b_orange)
plot(time(temp),i3(temp),'Color',color_2014b_green)
hold off
xlabel('Time [s]')
ylabel('Stator Current [A]')
legend('{\iti}_{s.1}','{\iti}_{s.2}','{\iti}_{s.3}')
grid on

subplot(2,2,3)
plot(time(temp),n(temp),'Color',color_2014b_blue)
grid on
xlabel('Time [s]')
ylabel('Rotor Speed [rpm]')
legend('{\it\Omega}_{r}')
grid on

subplot(2,2,4)
hold on
plot(time(temp),p_noload(temp)/1e3,'Color',color_2014b_blue)
plot(time(temp),q_noload(temp)/1e3,'Color',color_2014b_orange)
plot(time(temp),P_noload*ones(1,length(time(temp)))/1e3,'Color',color_2014b_blue,'LineStyle','--')
plot(time(temp),Q_noload*ones(1,length(time(temp)))/1e3,'Color',color_2014b_orange,'LineStyle','--')
hold off
xlabel('Time [s]')
ylabel('Power [kW] [kVar]')
legend('{\itp}_{s}','{\itq}_{s}','{\itP}_{s}','{\itQ}_{s}')
grid on

result_noload = {
                    'U_rms_noload' , U_rms_noload , 'V';   
                    'I_rms_noload' , I_rms_noload , 'A'; 
                    'P_noload'     , P_noload/1e3 , 'kW';    
                    'Q_noload'     , Q_noload/1e3 , 'kVar';     
                };

disp(result_noload);

%% Parameter Calculation

% -- To Get R_s from Multimeter -- %

 R_s =3.9/3;            % [Ohm] stator resistance

% -- To Get R_r , L_sl and L_rl from Locked-Rotor Test -- %

 R_series =P_lockedrotor/(3*I_rms_lockedrotor^2) ;           % [Ohm] stator & rotor resistance in series
 X_series =Q_lockedrotor/(3*I_rms_lockedrotor^2) ;           % [Ohm] stator & rotor leakage reactance in series
 L_series =X_series/(100*pi) ;       	% [H]   stator & rotor leakage inductance in series

 R_r  = R_series - R_s;	% [Ohm] rotor resistance
 L_sl = L_series / 2;      % [H]   stator leakage inductance
 L_rl = L_series / 2;      % [H]   rotor leakage inductance

% -- To Get L_Mag from No-Load Test -- %

 X_Mag =Q_noload/(3*I_rms_noload^2)-X_series ;              % [Ohm] magnetizing reactance
 L_Mag = X_Mag/(100*pi);              % [H]   magnetizing inductance

% -- To Display Parameters -- %

 result_parameter = {
                     'R_s'     , R_s       , 'Ohm';	% [Ohm] stator resistance
                     'R_r'     , R_r       , 'Ohm';    % [Ohm] rotor resistance
                     'L_sl'	, L_sl*1e3  , 'mH';     % [H]   stator leakage inductance
                     'L_rl'	, L_rl*1e3  , 'mH';     % [H]   rotor leakage inductance
                     'L_Mag'	, L_Mag*1e3 , 'mH';     % [H]   magnetizing inductance
                 };

 disp(result_parameter);

%% Task 5 Direct Start

data = load('Task_5_IM_Direct_Start');

data = data';
time = data(:,1);
u1   = data(:,2);
u2   = data(:,3);
u3   = data(:,4);
i1   = data(:,5);
i2   = data(:,6);
i3   = data(:,7);
ua   = data(:,8);
ia   = data(:,9);
If   = data(:,10);
n    = data(:,11);
i6   = data(:,12);

% For t < 0.047 s the machine is at standing still
temp = time < 0.047;
SpeedOffSet = mean(n(temp));
n = n - SpeedOffSet;
% For t > 0.4 s the machine is in steady-state at no-load operation
temp = time > 0.4;
SpeedGain = 1500/mean(n(temp));
n = n * SpeedGain;

% -- Plot -- %

figure('Name','Direct Start')

clf

subplot(2,2,1)
plot(time,u3,'Color',color_2014b_green)
xlabel('Time [s]')
ylabel('Stator Voltage (V)')
legend('{\itu}_{s.3}')
grid on

subplot(2,2,2)
plot(time,i6,'Color',color_2014b_green)
xlabel('Time [s]')
ylabel('Stator Current [A]')
legend('{\iti}_{s.3}')
grid on

subplot(2,2,3)
plot(time,n,'Color',color_2014b_blue)
xlabel('Time [s]')
ylabel('Rotor Speed [rpm]')
legend('{\it\Omega}_{r}')
grid on

subplot(2,2,4)
plot(time,ia,'Color',color_2014b_red)
grid on
xlabel('Time [s]')
ylabel('DC Armature Current [A]')
legend('{\iti}_{a}')
grid on

