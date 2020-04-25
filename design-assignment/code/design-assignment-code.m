%% ENM 056 : Machine Design Assignment
% V 1.0
% Responsible: ENM 056 Teaching group
% Contact: sharman@chalmers.se

% **************************************************************************
% This is a template file to help the studnets start up machine design 
% assignment in ENM 056. Please fill in the relevant sections one by one. 
% In order to run the code section by section, use "Run Section" option in
% MATLAB or press Ctrl + Enter. 
% Note: Please be careful when changing variable names in the template
% file. Also, the studen    ts may still need to write some sections by
% themselves.
% **************************************************************************

clear
close all
clc

%% Parameters of the reference machine
% Please fill the parameters of the machine below in SI units

mm = 1e-3; % mm to SI unit
OD_stator = 176 * mm; % Outer diameter of stator 
ID_stator = 124 * mm; % Inner diameter of stator
OD_rotor =122*mm  ; % Outer diameter of rotor
ID_rotor =  60*mm; % Inner diameter of rotor
L_stack = 100*mm  ; % Stack length
Hs0 =0.5*mm       ; % Slot opening height
Hs1 =0.5*mm       ; % Slot wedge height
Hs2 =14*mm       ; % Slot body height
w_tooth = 4.4*mm  ; % Tooth width
Rs =0.5*mm        ; % Slot bottom radius fillet
Bs0 =2*mm       ; % Slot opening
N_pole =8    ; % Number of poles
N_slot = 48   ; % Number of slots
w_layer = 4  ; % Number of winding layers
k_p =  5     ; % Coil pitch
D_wire =  0.72*mm  ; % Wire diameter
f_Cu_max =0.45  ; % Maximum Cu-fill factor
V_DC =400      ; % DC link voltage [V]
mu_0 = 4 * pi * 1e-7; % Magnetic permeability of vacuum [H/m]
Omega_base = 4000;    % Base speed [rpm]
Omega_max = 12000        ; % Maximum speed [rpm]
N_parallel =   4    ; % Number of parallel branch
%N_strand =   2        ; % Number of strands
rho_Cu = 1.72e-8; % Resistivity of Cu [ohm/m]
t_mag = 5 * mm ; % Thickness of magnet segment
w_mag =  20 * mm; % Width of magnet segment

%% Load line

% Start by assuming an air-gap flux density [T]
B_gap = 0.01: 0.001: 1 ;

% The effective air-gap cross section orthogonal to flux crossing the [sq.m]
% air-gap
A_gap = 2.5/48*pi*(ID_stator+OD_stator)/2*L_stack; 
l_gap=(ID_stator-OD_rotor)/2; % length of airgap
% Flux in the air-gap [Wb]
Phi_gap =B_gap*A_gap;

% Cross-section of tooth perpendicular to flux [sq.m]
A_tooth =w_tooth*2.5*L_stack ;

t_yoke =OD_stator-ID_stator-Hs2-Hs1-Hs0 ; % yoke thickness [m]
A_rotor=L_stack*w_mag;
% Cross-section of yoke perpendicul+ar to flux [sq.m]
A_yoke =((OD_stator-ID_stator)/2-(Hs0+Hs1+Hs2+Rs))*L_stack ;

% Output values
disp('Cross-section of different parts')
fprintf('Air-gap cross-section = % .2f [mm^2] \n',A_gap * 1e6)
fprintf('Stator tooth cross-section = % .2f [mm^2] \n',A_tooth * 1e6)
fprintf('Stator yoke cross-section = % .2f [mm^2] \n \n',A_yoke * 1e6)

% Flux densities in different parts of the machine [T]
B_tooth = Phi_gap/A_tooth; % Flux density in stator tooth
B_yoke = Phi_gap/A_yoke; % Flux density in stator yoke 
B_rotor =Phi_gap/A_rotor;

index = B_gap == 0.7;

% Output corresponding flux densities for air-gap flux density of 0.7 [T] in
% air-gap
fprintf('Flux densities in the different part of the machine when air-gap flux density is % .2f [T]\n', B_gap(index))
fprintf('Stator tooth flux density = % .2f [T] \n',B_tooth(index))
fprintf('Stator yoke flux density = % .2f [T] \n',B_yoke(index))
fprintf('Rotor yoke flux density = % .2f [T] \n \n',B_rotor(index))

% Import the B-H Curve of the M235-35A Steel from a TAB file
BH_data = importdata('BHCurve1.tab'); % import data
H_data = BH_data(:,1); % copy (Row All , Column One) as H
B_data = BH_data(:,2); % copy (Row All , Column Two) as B


% Calculated magnetic field intensity [A/m]
% Interpolation
method   = 'spline'; % 'linear' or 'spline' can be selected as interpolation method
H_rotor = interp1(B_data,H_data,B_rotor,method); % Interpolat rotor flux density, B_rotor to calculate corresponding H_rotor
H_stator_tooth =interp1(B_data,H_data,B_tooth ,method);  % Interpolat tooth flux density, B_tooth to calculate corresponding H_tooth
H_stator_yoke =interp1(B_data,H_data,B_yoke,method);  % Interpolat tooth flux density, B_yoke to calculate corresponding H_yoke
H_gap = B_gap/mu_0 ; % Calculate magnetic field intensity in the air-gap. Note: Air-gap does not contain iron.

% Output magnetic field intensities
fprintf('Magnetic field intensities in the different part of the machine when air-gap flux density is % .2f [T] \n', B_gap(index))
fprintf('Air-gap Magnetic field intensity = % .2f [A/m] \n',H_gap(index))
fprintf('Stator tooth Magnetic field intensity = % .2f [A/m] \n',H_stator_tooth(index))
fprintf('Stator yoke Magnetic field intensity = % .2f [A/m] \n',H_stator_yoke(index))
fprintf('Rotor yoke Magnetic field intensity = % .2f [A/m] \n \n',H_rotor(index))
figure(52)
plot(H_data,B_data)
hold on
plot(H_stator_tooth(index),B_tooth(index),'x') 
plot(H_stator_yoke(index),B_yoke(index),'*')
legend('B-H curve' , 'Stator tooth flux density','Stator Yoke flux density')
xlabel(' Flux density[T]')
ylabel('Magnetic field strength[A/m]')     
hold off
% Length of flux path [m]
l_rotor = 1/2*pi*(4/48*pi*OD_rotor)-2*t_mag; % Length of flux path in rotor yoke
l_stator_tooth =Hs0+Hs1+Hs2+Rs ; % Length of flux path in stator tooth
dy1=ID_stator+2*(Hs0+Hs1+Hs2+Rs);
a1=(dy1+OD_stator)/24;
a2=(OD_stator-dy1)/2;
l_stator_yoke =pi*a1+a2 ; % Length of flux pat;h in stator yoke

l_gap = ( ID_stator - OD_rotor ) / 2; % Length of flux path in air-gap

% Output values
disp('Length of flux path in different parts of the machine')
fprintf('Length of flux path in rotor = % .2f [mm] \n',l_rotor * 1e3)
fprintf('Length of flux path in stator yoke = % .2f [mm] \n',l_stator_yoke * 1e3)
fprintf('Length of flux path in stator tooth = % .2f [mm] \n',l_stator_tooth * 1e3)
fprintf('Length of flux path in air-gap = % .2f [mm] \n \n',l_gap * 1e3)

% MMF drop in different flux path [Aturn]
MMF_rotor = H_rotor * l_rotor; % MMF drop in rotor yoke
MMF_stator_yoke =H_stator_yoke * l_stator_yoke ; % MMF drop in stator yoke
MMF_stator_tooth =H_stator_tooth * l_stator_tooth ; % MMF drop in stator tooth
MMF_gap = H_gap * l_gap ; % MMF drop in air-gap

% Output values
fprintf('MMF drop in different parts of the machine when air-gap flux density is % .2f [T] \n', B_gap(index))
fprintf('MMF drop in rotor yoke = % .2f [A-turn] \n',MMF_rotor(index))
fprintf('MMF drop in stator tooth = % .2f [A-turn] \n',MMF_stator_tooth(index))
fprintf('MMF drop in stator yoke = % .2f [A-turn] \n',MMF_stator_yoke(index))
fprintf('MMF drop in air-gap = % .2f [A-turn] \n\n',MMF_gap(index))

% Total MMF drop
MMF_total =MMF_rotor+MMF_stator_yoke+ 2 *MMF_stator_tooth+ 2 *MMF_gap ;

MMF_iron =MMF_rotor+MMF_stator_yoke+MMF_stator_tooth ; % Total MMF drop in the iron parts
MMF_air =MMF_gap ; % Total MMF drop in air-gap

% Plot load line
 figure(1)               % create Figure 1
 clf                     % clear figure
 subplot(1,2,1)
 plot(MMF_total,Phi_gap*1e3, 'LineWidth', 2)
 xlabel('Total MMF drop [A-turn]')
 ylabel('Flux in air-gap[mWb]')
 legend('Load line')
 grid on

 subplot(1,2,2)
 hold on
 plot(MMF_iron,Phi_gap*1e3, 'LineWidth', 2)
 plot(MMF_air,Phi_gap*1e3, 'LineWidth', 2)
 xlabel('MMF drop [A-turn]')
 ylabel('Flux in air-gap[mWb]')
 legend('MMF drop in iron','MMF drop in air')
 grid on

%% Magnet dimension



% Magnet data 
B_mag = [0  0.5912   1.1824]; 
H_mag = [-902285    -451142     0];

MMF_mag = 2 * H_mag * t_mag; % MMF produced by magnet
Phi_mag = B_mag*w_mag*L_stack; % Flux produced by the magnet

% Magnet characteristic
figure(2)               % create Figure 1
clf                     % clear figure
plot(MMF_mag,Phi_mag * 1e3, 'LineWidth', 2)
xlim([min(MMF_mag) 0])
xlabel('MMF produced by magnet [A-turn]')
ylabel('Flux due to magnet [mWb]')
legend('Demagnetization characteristic')
grid on

%% Finding the intersection
% To avoid extrapolation of data, we can limit the x-axis of the load line
% to maximum MMF that can be produced by magnet
index = MMF_total < max(abs(MMF_mag));
MMF_total_con = MMF_total(index);
Phi_gap_con = Phi_gap(index);

% Interpolate magnet characteristic corresponding to total MMF drop to have
% same x-axis. Don't forget the negative sine to move it to second quadrant
Phi_mag_interp = interp1( abs(MMF_mag), Phi_mag , MMF_total_con , method);
[value, index] = min(abs(Phi_mag_interp - Phi_gap_con));

Phi_gap_no_load = Phi_gap_con(index);
MMF_drop_iron = MMF_total_con(index);
B_gap_no_load = Phi_gap_no_load / A_gap;
B_tooth_no_load = Phi_gap_no_load / A_tooth; % Flux density in stator tooth
B_yoke_no_load = Phi_gap_no_load / A_yoke; % Flux density in stator yoke



fprintf('Magnet thickness = % .2f [mm] \n',t_mag * 1e3)
fprintf('Magnet width = % .2f [mm] \n',w_mag * 1e3)
fprintf('No load flux in the air-gap = % .5f [Wb] \n',Phi_gap_no_load)
fprintf('No load flux density in the air-gap = % .2f [T] \n',B_gap_no_load)
fprintf('No load flux density in stator tooth = % .2f [T] \n',B_tooth_no_load)
fprintf('No load flux density in stator yoke = % .2f [T] \n\n',B_yoke_no_load)

%% Perform sensitivity analysis
% analysis 0: Do nothing
% analysis 1: Change in magnet thickness

analysis = 1;

switch analysis
    
    case 0
    % DO NOTHING    
    
    case 1       % Change in magnet thickness
        
        t_mag_1 = t_mag * 1.1; % 10% increase in magnet thickness
        t_mag_2 = t_mag * 0.9; % 10% decrease in magnet thickness
        
        MMF_mag_1 = 2 * H_mag * t_mag_1;
        MMF_mag_2 = 2 * H_mag * t_mag_2;
        
        % Magnet characteristic
        figure(4)               % create Figure 1
        clf                     % clear figure
        hold on
        plot(MMF_mag,Phi_mag * 1e3, 'LineWidth', 2)
        plot(MMF_mag_1,Phi_mag * 1e3, 'LineWidth', 2)
        plot(MMF_mag_2,Phi_mag * 1e3, 'LineWidth', 2)
        plot(-MMF_total,Phi_gap * 1e3, 'LineWidth', 2)
        hold off
        xlim([min(MMF_mag_1) 0])
        xlabel('MMF [A-turn]')
        ylabel('Flux [mWb]')
        legend('Magnet thickness 5 mm' , 'Magnet thickness 5.5 mm','Magnet thickness 4.5 mm', 'Load line')
        grid on
     
end

%% Slot area [sq.m]
% Tooth base [m]
Bs1=pi*(ID_stator+2*(Hs0+Hs1))/N_slot-w_tooth;
Bs2=pi*(ID_stator+2*(Hs0+Hs1+Hs2))/N_slot-w_tooth;
%Tooth_base = ;

% Tooth area [sq.m]
%A_tooth = ;

% Tooth and slot area [sq.m]
%A_tooth_slot = ;

% Slot area [sq.m]
A_slot = (Bs1+Bs2)/2*Hs2+(Bs2-2*Rs)*Rs+pi/2*Rs^2;

% Tooth area [sq.m]
A_tooth = w_tooth*(Hs2 + Rs + Hs1 + Hs0);
%% Stator winding design 
N_phase=3; %Antal faser
r=2;
omega_mech=Omega_base*2*pi/60;
omega_elec = omega_mech*N_pole/2; % Electrical angular frequency [rad/s]
f_elec = omega_elec/(2*pi); % Electrical frequency [Hz]

q =N_slot/(N_phase*N_pole) ; % Number of slots per pole per phase
alpha_electrical=360/(N_phase*q*2);
alpha_mechanical=alpha_electrical*2/N_pole;
tau=N_slot/N_pole;
y_1=tau-1;
k_w1=sind(q*alpha_electrical/2)/(q*sind(alpha_electrical/2))*sind(y_1/tau*90);
% Total number of turns per phase
%E_phase=400;
%N_turn_total=E_phase*N_parallel/(sqrt(2)*pi*omega_elec*Phi_gap_no_load*k_w1*q*r*N_pole/2)*N_slot/2;
N_turn_coil = floor( 400     /   (   sqrt(3)    *(N_pole/2)   *omega_mech*    Phi_gap_no_load*    q   *   k_w1    *   N_pole  /N_parallel  )    )   ; 


% Number of series branch
N_series =N_pole/N_parallel ;

% Number of coils in series per phase
N_coil_ph = q*r/2*N_series;
N_turn_ph = N_turn_coil*N_coil_ph; % Total turns per phase of the machine
% Number of turns per coil
%N_turn_coil = round (N_ph_turn / N_coil_ph);

N_strand = floor((A_slot*f_Cu_max)/(r*N_turn_coil*pi/4*D_wire^2));

% Number of conductors per slot
N_cond_slot = N_turn_coil*r*N_strand;



% Cu area
A_Cu =(D_wire/2)^2*pi*N_strand*N_turn_coil*r ;
% Fill factor
f_Cu = A_Cu / (A_slot);

% Maximum induced voltage at no-load
E_no_load =sqrt(2)*pi*f_elec*N_turn_coil*Phi_gap_no_load*k_w1*q*r*(N_pole/2)/N_parallel ;

% Maximum line to line induced voltage at no -load
E_LL_no_load =E_no_load*sqrt(3);

% Maximum line to line induced voltage at max speed
omega_mech_max=Omega_max*2*pi/60;
omega_elec_max = omega_mech_max*N_pole/2; % Electrical angular frequency [rad/s]
f_elec_max = omega_elec_max/(2*pi); % Electrical frequency [Hz]

E_LL_max_speed =sqrt(3)*N_turn_coil*sqrt(2)*pi*f_elec_max*Phi_gap_no_load*k_w1*q*r*(N_pole/2)/N_parallel ;

fprintf('Total number of turns per phase = % .0f \n',N_turn_ph )
fprintf('Total number of coils in series per phase = % .0f \n',N_coil_ph )
fprintf('Total number of turns per coil = % .0f \n',N_turn_coil )
fprintf('Total number of conductors per slot = % .0f \n',N_cond_slot )
fprintf('Slot area = % .1f [mm^2] \n',A_slot * 1e6)
fprintf('Cu area = % .1f [mm^2] \n',A_Cu * 1e6)
fprintf('Fill factor = % .2f \n\n',f_Cu)

fprintf('Maximum per phase induced voltage at no-load = % .1f [V] \n',E_no_load)
fprintf('Maximum line to line induced voltage at no-load = % .1f [V] \n',E_LL_no_load)
fprintf('Maximum line to line induced voltage at %.0f [rpm] = % .1f [V] \n',Omega_max, E_LL_max_speed)
fprintf('Ratio of maximum line to line induced and DC link voltage at max speed of %.0f [rpm] = % .1f \n\n',Omega_max, E_LL_max_speed / V_DC)

%% Resistance and inductance calculation
l_coil=2*L_stack+2*0.8*L_stack; %+pi*((Bs1+Bs2)/2+w_tooth)*y_1 ; % average length of a turn
%d_strand=sqrt(((Bs1+Bs2)/2+w_tooth)^2*y_1^2+((Hs2+Rs)/2)^2 ); % Strand diameter
R_phase =rho_Cu* (N_turn_coil*l_coil/(pi/4*D_wire^2)) *   (r*N_pole/2*q/(N_strand*N_parallel^2)) ;  % Resistance per phase [ohm]


Re_mag = abs(MMF_mag(2))/Phi_mag(2); %don't multiply by 2
Re_gap =MMF_gap(index)/Phi_gap_no_load ;
Re_stator_tooth =MMF_stator_tooth(index)/Phi_gap_no_load ;
Re_stator_yoke =MMF_stator_yoke(index)/Phi_gap_no_load ;
Re_rotor_yoke =MMF_rotor(index)/Phi_gap_no_load ;

Re_d = (MMF_total(index)/Phi_gap_no_load)+Re_mag;
Re_q =MMF_total(index)/Phi_gap_no_load ;
mu_3=abs(B_mag(2)/(H_mag(2)*mu_0));
N =N_turn_coil*k_w1*q*r ;

L_d = N^2*N_pole/2/(Re_d*N_parallel);
L_q = N^2*N_pole/(2*Re_q*N_parallel);

saliency = L_q / L_d;

fprintf('Resistance per phase = % .1f [mOhm] \n',R_phase*1e3)
fprintf('Reluctance of d-flux path = % .1f [H^-1] \n',Re_d)
fprintf('Reluctance of q-flux path = % .1f [H^-1] \n',Re_q)
fprintf('D-axis inductance = % .1f [mH] \n',L_d*1e3)
fprintf('Q-axis inductance = % .1f [mH] \n',L_q*1e3)
fprintf('Saliency ratio = % .1f \n',saliency)

%% Load calculation

%calculate d- and q-axis currents
I_amp= 100; % amplitude of current
omega_rotor_p4=3000*2*pi/60;
omega_elec_p4 = omega_rotor_p4*N_pole/2;
f_elec_p4 = omega_elec_p4/(2*pi);
Kh=156.201; %Hysteresis co-efficient
Kc=0.0204184; %Eddy current co-efficient
I_d=0; 
I_q=100;
I_s=sqrt(I_d^2+I_q^2);
Flux_link_PM=(Phi_gap_no_load*N_turn_coil*q*k_w1*r*N_pole/2)/(N_parallel); %
Flux_link_d=Flux_link_PM; %
Flux_link_q=I_q*L_q; %

T_em=3/2*N_pole/2*(Flux_link_d*I_q); 

U_d = (R_phase*I_d-omega_elec_p4*Flux_link_q);
U_q = (R_phase*I_q+omega_elec_p4*Flux_link_d);

U_s=sqrt(U_d^2+U_q^2);

Flux_link_s=sqrt(Flux_link_d^2+Flux_link_q^2);

a=Flux_link_s^2;
b=2*R_phase*(Flux_link_d*I_q-Flux_link_q*I_d);
c=R_phase^2*I_s-U_s^2;

omega_r_max=(-b+sqrt(b^2-4*a*c))/(2*a);
Omega_r_max=30*omega_r_max/(pi*N_pole/2);

P_em=T_em*omega_rotor_p4;
P_cu=3/2*I_s^2*R_phase;
Volume_yoke = ((OD_stator^2 - ID_stator^2)*pi/4 - (A_slot + A_tooth)*N_slot)*L_stack;
Volume_2 = A_tooth*L_stack*N_slot;
P_Iron_loss_yoke = (Kh * B_yoke_no_load^2 * f_elec_p4 + Kc * B_yoke_no_load^2 * f_elec_p4^2 ) * Volume_yoke;
P_Iron_loss_tooth= (Kh * B_tooth_no_load^2 * f_elec_p4 + Kc * B_tooth_no_load^2 * f_elec_p4^2 ) * Volume_2;
P_Iron_loss=P_Iron_loss_yoke+P_Iron_loss_tooth;
    
n=P_em/(P_em+P_cu+P_Iron_loss);
PF=(P_em+P_cu+P_Iron_loss)/(3/2*U_s*I_s);

 
Theta_test=linspace(10,180); 
I_d_test=I_amp*cosd(Theta_test);
I_q_test=I_amp*sind(Theta_test);    
Flux_link_d_test=Flux_link_PM+L_d*I_d_test; %
Flux_link_q_test=I_q_test*L_q; %
T_em_test=3/2*N_pole/2.*(Flux_link_d_test.*I_q_test-Flux_link_q_test.*I_d_test); %

U_s_test=sqrt((R_phase.*I_d_test-omega_elec_p4.*Flux_link_q_test).^2+(R_phase.*I_q_test+omega_elec_p4.*Flux_link_d_test).^2);
figure
grid on
hold on
plot(Theta_test,T_em_test)
xlabel('Current angle')
ylabel('Electromagnetic torque')
hold off
grid off

figure
grid on
hold on
plot(Theta_test,U_s_test)
xlabel('Current angle')
ylabel('Terminal voltage')
hold off
grid off

