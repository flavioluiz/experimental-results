clear all;
global Tfinal input output times T AMP1 AMP2 Kd control_input T_STOP;
%% launch daq session
% for some reason, this removes the offset problem of sensors during a
% simulink session:
s = daq.createSession('ni')
s.addAnalogInputChannel('dev1',0,'voltage')
s.DurationInSeconds = 1;
[data,time] = s.startForeground;figure; plot(time,data);
delete(s)

%%
% symmetric at 0.5542 Hz
% symmetric at 1.196 Hz
% symmetric at 1.912 Hz
% symmetric and assymetric at 6.925 Hz and 9.25 Hz

freq1 = 0.5542;
freq2 = 1.196;
freq3 = 1.912;
freq4 = 6.925;
freq5 = 9.25;

Tfinal = 380; % simulation time


tf = 220;
T = 0:0.001:Tfinal;
f0 = 0.4;
ff = 25;
k = ((ff/f0)^(1/tf));
kamp = ((10/100)^(1/tf));

Ts = 0.001; % 

%%
SAMPLES = (Tfinal/Ts);
times = (0:SAMPLES-1)*Ts;
Fs = 1/Ts;
windowSize = 4;
%% frequency sweep NUMOFTIMES: closed loop and open loop
load controllersPAPERconstrained
Ts = 0.001;
% CONTROLLER 1 ("robust")
Kd = c2d(Krob,Ts);

Fs = 1/Ts;
KON = 1;
%
%% symmetric exc

SINE_FREQ = freq1*2*pi;
T_EXC = 60;
T_STOP = 110;
AMP1 = 100;
AMP2 = 100;
%
run_save_test('symexc_1stmode_Krob');

SINE_FREQ = freq2*2*pi;
T_EXC = 60;
T_STOP = 110;
AMP1 = 20;
AMP2 = 20;
%
run_save_test('symexc_2ndmode_Krob');

SINE_FREQ = freq3*2*pi;
T_EXC = 60;
T_STOP = 110;
AMP1 = 100;
AMP2 = 100;
%
run_save_test('symexc_3rdmode_Krob');

SINE_FREQ = freq4*2*pi;
T_EXC = 30;
T_STOP = 60;
AMP1 = 20;
AMP2 = 20;
run_save_test('symexc_4thmode_Krob');
%%
SINE_FREQ = freq5*2*pi;
T_EXC = 30;
T_STOP = 60;
AMP1 = 10;
AMP2 = 10;
run_save_test('symexc_5thmode_Krob_10V');
% 
% SINE_FREQ = freq4*2*pi;
% T_EXC = 30;
% T_STOP = 60;
% AMP1 = 20;
% AMP2 = -20;
% run_save_test('asymexc_4thmode_Krob');
% 
% SINE_FREQ = freq5*2*pi;
% T_EXC = 30;
% T_STOP = 60;
% AMP1 = 20;
% AMP2 = -20;
% run_save_test('asymexc_5thmode_Krob');


%% CONTROLLER 2 ("MM")
Kd = c2d(Krob2,Ts);
Fs = 1/Ts;
KON = 1;
%
% symmetric exc

SINE_FREQ = freq1*2*pi;
T_EXC = 60;
T_STOP = 150;
AMP1 = 100;
AMP2 = 100;
%
run_save_test('symexc_1stmode_Krob2');

SINE_FREQ = freq2*2*pi;
T_EXC = 60;
T_STOP = 150;
AMP1 = 20;
AMP2 = 20;
%
run_save_test('symexc_2ndmode_Krob2');

SINE_FREQ = freq3*2*pi;
T_EXC = 60;
T_STOP = 150;
AMP1 = 100;
AMP2 = 100;
%
run_save_test('symexc_3rdmode_Krob2');

SINE_FREQ = freq4*2*pi;
T_EXC = 30;
T_STOP = 60;
AMP1 = 20;
AMP2 = 20;
run_save_test('symexc_4thmode_Krob2');
%
SINE_FREQ = freq5*2*pi;
T_EXC = 30;
T_STOP = 60;
AMP1 = 10;
AMP2 = 10;
run_save_test('symexc_5thmode_Krob2_10V');

%%

%% no control
KON = 0;
SINE_FREQ = freq1*2*pi;
T_EXC = 60;
T_STOP = 250;
AMP1 = 100;
AMP2 = 100;
%%
run_save_test('symexc_1stmode');

SINE_FREQ = freq2*2*pi;
T_EXC = 60;
T_STOP = 250;
AMP1 = 20;
AMP2 = 20;
%
run_save_test('symexc_2ndmode');
%%
SINE_FREQ = freq3*2*pi;
T_EXC = 60;
T_STOP = 110;
AMP1 = 100;
AMP2 = 100;

run_save_test('symexc_3rdmode');
%%
SINE_FREQ = freq4*2*pi;
T_EXC = 30;
T_STOP = 60;
AMP1 = 20;
AMP2 = 20;
run_save_test('symexc_4thmode');
%%
SINE_FREQ = freq5*2*pi;
T_EXC = 30;
T_STOP = 60;
AMP1 = 10;
AMP2 = 10;
run_save_test('symexc_5thmode_10V');


