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
% symmetric at 1.192 Hz
% symmetric and assymetric at 9.521 Hz and 10.1 Hz

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
% symmetric exc

SINE_FREQ = 1.192*2*pi;
T_EXC = 100;
T_STOP = 150;
AMP1 = 20;
AMP2 = 20;
%
run_save_test('symexc_1d192Hz_Krob');

SINE_FREQ = 9.521*2*pi;
T_EXC = 30;
T_STOP = 60;
AMP1 = 20;
AMP2 = 20;
run_save_test('symexc_9d521Hz_Krob');

SINE_FREQ = 10.1*2*pi;
T_EXC = 30;
T_STOP = 60;
AMP1 = 20;
AMP2 = 20;
run_save_test('symexc_10d1Hz_Krob');

SINE_FREQ = 9.521*2*pi;
T_EXC = 30;
T_STOP = 60;
AMP1 = 20;
AMP2 = -20;
run_save_test('asymexc_9d521Hz_Krob');

SINE_FREQ = 10.1*2*pi;
T_EXC = 30;
T_STOP = 60;
AMP1 = 20;
AMP2 = -20;
run_save_test('asymexc_10d1Hz_Krob');


% CONTROLLER 2 ("MM")
Kd = c2d(Krob2,Ts);
Fs = 1/Ts;
KON = 1;
%
% symmetric exc

SINE_FREQ = 1.192*2*pi;
T_EXC = 100;
T_STOP = 150;
AMP1 = 20;
AMP2 = 20;
run_save_test('symexc_1d192Hz_Krob2');

SINE_FREQ = 9.521*2*pi;
T_EXC = 30;
T_STOP = 60;
AMP1 = 20;
AMP2 = 20;
run_save_test('symexc_9d521Hz_Krob2');

SINE_FREQ = 10.1*2*pi;
T_EXC = 30;
T_STOP = 60;
AMP1 = 20;
AMP2 = 20;
run_save_test('symexc_10d1Hz_Krob2');

SINE_FREQ = 9.521*2*pi;
T_EXC = 30;
T_STOP = 60;
AMP1 = 20;
AMP2 = -20;
run_save_test('asymexc_9d521Hz_Krob2');

SINE_FREQ = 10.1*2*pi;
T_EXC = 30;
T_STOP = 60;
AMP1 = 20;
AMP2 = -20;
run_save_test('asymexc_10d1Hz_Krob2');


% no control
KON = 0;

SINE_FREQ = 1.192*2*pi;
T_EXC = 100;
T_STOP = 350;
AMP1 = 20;
AMP2 = 20;
run_save_test('symexc_1d192Hz');
%
SINE_FREQ = 9.521*2*pi;
T_EXC = 30;
T_STOP = 60;
AMP1 = 20;
AMP2 = 20;
run_save_test('symexc_9d521Hz');
%
SINE_FREQ = 10.1*2*pi;
T_EXC = 30;
T_STOP = 60;
AMP1 = 20;
AMP2 = 20;
run_save_test('symexc_10d1Hz');

SINE_FREQ = 9.521*2*pi;
T_EXC = 30;
T_STOP = 60;
AMP1 = 20;
AMP2 = -20;
run_save_test('asymexc_9d521Hz');

SINE_FREQ = 10.1*2*pi;
T_EXC = 30;
T_STOP = 60;
AMP1 = 20;
AMP2 = -20;
run_save_test('asymexc_10d1Hz');