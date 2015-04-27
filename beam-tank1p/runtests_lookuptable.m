clear all;
global Tfinal input output times T AMP1 AMP2 Kd control_input;
%% launch daq session
% for some reason, this removes the offset problem of sensors during a
% simulink session:
s = daq.createSession('ni')
s.addAnalogInputChannel('dev1',0,'voltage')
s.DurationInSeconds = 1;
[data,time] = s.startForeground;figure; plot(time,data);
delete(s)

%%


Tfinal = 240; % simulation time
tf = 220;
T = 0:0.001:Tfinal;
f0 = 0.4;
ff = 25;
k = ((ff/f0)^(1/tf));
kamp = ((10/100)^(1/tf));

Ts = 0.001; % 



diff_phase = 0;
PHASE_1 = pi/2;
PHASE_2 = pi/2 + diff_phase*pi/180;

AMP1 = 100*kamp.^T.*cos(PHASE_1+(2*pi*f0.*(k.^T-1)/log(k)) );
AMP2 = 100*kamp.^T.*cos(PHASE_2+(2*pi*f0.*(k.^T-1)/log(k)) );
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

% symmetric exc
KON = 1;

diff_phase = 0;
PHASE_1 = pi/2;
PHASE_2 = pi/2 + diff_phase*pi/180;
AMP1 = 100*kamp.^T.*cos(PHASE_1+(2*pi*f0.*(k.^T-1)/log(k)) );
AMP2 = 100*kamp.^T.*cos(PHASE_2+(2*pi*f0.*(k.^T-1)/log(k)) );
run_save_test('symexc_1_Krob');
run_save_test('symexc_2_Krob');
run_save_test('symexc_3_Krob');
run_save_test('symexc_4_Krob');
%
diff_phase = 180;
PHASE_1 = pi/2;
PHASE_2 = pi/2 + diff_phase*pi/180;
AMP1 = 100*kamp.^T.*cos(PHASE_1+(2*pi*f0.*(k.^T-1)/log(k)) );
AMP2 = 100*kamp.^T.*cos(PHASE_2+(2*pi*f0.*(k.^T-1)/log(k)) );
run_save_test('asymexc_1_Krob');
run_save_test('asymexc_2_Krob');
run_save_test('asymexc_3_Krob');
run_save_test('asymexc_4_Krob');

% CONTROLLER 2 (MM)
Kd = c2d(Krob2,Ts);
Fs = 1/Ts;
% symmetric exc
KON = 1;
diff_phase = 0;
PHASE_1 = pi/2;
PHASE_2 = pi/2 + diff_phase*pi/180;
AMP1 = 100*kamp.^T.*cos(PHASE_1+(2*pi*f0.*(k.^T-1)/log(k)) );
AMP2 = 100*kamp.^T.*cos(PHASE_2+(2*pi*f0.*(k.^T-1)/log(k)) );
run_save_test('symexc_1_Krob2');
run_save_test('symexc_2_Krob2');
run_save_test('symexc_3_Krob2');
run_save_test('symexc_4_Krob2');
%
diff_phase = 180;
PHASE_1 = pi/2;
PHASE_2 = pi/2 + diff_phase*pi/180;
AMP1 = 100*kamp.^T.*cos(PHASE_1+(2*pi*f0.*(k.^T-1)/log(k)) );
AMP2 = 100*kamp.^T.*cos(PHASE_2+(2*pi*f0.*(k.^T-1)/log(k)) );
run_save_test('asymexc_1_Krob2');
run_save_test('asymexc_2_Krob2');
run_save_test('asymexc_3_Krob2');
run_save_test('asymexc_4_Krob2');

% NO CONTROLLER
KON = 0;
diff_phase = 0;
PHASE_1 = pi/2;
PHASE_2 = pi/2 + diff_phase*pi/180;
AMP1 = 100*kamp.^T.*cos(PHASE_1+(2*pi*f0.*(k.^T-1)/log(k)) );
AMP2 = 100*kamp.^T.*cos(PHASE_2+(2*pi*f0.*(k.^T-1)/log(k)) );
run_save_test('symexc_1');
run_save_test('symexc_2');
run_save_test('symexc_3');
run_save_test('symexc_4');
%
diff_phase = 180;
PHASE_1 = pi/2;
PHASE_2 = pi/2 + diff_phase*pi/180;
AMP1 = 100*kamp.^T.*cos(PHASE_1+(2*pi*f0.*(k.^T-1)/log(k)) );
AMP2 = 100*kamp.^T.*cos(PHASE_2+(2*pi*f0.*(k.^T-1)/log(k)) );
run_save_test('asymexc_1');
run_save_test('asymexc_2');
run_save_test('asymexc_3');
run_save_test('asymexc_4');

return;


%%
diff_phase = (rand-0.5)*180;
PHASE_1 = pi/2;
PHASE_2 = pi/2 + diff_phase*pi/180;
run_save_test('random_phase_exc_1');
save('diff_phase1','diff_phase');

diff_phase = (rand-0.5)*180;
PHASE_1 = pi/2;
PHASE_2 = pi/2 + diff_phase*pi/180;
AMP1 = 100*kamp.^T.*cos(PHASE_1+(2*pi*f0.*(k.^T-1)/log(k)) );
AMP2 = 100*kamp.^T.*cos(PHASE_2+(2*pi*f0.*(k.^T-1)/log(k)) );
run_save_test('random_phase_exc_2');
save('diff_phase2','diff_phase');



%%
 load symexc_1
    figure(1);
    respfreq = fft(outputvalues(:,2))./fft(inputvalues(:,1));
%    
    NFFT = length(respfreq);   
    f = Fs/2*linspace(0,1,NFFT/2+1);
    respfreqf = filter(ones(1,windowSize)/windowSize,1,respfreq);
    subplot(2,1,1);
    loglog(f,abs(respfreqf(1:NFFT/2+1))); hold all;
    axis([0.3 100 0.00001 1]) 
    subplot(2,1,2)
    semilogx(f,-2500+phase(respfreqf(1:NFFT/2+1))*180/pi); hold all;
    axis([0.3 100 -2000 0]) 
   % save('dataopenloopfullone','inputvalues','outputvalues','times')


%% load open loop data and make graph
load dataopenloopfull
figure(2)
respfreq = fft(outputvalues(:,1))./fft(inputvalues);
for i = 1:1    
    respfreq = fft(outputvalues(:,i))./fft(inputvalues) + respfreq;
end
respfreq = respfreq/1;
    NFFT = length(respfreq);    
    f = Fs/2*linspace(0,1,NFFT/2+1);
    windowSize = 1;
    respfreqf = filter(ones(1,windowSize)/windowSize,1,respfreq);
    loglog(f,abs(respfreqf(1:NFFT/2+1))); hold all;

%% load closed loop data and make graph
load dataclosedloopfull
figure(2)
respfreq = fft(outputvalues(:,1))./fft(inputvalues);
NUMFF = 1;
for i = 1:NUMFF
    respfreq = fft(outputvalues(:,i))./fft(inputvalues) + respfreq;
end
respfreq = respfreq/NUMFF;
    NFFT = length(respfreq);    
    f = Fs/2*linspace(0,1,NFFT/2+1);
    windowSize = 1;
    respfreqf = filter(ones(1,windowSize)/windowSize,1,respfreq);
    loglog(f,abs(respfreqf(1:NFFT/2+1)));
    axis([0.3 10 0.00001 10])
    
    
%%

%%
[mag,phase,wout] = bode(Ghigh.NominalValue(2,1)+Ghigh.NominalValue(2,2));
%[mag,phase,wout] = bode(beammass(1,1))
loglog(wout/2/pi, mag(:)); 
%%
[magcl,phasecl,woutcl] = bode(feedback(beammass,-1*K))
loglog(woutcl/2/pi, magcl(:)); 
legend('open loop', 'closed loop', 'theory open loop', 'theory closed loop')