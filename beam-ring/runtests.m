clear all;
global Tfinal input output times;
%% launch daq session
% for some reason, this removes the offset problem of sensors during a
% simulink session:
s = daq.createSession('ni')
s.addAnalogInputChannel('dev1',7,'voltage')
s.DurationInSeconds = 1
[data,time] = s.startForeground; plot(time,data);
delete(s)

%%
Tfinal = 240; % simulation time
Ts = 0.001; % 
diff_phase = 0;
PHASE_1 = pi/2;
PHASE_2 = pi/2 + diff_phase*pi/180;
AMP1 = 10;
AMP2 = 10;

%%
SAMPLES = (Tfinal/Ts);
times = (0:SAMPLES-1)*Ts;
Fs = 1/Ts;
windowSize = 4;
%% frequency sweep NUMOFTIMES: closed loop and open loop
% symmetric exc
diff_phase = 0;
PHASE_1 = pi/2;
PHASE_2 = pi/2 + diff_phase*pi/180;
run_save_test('symexc_1');
run_save_test('symexc_2');

diff_phase = 180;
PHASE_1 = pi/2;
PHASE_2 = pi/2 + diff_phase*pi/180;
run_save_test('asymexc_1');
run_save_test('asymexc_2');

diff_phase = (rand-0.5)*180;
PHASE_1 = pi/2;
PHASE_2 = pi/2 + diff_phase*pi/180;
run_save_test('random_phase_exc_1');
save('diff_phase1','diff_phase');

diff_phase = (rand-0.5)*180;
PHASE_1 = pi/2;
PHASE_2 = pi/2 + diff_phase*pi/180;
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