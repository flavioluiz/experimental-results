clear all; clc;
%%
addpath('..\mimo_frf\');
load asymexc_1;
time_input{1} = inputvalues;
time_output{1} = outputvalues;
load asymexc_2;
time_input{2} = inputvalues;
time_output{2} = outputvalues;
load symexc_1;
time_input{3} = inputvalues;
time_output{3} = outputvalues;
load symexc_2;
time_input{4} = inputvalues;
time_output{4} = outputvalues;
%load random_phase_exc_1;
%time_input{5} = inputvalues;
%time_output{5} = outputvalues;
%load random_phase_exc_2;
%time_input{6} = inputvalues;
%time_output{6} = outputvalues;
ts = times(2) -times(1);

%%
[freq,frff] = mimo_frf_H1(time_input, time_output, ts)
%%
figure;
loglog(freq, abs(frff{1,1}));  hold all;
loglog(freq, abs(frff{1,2}));
loglog(freq, abs(frff{2,1}));
loglog(freq, abs(frff{2,2}))

%%
figure;
loglog(freq, abs(frff{1,1})); hold all;
load symexc_1;
freqresp = fft(outputvalues(:,1))./fft(inputvalues(:,1));
freqresp = freqresp(1:length(frff{1,1}));
loglog(freq, abs(freqresp)); 
%[MAG,PHASE] = bode(sys(4,1),freq*2*pi);
%loglog(freq,MAG(:));
axis([5e-1 60 1e-3 1e3])
legend('H1', 'theo');

return;
%%
figure;
semilogx(freq, phase(frff{1,1}));  hold all;
semilogx(freq, phase(frff{1,2}));
semilogx(freq, phase(frff{2,1}));
semilogx(freq, phase(frff{2,2}))
