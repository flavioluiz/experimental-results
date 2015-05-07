clear all;
%load Ktor
load controllersPAPERconstrained
Tf = 1200;
Ts = 0.0125;
Kd = c2d(Krob,Ts);
Fs = 1/Ts;
%%
    set_param(gcs,'SimulationCommand','connect')
    set_param(gcs,'SimulationCommand','start')