%clear all;
%%
%addpath('..\..\..\codes\structures');
%%
[~, params,~] = dataexperiment(0.0001)
beam = beammodelAPH(3,1,0.01,params);

% force and speed
Cv = beam.getCvelocity(1.36,0);
Bf = beam.getBforce(1.36,0);

% voltage and its conjugate
Bv1 = beam.getBpiezo(1);
Bv2 = beam.getBpiezo(-1);
Bvoltage = Bv1+Bv2;
%
A = beam.getA;

Cp = beam.getCpiezo(1)+beam.getCpiezo(-1);
%%
sys = ss(A, [Bvoltage Bf],[Cp;Cv],0);
figure;rlocus(sys(1,1))

%%
figure; bode(feedback(sys(1,1),0)); hold all;bode(feedback(sys(1,1),1000000));

%%
Q = eye(size(A));
R = 100;
L = lqr(A',Cv',Q,R)';
%%
k=1000000;
obsv = ss(A-L*Cv-k*Bvoltage*Cp, L, Cp, 0);

controller = k*obsv;
sysf = feedback(sys(2,1),controller);
 
%figure; impulse(sys(2,1)/tf([1 0],1)); hold all;impulse(sysf/tf([1 0],1));
%figure; bode(sys(2,1)); hold all;bode(sysf);

%% test on higher order dynamics
beam = beammodelAPH(10,5,0.01,params);

% force and speed
Cvho = beam.getCvelocity(1.36,0);

% voltage and its conjugate
Bv1 = beam.getBpiezo(1);
Bv2 = beam.getBpiezo(-1);
Bvoltageho = Bv1+Bv2;
%
Aho = beam.getA;

sysho = ss(Aho, [Bvoltageho],[Cvho],0);

%%
sysf = feedback(sysho,controller);
 
%figure; impulse(sysho/tf([1 0],1)); hold all;impulse(sysf/tf([1 0],1));
%figure; bode(sys(2,1)); hold all;bode(sysf);
