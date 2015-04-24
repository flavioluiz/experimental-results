function coupledsystem = loadcoupled(numslosh,numflex,numtorsion,damping, filling)
% load coupled system
    %addpath(genpath('../../codes/'));

    % this files couples the analytical beam model with sloshing
    % in a rigid tank

    %% load all experimental description data from dataexperiment.m
    [pslosh pbeam prigid] = dataexperiment(filling);

    %% load sloshing model
    %sl =  sloshTF(numslosh, damping, pslosh);  % args:(number of modes, damping, parameters)

    %% load structural model
    plate = beammodelA(numflex,numtorsion,damping,pbeam); % load analytical beam model
    %args:(number of bending, torsion modes, damping, properties)

    % get A,Bv,BFM, C, Dv, DFM matrices:
    A = getA(plate);
    Bv = [getBpiezo(plate,1) getBpiezo(plate,-1)]; %cuidado aqui
    xr = prigid.xrb; yr = 0;
    BF = getBforce(plate,xr,yr);
    BM = getBmoment(plate,xr,yr);
    BMb = getBbmoment(plate,xr,yr);
    BFM = [BF BM BMb];
    Cv = getCvelocity(plate,xr,yr);
    Cw = getCangrate(plate,xr,yr);
    Cw2 = getCangrateb(plate,xr,yr);
    C = [Cv;Cw; Cw2; getCvelocity(plate,1.36,0.073); getCvelocity(plate,1.36,-0.073); getCvelocity(plate,1.105,0.073); getCvelocity(plate,1.105,-0.073); getCvelocity(plate,0.802,0.073); getCvelocity(plate,0.802,-0.073); getCvelocity(plate,0.451,0.073); getCvelocity(plate,0.451,-0.073)];

    %% coupling
    %rigid inertias:
    sl2.D = -[prigid.mrb 0 0; 0 prigid.Irb 0;0 0 prigid.Irb+prigid.Ifluid ];
    Q = [[1 0 0; 0 1 0; 0 0 1], zeros(3,size(C,1)-3)];
    %converting sloshing SS to 3x3 (3 inputs, 3 outputs, including "yaw"
    %using rigid body fluid+tank inertia):
    sln = ss([], [], [], zeros(size(sl2.D)));
    sln.D = sln.D + sl2.D;
    %function that couples structural dynamics and fluid dynamics:
    [ coupledsystem ] = coupleslosh( A,BFM,Bv,C,Q,sln);

end