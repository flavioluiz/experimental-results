function [ coupledsystem ] = coupleslosh( A,BFM,Bv,Co,Q,sl)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    C = Co*A;
    Dv = Co*Bv;
    DFM = Co*BFM;       
       
    k1 = inv(eye(size(sl.D*Q*DFM))-sl.D*Q*DFM);
    k2 = inv(eye(size(C,1))-DFM*sl.D*Q);
    
    if size(sl.A)==[0,0] %% in case of "static sloshing" / masses/inertias etc.
        AA = [A+BFM*k1*sl.D*Q*C];
        BB = [Bv+BFM*k1*sl.D*Q*Dv];
        CC = k2*[C];
    else
        AA = [A+BFM*k1*sl.D*Q*C, BFM*k1*sl.C; sl.B*Q*k2*C, sl.A+sl.B*Q*k2*DFM*sl.C];
        BB = [Bv+BFM*k1*sl.D*Q*Dv; sl.B*Q*k2*Dv];
        CC = k2*[C DFM*sl.C];
    end
    DD = k2*Dv;
  
    coupledsystem = ss(AA,BB,CC,DD);
end

