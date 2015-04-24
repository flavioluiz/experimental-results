%% Analytical beam model with piezoelectric patches as actuators
% this class creates a beam model (solution of uniform
% Euler-Bernoulli for bending, shaft dynamics for torsion)
% 
% constructor arguments:
%       -nben: number of bending modes
%       -ntor: number of torsion modes
%       -damp: modal damping
%       -params: structural description parameters:
%          beam data:
%           p.E (Young Modulus); p.mu (Poisson ratio); p.L (Length);
%           p.rho (density); p.h (thickness);
%           p.l (width);
%          piezoelectric material data:
%           p.Ep (Young modulus); p.tp (thickness); p.d31 (piezoelectric
%           constant);
%
%
%   after object construction, functions getA, getB*, getC* allows
% finding A,B,C matrices for different purposes (piezoelectric patch
% excitation, point force/moment in different positions (getB*), position, velocity
% in different positions (getC*).
%

classdef beammodelA
    properties
        NUM;
        M;
        K;
        ff;
        autovet;
        freord;
        ord;
        autoval;
        omeganTo;
        omeganB;
        Bvt;
        Bvb;
        prop;
        betan;
        NUMBEN;
        NUMTOR;
        damp;
    end
    methods
        function beam = beammodelA(nben,ntor, damp, params)            
                beam.NUMBEN = nben;
                beam.NUMTOR = ntor;
                beam.damp = damp;
                beam = getAB(beam,params);            
        end
        
        function beam = getAB(beam,p)
            path(path, './temp');
            % data    
            E = p.E; mu = p.mu; L = p.L; rho = p.rho; h = p.h;
            l = p.l;

            S = h*l;
            I = l*h^3/12;
            I2 = l^3*h/12;    
            J = 0.33*l*h^3;
            G = E/2/(1+mu);
            Ip = I+I2;
            
            beam.prop.S = S;
            beam.prop.rho = rho;
            beam.prop.L = L;
            beam.prop.Ip = Ip;
            beam.prop.G = G;
            beam.prop.J = J;

            Ep = p.Ep; tp = p.tp; d31 = p.d31;
    
            mup = mu;
            beta = Ep/E;
            eta = tp/h;

            betan0 = [0.1 3 6 8 10 13 15 17 19.5 22 24 26 29 31 33.3 35.6];
            for i = 1:beam.NUMBEN
                betan(i) = fsolve(@(x) cos(x*L)*cosh(x*L) + 1,betan0(i))            
            end
            beam.betan = betan;            
            omeganB = (betan).^2 * sqrt(E*I/(rho*S));
            
            n = length(omeganB);
            
            delta = 0.0000001;
            
            %K1a = E*d31*(tp+h)/4/(1-mu); 
            eta = tp/h;                
            K2a = eta*(eta+1)/(1/(1-mu) + beta/(1-mup)*(6*eta +12*eta^2+8*eta^3))*E/(1-mu)*d31/(1-mup)*beta*h^2 /tp /2;
            K1a = K2a;
            for i = 1:n
                Bv(i) = K1a* l / 2 * (beam.WBn(0.14+delta,i)-beam.WBn(0.14-delta,i))/2/delta;
            end
            
            %% torsion part

            omeganT = @(n) (2*n-1)*pi/2*sqrt(G*J/(Ip*L^2*rho));                      
            
            omeganTo = omeganT(1:beam.NUMTOR);
            beam.omeganB = omeganB;
            beam.omeganTo = omeganTo;
            n = length(omeganTo);
            
            delta = 0.0000001;
            

            for i = 1:n
                Bvt(i) = K1a* l^2 / 8 * (beam.WTn(0.14+delta,i)-beam.WTn(0.14-delta,i))/2/delta;
            end

            beam.Bvb = Bv;
            beam.Bvt = Bvt;
           
        end        
        function wnb = WBn(beam,x,i)
            betan = beam.betan;
            S = beam.prop.S;
            rho = beam.prop.rho;
            L = beam.prop.L;
            
            WB = @(x,i) - ((sin(betan(i)*x)-sinh(betan(i)*x))*(sinh(betan(i)*L)-sin(betan(i)*L))+...
                           (cosh(betan(i)*x) - cos(betan(i)*x))* (cos(betan(i)*L) +cosh(betan(i)*L)))./(sinh(betan(i)*L)-sin(betan(i)*L));
          
            xl = linspace(0,L,100);
            norma = @(i) sqrt(trapz(xl,rho*S*WB(xl,i).^2));
            %norma =@(i)sqrt(rho*S*( -(-4*L*betan(i)*cos(L*betan(i))*cosh(L*betan(i)) - L*betan(i)*cos(2*L*betan(i)) - 2*L*betan(i)*cosh(L*betan(i))^2 - L*betan(i) + 2*sin(L*betan(i))*cosh(L*betan(i)) + sin(2*L*betan(i))*cosh(L*betan(i))^2 + 2*cos(L*betan(i))*sinh(L*betan(i)) + cos(2*L*betan(i))*sinh(L*betan(i))*cosh(L*betan(i)) + sinh(L*betan(i))^3*cosh(L*betan(i)) - sinh(L*betan(i))*cosh(L*betan(i))^3 + 2*sinh(L*betan(i))*cosh(L*betan(i)))/(2*betan(i)*(sin(L*betan(i)) - sinh(L*betan(i)))^2)));
            
            wnb = WB(x,i)/norma(i);
            
            %wnb = -sqrt(2)*sqrt(betan(i))*(-sin(L*betan(i))*sin(betan(i)*x) + sin(L*betan(i))*sinh(betan(i)*x) + sin(betan(i)*x)*sinh(L*betan(i)) - cos(L*betan(i))*cos(betan(i)*x) + cos(L*betan(i))*cosh(betan(i)*x) - cos(betan(i)*x)*cosh(L*betan(i)) + cosh(betan(i)*(L - x)))/(sqrt(S)*sqrt(rho)*sqrt(4*L*betan(i)*cos(L*betan(i))*cosh(L*betan(i)) + L*betan(i)*cos(2*L*betan(i)) + 2*L*betan(i)*cosh(L*betan(i))^2 + L*betan(i) - 2*sin(L*betan(i))*cosh(L*betan(i)) - sin(2*L*betan(i))*cosh(L*betan(i))^2 - 2*cos(L*betan(i))*sinh(L*betan(i)) - cos(2*L*betan(i))*sinh(L*betan(i))*cosh(L*betan(i)) - sinh(L*betan(i))^3*cosh(L*betan(i)) + sinh(L*betan(i))*cosh(L*betan(i))^3 - 2*sinh(L*betan(i))*cosh(L*betan(i))));
        end
        
        function wtn = WTn(beam,x,i)
            Ip = beam.prop.Ip;
            G = beam.prop.G;
            J = beam.prop.J;
            rho = beam.prop.rho;
            L = beam.prop.L;
            
            omeganT = beam.omeganTo;
            WT = @(x,i) sin(omeganT(i)*sqrt(Ip/(G*J)*x));            
            xl = linspace(0,L,100);            
            normaT = @(i) sqrt(trapz(xl,rho*Ip*WT(xl,i).^2));            
            %normaT = @(i) Ip^(1/4)*sqrt(rho)*sqrt(-sqrt(G*J)*sin(2*sqrt(Ip)*L*omeganT(i)/sqrt(G*J)) + 2*sqrt(Ip)*L*omeganT(i))/(2*sqrt(omeganT(i)));
            wtn = WT(x,i)/normaT(i);
        end
        
        function A = getA(beam)            
            % get A matrix:
            dampingratio = beam.damp;
            Damp = diag(2*dampingratio*[beam.omeganB, beam.omeganTo]);
            Rigid = diag([beam.omeganB.^2 beam.omeganTo.^2]);
            A = [-Damp -Rigid; eye(size(Damp)) zeros(size(Damp))];
        end      
        
        function Bv = getBpiezo(beam,pos) %pos = +1,-1
            % calculate B matrix (piezoelectric effect):
            By = beam.Bvb;
            Bt = beam.Bvt;
            BB = [By pos*Bt];
            Bv = [BB BB*0]';        
        end
        
        function BF = getBforce(beam,xr,yr)
            % calculate BF matrix (external force at x=xr)
            for i = 1:length(beam.omeganB)
                By(i) = beam.WBn(xr,i);
            end
            for i = 1:length(beam.omeganTo)
                Bt(i) = yr*beam.WTn(xr,i);
            end
                
            BFi = [By Bt]';
            BF = [BFi; BFi*0];
        end
        function BM = getBmoment(beam,xr,yr)
            % calculate BM matrix (external moment at x=xr)
            for i = 1:length(beam.omeganB)
                    By(i) = 0; 
            end
            for i = 1:length(beam.omeganTo)
                    Bt(i) = beam.WTn(xr,i);            
            end
            BM = [By Bt]';
            BM = [BM; BM*0];
        end
        function BM = getBbmoment(beam,xr,yr)
            % calculate BM matrix (external moment at x=xr)
            delta = 1e-6;
            for i = 1:length(beam.omeganB)
                    By(i) = (beam.WBn(xr+delta,i)-beam.WBn(xr-delta,i))/2/delta; 
            end
            for i = 1:length(beam.omeganTo)
                    Bt(i) = 0;            
            end
            BM = [By Bt]';
            BM = [BM; BM*0];
        end
        function Cv = getCpos(beam,xr,yr)
            for i = 1:length(beam.omeganB)
                Cy(i) = beam.WBn(xr,i);
            end
            for i = 1:length(beam.omeganTo)
                Ct(i) = yr*beam.WTn(xr,i);
            end
            CC = [Cy Ct];
            Ctt = [CC*0 CC]; 
            Cv = Ctt;
        end        
        function Cv = getCvelocity(beam,xr,yr)            
            for i = 1:length(beam.omeganB)
                Cy(i) = beam.WBn(xr,i);
            end
            for i = 1:length(beam.omeganTo)
                Ct(i) = yr*beam.WTn(xr,i);
            end
            CC = [Cy Ct];
            Ctt = [CC CC*0]; 
            Cv = Ctt;
        end
        function Cvr = getCangrateb(beam,xr,yr)
            delta = 1e-6;
            for i = 1:length(beam.omeganB)
                Cy(i) = (beam.WBn(xr+delta,i)-beam.WBn(xr-delta,i))/2/delta; 
            end
            for i = 1:length(beam.omeganTo)
                Ct(i) = 0;
            end
            CCr = [Cy Ct];            
            Ctr = [CCr CCr*0]; 
            Cvr = Ctr;
        end
        
        function Cvr = getCangrate(beam,xr,yr)
            
            for i = 1:length(beam.omeganB)
                Cy(i) = 0;
            end
            for i = 1:length(beam.omeganTo)
                Ct(i) = beam.WTn(xr,i);
            end
            CCr = [Cy Ct];            
            Ctr = [CCr CCr*0]; 
            Cvr = Ctr;
        end
        
    end
end

