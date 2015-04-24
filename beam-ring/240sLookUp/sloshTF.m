%% sloshTF: transfer function for rectangular tank force/moment under acceleration
% arguments: NUM- number of elements in the TF series (n of states)
%            dampi - modal damping
%            params - structure with tank dimensions, fluid density and
%                     gravity acceleration: p.a,p.b,p.h, p.g, p.rho
%            solvepitch - FLAG true by default. if false: ignore pitch eqs.
%
% To do list: some problems appears in the moment calculations!
%        -  I've just realized that the first (static) term of moment due
%        to rotation should be multiplied by TWO! then, static inertia will
%        be equal to rigid body one when filling = 100%. we need to verify
%        WHY this happen;
%        - the moment calculation is relative to the half height of fluid IN THE RECTANGULAR TANK.
%        it should not be a problem for partially filled (<70% for
%        example), since in this case the moment due to pressure
%        distributions along length will be much bigger than the moment due
%        to pressure distributions along lateral walls. However, for almost
%        full tank, equivalent height is very big and moments become
%        unrealistic big. VERIFY how to solve it (I tried to change the h/2
%        parameter but I am not sure if this is what we have to do... We
%        should verify the equations dedution to check what is the RIGHT
%        way to do it).
%
function sys = sloshTF(NUM, dampi, params, solvepitch)
    if nargin < 3
        dampi = 0.0;
        solvepitch = 1;
        params.g = 9.8;
        params.rho = 1000;
        
        L = 0.5; % tank length
        R = 0.105/2; % tank radius
        e = 0.7;
        % hs/2R tank filling level
        theta = 2*acos(1-2*e);
        %ls = R/2*sqrt(1-(2*e-1)^2);
        ls = 2*R*sqrt(1-(2*e-1)^2);
        Vcyl = R^2/2*(theta-sin(theta))*L;

        params.a = L;
        params.b = ls;
        params.h = Vcyl/ls/L;        
        
    elseif nargin == 3
        solvepitch = 1;
    end
        
        a = params.a;
        b = params.b;
        h = params.h;
        g = params.g;
        rho = params.rho;
        
        wn2 = @(n) (g*(2*n+1)*pi/a.*tanh((2*n+1)*pi *h/a));

        % forca e momento devido a aceleracao lateral
        Fsumn = @(n) wn2(n) / g * 8 * rho * a^3 * b / (pi^4 * (2*n+1)^4) * tf ([1 0 0],[1 2*sqrt(wn2(n))*dampi wn2(n)]);

        Flatsum = -rho*a*b*h;

        for i = 0:NUM
            Flatsum = Flatsum + Fsumn(i);
        end
        
        Msumn = @(n) wn2(n)/g * (8*rho*a^3*b)/(pi^4 * (2*n+1)^4)*(h/2 - (2*a*tanh((2*n+1)*pi*h/2/a)/((2*n+1)*pi)) + g/wn2(n))* tf([1 0 0], [1 2*sqrt(wn2(n))*dampi wn2(n)]);
        
        Mlatsum = -rho*a^3*b/12;
        
        for i = 0:NUM
            Mlatsum = Mlatsum + Msumn(i);
        end
        
         % forca e momento devido a rotacao
         Frotsumn = @(n) wn2(n) / g * 8 * rho * a^3 * b / (pi^4 * (2*n+1)^4) * (h/2 - 2*a*tanh((2*n+1)*pi/a*h/2)/(2*n+1)/pi + g/wn2(n)) *  tf ([1 0 0],[1 2*sqrt(wn2(n))*dampi wn2(n)]);

         Frotsum = -rho*a^3*b/12;

         for i = 0:NUM
             Frotsum = Frotsum + Frotsumn(i);
         end
        % 
        TEMP = h/2;
        %R = 0.105/2; e = (sqrt(1-b^2/4/R^2)+1)/2; TEMP = 2*R*e-R;
        Mrotsumn = @(n) ...
            -8 * rho * a^3 * b / (pi^4 * (2*n+1)^4) * (h/2 - a*tanh((2*n+1)*pi/a*h/2)/(2*n+1)/pi + g/wn2(n)) + ... % ATENCAO: multipliquei por 2 para ficar coerente.. no artigo original, (2n+1)^1 esta diferente! aparentemente faltou ainda multiplicar por 2!!
            -8 * rho * h^3 * b / (pi^4 * (2*n+1)^4) *(a/2 - 3*h*tanh((2*n+1)*pi/a*h/2)/(2*n+1)/pi)*0  + ...
            8 * rho * a^3 * b / (pi^4 * (2*n+1)^4) *wn2(n)/g* (TEMP - 2*a*tanh((2*n+1)*pi/a*h/2)/(2*n+1)/pi + g/wn2(n))^2 *tf ([1 0 0],[1 2*sqrt(wn2(n))*dampi wn2(n)]) ;
        
        Mrotsum = rho*a^3*b*g/12*tf(1,[1 0 0])*solvepitch*0; 
        
        for i = 0:NUM
            Mrotsum = Mrotsum + Mrotsumn(i);
        end

        TFM = [Flatsum, Frotsum; Mlatsum, Mrotsum];
        
        sys = minreal(ss(TFM));        
end
