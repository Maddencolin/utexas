
%%
data = table2array(Fatiguedata);
%%
diameter = .192; %in
radius = diameter/2;
area = pi*radius^2;
%%
load = data(:,1);
N_sandblasted = data(:,2); %number of cycles
N_polished = data(:,3); 
%%
L=4; %in
stress_sandblasted = 16*load*L/(pi*diameter^3); %psi
stress_polished = 16*load*L/(pi*diameter^3); %psi
%%
%%
f_sand = fit(N_sandblasted,stress_sandblasted,'power1')
figure();
plot(f_sand,N_sandblasted, stress_sandblasted);
title('Number of cycles vs Stress - Sandblasted');
xlabel('Number of cycles');
ylabel('Stress (psi)');

%%
f_polished = fit(N_polished, stress_polished,'power1')
figure();
plot(f_polished,N_polished, stress_polished);
title('Number of cycles vs Stress - Polished');
xlabel('Number of cycles');
ylabel('Stress (psi)');
%% problem 3
N_cycles = [0, 1540, 2740, 3440, 3900, 4190, 4300, 4330];
cracklength = [1.25, 1.375, 1.5, 1.625, 1.75, 1.875, 2, 2.125]; %in
f = fit(N_cycles', cracklength', 'exp1')
figure();
plot(f,N_cycles, cracklength);
title('Crack Length vs Number of Cycles');
xlabel('N - Number of Cycles');
ylabel('Crack Length (in)');
%%
fittedCrack = 1.143*exp(0.000124*N_cycles)
%%
%fittedCrack = 2.886e-11*N_cycles.^3-1.393e-7*N_cycles.^2+.000244*N_cycles+1.246;
%vpa N
%funct = 2.886e-11*N.^3-1.393e-7*N.^2+.000244*N+1.246;
%diffDADN = diff(funct)
%dadN2 = diff(cracklength)./diff(N_cycles);
dadN = 1.399e-4*exp(.0001224*N_cycles);
%dadN2 = 8.658e-11*N_cycles.^2-2.786e-7*N_cycles+.000244;
f2 = fit(cracklength',dadN','poly3')
figure();
plot(f2,cracklength,dadN)
title('Fitted curve of dadN vs crack length');
xlabel('Crack Length (in)');
ylabel('dadN')
%%
figure();
plot(log(N_cycles),log(cracklength))
%%
aW = [.2,.3,.4,.5,.6,.7,.8];
QaW = [2.070,2.090,1.836,1.634,1.46,1.351,1.314];
f3 = fit(aW',QaW','poly1')
plot(f3,aW,QaW)
%%
cracklength./3.25;
Q(cracklength./3.25);
deltaK = (1999.8./(0.5.*sqrt(3.25.*pi))).*((2+cracklength./3.25)./(1-cracklength./3.25).^(3/2)).*Q(cracklength./3.25);
%dadN_new = 5.432e-5*cracklength.^3-.0004397*cracklength.^2+.001122*cracklength-.0006817;
figure();
plot(deltaK,dadN)
title('deltaK vs da/dN')
ylabel('da/dN')
xlabel('Delta K (psi sqrt(in)')
figure();
plot(log(deltaK),log(dadN));
title('log(deltaK) vs log(da/dN)')
ylabel('log(da/dN)')
xlabel('log(Delta K) (psi sqrt(in)')
%%
f4 = fit(deltaK',dadN','power1');
% deltaK_new = linspace(18350,23690);
% dadN_new = 2.05e-10*deltaK_new.^(1.415);
% f5 = fit(deltaK(1:6)',dadN(1:6)','power1')
%%
Q(.01)
%%
fun = @(a) 1./(2.05*10.^-10.*(1.12.*8000*sqrt(pi.*a)).^1.415);
q = integral(fun,.01,2.2838)
%%
fun2 = @(a) 1./(10^-20.*((5*2932.9367*9*sqrt(pi*a))/(2*(9-2.5^2))).^4);
q2 = integral(fun2, .05,.5)
%% funciton

function Qfun = Q(a)
Qfun = -1.472*a + 2.415;
end




