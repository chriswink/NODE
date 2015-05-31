% description: Pendel script für Übungsblatt 3, NODE
%Untersuchung Stabilitätsverhalten verschiedener Verfahren am mathemat.
%Pendel

% author: Christian Winkler, Alex Blech. 
% mail: christian.winkler@uni.kn, alexander.blech@uni.kn

function pendel()
%PARAMETER
ZEITPLOT = true; PHASENPLOT = true; FSOLVE = false;
g   = 9.81;%Erdbeschl.
l   = 1.0; %Länge des Pendels
t0  = 0;
t1  = 10; %intervall, groß, da langzeitverhalten betrachtet wird

h = 1.5e-2; %Schrittweite
N = round((t1-t0)/h); %Anzahl Schritte
%rechte Seite
f   = @(t,x)[x(2);-g/l*x(1)];
x0  = [pi/20;0]; %Anfangswert: [Auslenkung(in Grad);Anfangswinkelgeschw.]

%SOLVER
%Jetzt alle infos in struct packen für die solver
BT_RK = [0,0,0,0,0;0.5,0.5,0,0,0;0.5,0,0.5,0,0;1,0,0,1,0;0,1/6.,1/3.,1/3.,1/6.]; % klass. RK Stufe 4
BT_Ee = [0,0;0,1]; %Euler_expl
In = struct('d',2,'xstart',x0,'grid',linspace(t0,t1,N),'BT',zeros(0,0));
if FSOLVE %Verfahren zur Nullstellenbestimmung
    In.zerosolver = @fsolve; 
else 
    maxIt = 10000;
    eps = 1.e-4;
    nullit = @(phi,x0) zeroIterate(phi,x0,maxIt,eps);
    In.zerosolver = nullit;
end
R.F = f; %Rechte Seite

%Diverse solver anwerfen
%explizite
In.BT = BT_RK;
L_RK = explRK(R,In);
L_RK.name = 'kl. Runge-Kutta';
In.BT = BT_Ee;
L_Ee = explRK(R,In);
L_Ee.name = 'Euler expl.';
%implizite
L_Ei = implEuler(R,In);
L_Mp = implMipu(R,In);




%PLOT
if ZEITPLOT
    plot_time_phi([L_Ei,L_RK,L_Ee,L_Mp]);
end
if PHASENPLOT
    figure(2);
    plot_phasenraum([L_Ei,L_RK,L_Mp,L_Ee]);
end
% 
end

function ax = plot_time_phi(L)
%Plotte Zeit gegen Auslenkung Phi des Pendels.
%L ...  Liste mit Lösungen aus verschiedenen Solvern: L=[L_1,...,L_n]
%return axes instance
hold on;
for j=1:length(L)
    h   = L(j).grid(2)-L(j).grid(1);
    ax  = plot(L(j).grid,L(j).x(1,:),'DisplayName',sprintf('%s, h:%.2e',L(j).name,h));
end
tit = sprintf('Math. Pendel, t=%.1e - %.1e',L(1).grid(1),L(1).grid(end));
title(tit);
xlabel('Zeit');
ylabel('Auslenkung in Bogenmaß');   
legend('Location','best');
hold off;
end

function ax = plot_phasenraum(L)
%Plotte Zeit gegen Auslenkung Phi des Pendels.
%L ...  Liste mit Lösungen aus n verschiedenen Solvern: L=[L_1,...,L_n]
%return axes instance
hold on;
for j=1:length(L)
    h   = L(j).grid(2)-L(j).grid(1);
    ax  = plot(L(j).x(1,:),L(j).x(2,:),'DisplayName',sprintf('%s, h:%.2e',L(j).name,h));
end
tit = sprintf('Math. Pendel, t=%.1f - %.1f, Phasenraum',L(1).grid(1),L(1).grid(end));
title(tit);
xlabel('$\varphi$','Interpreter','LaTex');
ylabel('$\dot\varphi(t)$','Interpreter','LaTex');   
legend('Location','best');
hold off;
end