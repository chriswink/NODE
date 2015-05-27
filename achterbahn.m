% description: Achterbahn script für Übungsblatt 2, NODE

% author: Christian Winkler, Alex Blech. 
% mail: christian.winkler@uni.kn, alexander.blech@uni.kn

PLOT=true;
FSOLVE = false; %benutze fsolve zur Nullstellenbestimmung oder 
                %eigene Fixpunktiteration
m 	= 1;%Masse 1kg
gr 	= [0;-10];%Erdbeschleunigung

%Frage User nach Bahn
bahn = input('Wähle Bahn (1-5):');

%Definiere die Bahnen
if 		bahn==1
	g		= @(x) [x;0.8-0.8*x]; %Lineare Bahn
	dg		= @(x) [1;-0.8]; %1.Ableitung der Bahn
	ddg		= @(x) [0;0]; %2.Ableitung der Bahn
	nam		= 'Gerade';
elseif 	bahn==2
	g 		= @(x) [x;0.8*(x-1).^2]; %parabel ast mit scheitel am ziel
	dg		= @(x) [1;1.6*(x-1)]; %1.Ableitung der Bahn
	ddg		= @(x) [0;1.6]; %2.Ableitung der Bahn
	nam		= 'Parabel';
elseif 	bahn==3
	g 		= @(x) [x;4/3*x.^2-1.6*4/3*x+0.8]; %parabel ast mit scheitel links vom Ziel
	dg		= @(x) [1;8/3*x-1.6*4/3]; %1.Ableitung der Bahn
	ddg		= @(x) [0;8/3]; %2.Ableitung der Bahn
	nam		= 'Parabel verschoben';
elseif 	bahn==4
	g 		= @(x) [x;0.8-0.8*x-0.05*sin(8*pi*x)]; %wellenbahn. Diese ist schwer für die solver!
	dg		= @(x) [1;-0.8-0.4*pi*cos(8*pi*x)]; %1.Ableitung der Bahn
	ddg		= @(x) [0;3.2*pi*sin(8*pi*x)]; %2.Ableitung der Bahn
	nam		= 'Wellen';
elseif 	bahn==5
    %Diese Bahn ist nicht stetig diffbar bei 0, deshalb verschoben um 0.01
    %nach links
	g 		= @(x) [0.8*(1.35*pi/2*(x+0.01)-sin(1.35*pi/2*(x+0.01)));...
                    0.8 *(-1+cos((x+0.01)*pi/2)) + 0.8]; %Zykloid
    
    dg		= @(x) [0.8*1.35*pi/2-0.8*1.35*pi/2*cos(1.35*pi/2*(x+0.01));...
                    -sin((x+0.01)*pi/2)*0.8*pi/2]; %1.Ableitung der Bahn
	ddg		= @(x) [0.8*(1.35*pi/2)^2*sin(1.35*pi/2*(x+0.01));...
                    -0.8*(pi/2)^2*cos((x+0.01)*pi/2)]; %2.Ableitung der Bahn
	nam		= 'Zykloid';
end

%plotte Bahn
if PLOT
	x = linspace(0,1,50);
	b=g(x);
	plot(b(1,:),b(2,:));
	title(nam);
	xlabel('x');
	ylabel('y');
end

%Definiere 'Rechte Seite' -> AWP vom ÜB auf System 1. Ordnung umgewandelt: y(1):=sigma(t),y(2) = sigma'(t)
F = @(t,y) [y(2); (1./(norm(dg(y(1))))^2)*(m*gr'*dg(y(1)) - y(2)^2*(ddg(y(1)))'*dg(y(1)))];
%Definiere Startwert [s(0),y(0)]: s(0) ist der Ort zum Zeitpunkt 0,
%also Start der Bahn = g(0). y(0): Startgeschw. = 0
x0 = [0.0;0.0];
%  x0 = [0;0.5];

%Definiere Intervall zur Lsg und Anzahl Zeitschritte
t0	=0;
t1	=1;
N 	= 500;
% N= input('Anzahl Zeitschritte:');

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

R.F = F; %Rechte Seite

%Diverse solver anwerfen
%explizite
In.BT = BT_RK;
L_RK = explRK(R,In);
In.BT = BT_Ee;
L_Ee = explRK(R,In);
%implizite
L_Ei = implEuler(R,In);
L_Mp = implMipu(R,In);
%plot height vs. time für verschiedene solver
figure();
hold on;
p	 = plot_time_height(L_RK,nam,'RK expl', g,'-r');
p	 = plot_time_height(L_Ee,nam,'Euler expl', g, ':r');
p	 = plot_time_height(L_Ei,nam,'Euler impl', g, '.-g');
p	 = plot_time_height(L_Mp,nam,'impl. Mittelpunkt', g, '--g');
plot(gca(),L_RK.grid,zeros(length(L_RK.grid)),':b'); %Einen 'Boden' einzeichnen
h = (t1-t0)/N;
leg1 = sprintf('Solver:RK-4, Schrittweite:%.2e',h);
leg2 = sprintf('Solver:Euler explizit, Schrittweite:%.2e',h);
leg3 = sprintf('Solver:Euler impl, Schrittweite:%.2e',h);
leg4 = sprintf('Solver:impl. Mittelpunkt, Schrittweite:%.2e',h);

legend(leg1,leg2,leg3,leg4,'Boden','Location','southwest');
hold off;
