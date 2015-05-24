PLOT=true;
bahn = 1;
m 	= 1;%Masse 1kg
gr 	= [0;-10];%Erdbeschleunigung

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
	g 		= @(x) [x;0.8-0.8*x-0.1*sin(8*pi*x)]; %wellenbahn. Diese ist schwer für die solver!
	dg		= @(x) [1;-0.8-0.8*pi*cos(8*pi*x)]; %1.Ableitung der Bahn
	ddg		= @(x) [0;6.4*pi*sin(8*pi*x)]; %2.Ableitung der Bahn
	nam		= 'Wellen';
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
F = @(y) [y(2); (1./(norm(dg(y(1))))^2)*(m*gr'*dg(y(1)) - y(2)^2*(ddg(y(1)))'*dg(y(1)))];

%Definiere Startwert [s(0),y(0)]: s(0) ist der Ort zum Zeitpunkt 0, also Start der Bahn = g(0). y(0): Startgeschw. = 0
x0 = [0;0];
%  x0 = [0;0.5];

%Definiere Intervall zur Lsg und Anzahl Zeitschritte
t0	=0;
t1	=1.2;
m 	=100;

%Jetzt alle infos in struct packen für den solver
BT_RK = [0,0,0,0,0;0.5,0.5,0,0,0;0.5,0,0.5,0,0;1,0,0,1,0;0,1/6.,1/3.,1/3.,1/6.]; % klass. RK Stufe 4
BT_Ee = [0,0;0,1]; %Euler_expl
In = struct('d',2,'xstart',x0,'grid',linspace(t0,t1,m),'BT',zeros(0,0));
R.F = F;
%Diverse solver anwerfen
In.BT = BT_RK;
L_RK = explRK(R,In);
In.BT = BT_Ee;
L_Ee = explRK(R,In);

%plot height vs. time für verschiedene solver
figure();
hold on;
p	 = plot_time_height(L_RK,nam,'sol', g);
p	 = plot_time_height(L_Ee,nam,'expl', g);

plot(gca(),L_RK.grid,zeros(length(L_RK.grid))); %Einen 'Boden' einzeichnen
leg1 = sprintf('Bahn: %s, Solver:RK-4',nam);
leg2 = sprintf('Bahn: %s, Solver:Euler expl.',nam);
legend(leg1,leg2,'Boden');
hold off;
