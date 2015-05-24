PLOT=false;
m = 1;%Masse 1kg
g = [0;-10];%Erdbeschleunigung
%Definiere die Bahnen
g1	= @(x) [x;0.8-0.8*x]; %Lineare Bahn
dg1	= @(x) [1;-0.8]; %1.Ableitung der Bahn
ddg1= @(x) [0;0]; %2.Ableitung der Bahn
%  g2 = @(x) [x.;0.8*(x.-1)^2];
%  g3 = @(x) [x.;0.8-0.8*x];

%plotte Bahnen
if PLOT
x = linspace(0,1,10);
b=g1(x);
plot(b(1,:),b(2,:));
legend('BAHN');
xlabel('x');
ylabel('y');
end

%Definiere 'Rechte Seite' -> AWP vom ÜB auf System 1. Ordnung umgewandelt: y(1):=sigma(t),y(2) = sigma'(t)
F = @(y) [y(2); (1./(norm(dg1(y(1))))^2)*(m*g'*dg1(y(1)) - y(2)^2*(ddg1(y(1)))'*dg1(y(1)))];
%Definiere Startwert [s(0),y(0)]: s(0) ist der Ort zum Zeitpunkt 0, also Start der Bahn = g1(0). y(0): Startgeschw. = 0
x0 = [0;0];
%Definiere Intervall zur Lsg und Anzahl Zeitschritte
t0	=0;
t1	=0.7;
m 	=50;
%Jetzt alle infos in struct packen
BT_RK = [0,0,0,0,0;0.5,0.5,0,0,0;0.5,0,0.5,0,0;1,0,0,1,0;0,1/6.,1/3.,1/3.,1/6.]; % klass. RK Stufe 4
In = struct('d',2,'xstart',x0,'grid',linspace(t0,t1,m),'BT',zeros(0,0));
R.F = F;
In.BT = BT_RK;
L_RK = explRK(R,In);