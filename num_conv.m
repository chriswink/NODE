%Teste  numerische Konvergenz einiger Verfahren
%anhand approx. Lsg von y'=2t(y-a)^2, a>0
clear all
a = 0.1;%Parameter in DGL
y0 = 1;%Startwert
t0 = 0;
t1 = 0.5/sqrt(1+a); %Intervall Anfang und Ende
f = @(t,x) 2*t*(x-a)^2; %rechte Seite
sol = @(t) 1./(1/(1+a)-t^2) - a;
%Definition der Butcher Tableaus:
BT_Ee = [0,0;0,1]; %Euler_expl
BT_Ei = [1,1;0,1];%Euler_impl
BT_Mi = [0.5,0.5;0,1]; % Mittelpkt. implizit
BT_RK = [0,0,0,0,0;0.5,0.5,0,0,0;0.5,0,0.5,0,0;1,0,0,1,0;0,1/6.,1/3.,1/3.,1/6.]; % klass. RK Stufe 4

%Berechne approx. Lösungen
im = 1;
printf('\n');
%  for m=1:50:1000 %Anzahl Zeitschritte
for i=0:1:3%Anzahl Zeitschritte
	m = 10^i;
	printf('\n expl: errAnz. Zeitschritte: %d',m);
	In = struct('d',1,'xstart',y0,'grid',linspace(t0,t1,m),'BT',zeros(1,1));
	R.F = f;
	
	In.BT = BT_Ee;
	L_Ee = explRK(R,In);
	err_Ee(im) = greatestError(sol,L_Ee);
	In.BT = BT_RK;
	L_RK = explRK(R,In);
	err_RK(im) = greatestError(sol,L_RK);
	M(im) = m;
	im = im+1;
end
im=1;
for i=0:1:3%Anzahl Zeitschritte
	m = 10^i;
	printf('\n impl: Anz. Zeitschritte: %d',m);
	In = struct('d',1,'xstart',y0,'grid',linspace(t0,t1,m));
	R.F = f;
	L_Ei = implEuler(R,In);
	err_Ei(im) = greatestError(sol,L_Ei);
	L_Mi = implMipu(R,In);
	err_Mi(im) = greatestError(sol,L_Mi);
	M(im) = m;
	im = im+1;
end
%Plot der Lösungen
%  hold on;
%  t=linspace(t0,t1,100);
%  for i=1:length(t)
%  x(i) = sol(t(i));
%  end
%  for i=1:5
%  x_5(i) = sol(L_Ei.grid(i));
%  end
%  x_5(i)
%  plot(t,x);
%  plot(L_Ei.grid,abs(L_Ei.x-x_5),'r');
%  plot(L_Mi.grid,abs(L_Mi.x-x_5),'g');
%  plot(L_Ee.grid,abs(L_Ee.x-x_5),'y');
%  plot(L_RK.grid,abs(L_RK.x-x_5),'c');
%  legend('exakte Loesung','impl. Euler','impl. Mitt.','expl. Euler','kl. RK');
%  xlim([0,0.6]);
%  ylim([0.9,1.5]);
%  xlabel('time');
%  ylabel('y(t)');
%  print -djpg ddd.jpg
%  hold off;

%Bestimme Regressionsgeraden
lM = log10(M);
l_Ee=log10(err_Ee);
l_Ei=log10(err_Ei);
l_Mi=log10(err_Mi);
l_RK=log10(err_RK);
%  [r_Ee,m,b] = regression(lM,l_Ee)
%  [r_Ee,m,b] = regression(lM,l_Ee)
%  [r_Ee,m,b] = regression(lM,l_Ee)
%  [r_Ee,m,b] = regression(lM,l_Ee)

