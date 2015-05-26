%Teste  numerische Konvergenz einiger Verfahren
%anhand approx. Lsg von y'=2t(y-a)^2, a>0

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
%  for m=1:10:500 %Anzahl Zeitschritte
%  	printf('\rAnz. Zeitschritte: %d',m);
%  	In = struct('d',1,'xstart',y0,'grid',linspace(t0,t1,m),'BT',zeros(1,1));
%  	R.F = f;
%  	
%  	In.BT = BT_Ee;
%  	L = explRK(R,In);
%  	err_Ee(im) = greatestError(sol,L);
%  	In.BT = BT_RK;
%  	L = explRK(R,In);
%  	err_RK(im) = greatestError(sol,L);
%  	M(im) = m;
%  	im = im+1;
%  end
for m=5:5:25 %Anzahl Zeitschritte
	printf('\rAnz. Zeitschritte: %d',m);
	In = struct('d',1,'xstart',y0,'grid',linspace(t0,t1,m));
	R.F = f;
	L = implEuler(R,In);
	err_Ei(im) = greatestError(sol,L);
	M(im) = m;
	im = im+1;
end

%Plot (logarithmisch)