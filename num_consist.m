%Teste  numerische Konvergenz einiger Verfahren
%anhand approx. Lsg von y'=2t(y-a)^2, a>0

%Optionen
FSOLVE = true;
a = 0.1;%Parameter in DGL
y0 = 1;%Startwert
t0 = 0;
t1 = 0.5/sqrt(1+a); %Intervall Anfang und Ende
f = @(t,x) 2*t*(x+a).^2; %rechte Seite
sol = @(t) 1./(1/(1+a)-t.^2) - a;

%Definition der Butcher Tableaus:
BT_Ee = [0,0;0,1]; %Euler_expl
BT_RK = [0,0,0,0,0;0.5,0.5,0,0,0;0.5,0,0.5,0,0;1,0,0,1,0;0,1/6.,1/3.,1/3.,1/6.]; % klass. RK Stufe 4
In = struct('d',1,'xstart',y0,'grid',linspace(t0,t1,5),'BT',zeros(0,0));
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
for im = 1:1:8
    N = 2^(im+2); %anzahl Zeitschritte
    fprintf('\rAnz. Zeitschritte: %d',N);
    In.grid = linspace(t0,t1,N);
    In.BT = BT_RK;
    L_RK = explRK(R,In);
    L_RK.name = 'kl. Runge-Kutta';
    In.BT = BT_Ee;
    L_Ee = explRK(R,In);
    L_Ee.name = 'Euler expl.';
    %implizite
    L_Ei = implEuler(R,In);
    L_Mp = implMipu(R,In);
    err_Ei(im) = greatestError(sol,L_Ei);
    err_Mp(im) = greatestError(sol,L_Mp);
    err_RK(im) = greatestError(sol,L_RK);
    err_Ee(im) = greatestError(sol,L_Ee);
	M(im) = N;
end
%Plot (logarithmisch)
lM = log10(M);
l_Ee=log10(err_Ee);
l_Ei=log10(err_Ei);
l_Mp=log10(err_Mp);
l_RK=log10(err_RK);
 
[r_Ei,m,b] = regression(lM,l_Ei);
fprintf('\nKonvergenzordnung %s:%.4f',L_Ei.name,abs(m));
[r_Ee,m,b] = regression(lM,l_Ee);
fprintf('\nKonvergenzordnung %s:%.4f',L_Ee.name,abs(m));
[r_RK,m,b] = regression(lM,l_RK);
fprintf('\nKonvergenzordnung %s:%.4f',L_RK.name,abs(m));
[r_Mp,m,b] = regression(lM,l_Mp);
fprintf('\nKonvergenzordnung %s:%.4f',L_Mp.name,abs(m));

 %  plotregression(lM,l_Ei);
%  plotregression(lM,l_Ee);
% plotregression(lM,l_RK);
%  plotregression(lM,l_Mp);
 