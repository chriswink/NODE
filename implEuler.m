function L = implEuler(R,In)
% description: DGL Solver mit implizitem Euler Verfahren
% 
% input:
% R ... struct mit Infos zu rechter Seite der DGL
% R.F Funktion F(t,x) = x'; F: RxR^d -> R^d
% R.DF Jacobian
% In ... struct mit
% In.d Dimension des Systems
% In.xstart Startwert x0 in R^d
% In.grid diskr. Zeitgitter in R^(1xm): [t0,t1,...,t_(m-1)]
% In.newton...struct mit Parametern zum Newtonverfahren
% In.newton.eps_rel ... relative Fehlerschranke
% In.newton.eps_abs ... absolute Fehlerschranke
% In.newton.maxIt ... maximale Anzahl Iterationen im Newtonverfahren

%
% output:
% L.grid Zeitgitter (Zeitpunkte der berechneten Lösungen)
% L.x Matrix mit Lösungen x(t_i) in R^(dxm), Jede Spalte: ein x_i, i=0..m-1
%
% author: Christian Winkler

m = length(In.grid); %Anzahl der Zeitschritte

L.grid = In.grid;
L.x = zeros(In.d,m); 
L.x(:,1) = In.xstart;

%%%%%%%%%%%beginne mit berechnung der lösungen%%%%%%%%%%%%%%%%%%%%%%%
for it=1:1:m-1
	t1 = In.grid(it+1);
	h = In.grid(it+1)-In.grid(it); %breite Zeitschritt
	x0 = L.x(:,it);
	Phi = @(x) x0-x+h*R.F(t1,x);
%     DPhi = h*R.DF - eye('like',R.DF);
	% Berechnung der approx. Lösung am nächsten Gitterzeitpunkt mithilfe fsolve (sucht nullstelle von Phi)
	x1 = fsolve(Phi,x0);
%     x1 = newton(Phi,DPhi,x0,R.newton.eps_rel,R.newton.eps_abs,R.newton.maxIt);
	%Speichern in der Outputvariablen
	L.x(:,it+1) = x1;
end