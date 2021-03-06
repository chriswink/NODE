function L = implEuler(R,In)
% description: DGL Solver mit implizitem Euler Verfahren
% 
% input:
% R ... struct mit Infos zu rechter Seite der DGL
% R.F Funktion F(t,x) = x'; F: RxR^d -> R^d
%
% In ... struct mit
% In.d Dimension des Systems
% In.xstart Startwert x0 in R^d
% In.grid diskr. Zeitgitter in R^(1xm): [t0,t1,...,t_(m-1)]
% In.zerosolver function handle x = zerosolver(f,x0) mit f(x)=0
% output:
% L.grid Zeitgitter (Zeitpunkte der berechneten Lösungen)
% L.x Matrix mit Lösungen x(t_i) in R^(dxm), Jede Spalte: ein x_i, i=0..m-1
% L.name string mit Name des Verfahrens
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
	Phi = @(x) x0 + h*R.F(t1,x) - x;
	% Berechnung der approx. Lösung am nächsten Gitterzeitpunkt mithilfe fsolve (sucht nullstelle von Phi)
	x1 = In.zerosolver(Phi,x0);
%     x1 = newton(Phi,DPhi,x0,R.newton.eps_rel,R.newton.eps_abs,R.newton.maxIt);
	%Speichern in der Outputvariablen
	L.x(:,it+1) = x1;
end
L.name = 'Euler impl.';