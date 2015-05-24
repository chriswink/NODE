function L = explRK(R,In)
% description: DGL Solver mit expliziten RungeKutta Verfahren
% 
% input:
% R ... struct mit Infos zu rechter Seite der DGL
% R.F Funktion F(t,x) = x'; F: RxR^d -> R^d
% In ... struct mit
% In.d Dimension des Systems
% In.xstart Startwert x0 in R^d
% In.grid diskr. Zeitgitter in R^(1xm): [t0,t1,...,t_(m-1)]
% In.BT Butcher Tableau des s-stufigen expliziten Verfahrens in R^(s+1xs+1) 
%      (entsprechend  http://de.wikipedia.org/wiki/Runge-Kutta-Verfahren#Butcher-Tableau)
%
% output:
% L.grid Zeitgitter (Zeitpunkte der berechneten Lösungen)
% L.x Matrix mit Lösungen x(t_i) in R^(dxm), Jede Spalte: ein x_i, i=0..m-1
%
% author: Christian Winkler
[s,s] = size(In.BT); s = s-1;%Verfahrensstufe
m = length(In.grid); %Anzahl der Zeitschritte
%teste, ob Butcher Tableau tatsächlich explizites Verfahren beschreibt
if triu(In.BT(1:s,2:s+1)) ~= zeros(s,s)
	error('Geg. Butcher Tabl. kein expl. Verfahren! (Koeff nicht strikte untere Dreiecksmat.)');
end

L.grid = In.grid;
L.x = zeros(In.d,m); 
L.x(:,1) = In.xstart;

%Lese koeefizienten des butcher tableau in neue struct BT ein
BT.c = In.BT(1:s,1);
BT.A = In.BT(1:s,2:s+1);
BT.b = In.BT(s+1,2:s+1);
%%%%%%%%%%%beginne mit berechnung der lösungen%%%%%%%%%%%%%%%%%%%%%%%
for it=1:1:m-1
	t0 = In.grid(it);
	h = In.grid(it+1)-In.grid(it); %breite Zeitschritt
	x0 = L.x(:,it);
	
	K = zeros(length(x0),s); %Platz für die Zwischenrechnungen der K^(j)=F(y0+sum_i^j(h*A(j,i)*K^(j-1)))

	%jetzt berechne die K(j)
	for j=1:1:s
		% erst das Argument (für bessere Übersicht), dann: K(j) = F(arg)
		arg = x0;
		for i=1:1:j-1
			arg = arg + h * BT.A(j,i) * K(1:end,j-1);
		end
		K(1:end,j) = R.F(arg);
	end
	% Jetzt die Verfahrensfunktion
	V = zeros(length(x0),1);
	for j=1:1:s
		V = V + BT.b(j) * K(1:end,j);
	end
	
	% Berechnung der approx. Lösung am nächsten Gitterzeitpunkt
	x1 = x0 + h*V;
	%Speichern in der Outputvariablen
	L.x(:,it+1) = x1;
end
