function L = thetaVerfahren(R,In)
% description: DGL Solver mit theta-Verfahren: für theta = 1:
% impl. Euler, für theta = 0: expl. Euler, sonst: Mischformen
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
% In.theta ... Parameter in [0,1], siehe Blatt 4 Aufg. 3
% output:
% L.grid Zeitgitter (Zeitpunkte der berechneten Lösungen)
% L.x Matrix mit Lösungen x(t_i) in R^(dxm), Jede Spalte: ein x_i, i=0..m-1
% L.name string mit Name des Verfahrens
%
% author: Christian Winkler, Alexander Blech

m = length(In.grid); %Anzahl der Zeitschritte

%Den Fall theta = 0 kann man durch expl. Euler lösen, dann kann man sich
%das implizite Lösen sparen ...
if In.theta == 0
    In.BT = [0,0;0,1]; %Euler_expl
    L = explRK(R, In); %damit sind wir schon fertig ... 
else
    L.grid = In.grid;
    L.x = zeros(In.d,m); 
    L.x(:,1) = In.xstart;

    %%%%%%%%%%%beginne mit berechnung der lösungen%%%%%%%%%%%%%%%%%%%%%%%
    for it=1:1:m-1
        theta = In.theta; %Parameter im Verfahren, 'mischt' Euler expl. und impl.
        t0 = In.grid(it+1);
        h = In.grid(it+1)-In.grid(it); %breite Zeitschritt
        x0 = L.x(:,it);
        Phi = @(x) x0 + h*(1-theta)*R.F(t0,x0) + h*theta*R.F(t0+h,x) - x;
        % Berechnung der approx. Lösung am nächsten Gitterzeitpunkt mithilfe fsolve (sucht nullstelle von Phi)
        x1 = In.zerosolver(Phi,x0);
    %     x1 = newton(Phi,DPhi,x0,R.newton.eps_rel,R.newton.eps_abs,R.newton.maxIt);
        %Speichern in der Outputvariablen
        L.x(:,it+1) = x1;
    end
end
L.name = sprintf('Theta-Verfahren:theta=%.2e',In.theta);
end