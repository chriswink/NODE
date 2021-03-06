function L = gauss4(R,In)
% description: DGL Solver mit implizitem Gauss Verfahren Stufe 4
% 
% input:
% R ... struct mit Infos zu rechter Seite der DGL
% R.F Funktion F(t,x) = x'; F: RxR^d -> R^d
% R.dF ... Ableitung von F für Newton

% In ... struct mit
% In.d Dimension des Systems
% In.xstart Startwert x0 in R^d
% In.grid diskr. Zeitgitter in R^(1xm): [t0,t1,...,t_(m-1)]
% In.zerosolver function handle x = zerosolver(f,x0) mit f(x)=0
% In.newton ... sofern das attribut Newton existiert, wird die Nullstelle 
%               mit Newton gesucht
% output:
% L.grid Zeitgitter (Zeitpunkte der berechneten Lösungen)
% L.x Matrix mit Lösungen x(t_i) in R^(dxm), Jede Spalte: ein x_i, i=0..m-1
% L.name string mit Name des Verfahrens
%
% author: Christian Winkler, Alexander Blech

m = length(In.grid); %Anzahl der Zeitschritte

L.grid = In.grid;
L.x = zeros(In.d,m); 
L.x(:,1) = In.xstart;

%Die Werte aus dem Butcher Tableau:
a11 = 0.25; a12 = 0.25 - sqrt(3)/6; a21 = 0.25 + sqrt(3)/6; a22 = 0.25;
b1 = 0.5 - sqrt(3)/6; b2 = 0.5 + sqrt(3)/6;
c1 = 0.5; c2 = 0.5;

%%%%%%%%%%%beginne mit berechnung der lösungen%%%%%%%%%%%%%%%%%%%%%%%
for it=1:1:m-1
    t0 = In.grid(it);
	t1 = In.grid(it+1);
	h = In.grid(it+1)-In.grid(it); %breite Zeitschritt
	x0 = L.x(:,it);
    %Implizites Gleichungssystem für K1,K2:
	Phi = @(K) [R.F(t0+b1*h, x0+h*a11*K(1)+h*a12*K(2)) - K(1);...
                R.F(t0+b2*h, x0+h*a21*K(1)+h*a22*K(2)) - K(2)];
    K0 = [R.F(t0,x0);R.F(t0,x0)];%startwert für Iteratiion/Newton
	% Berechnung der approx. Lösung am nächsten Gitterzeitpunkt (suche nullstelle von Phi)
    if any(strcmp('newton',fieldnames(In)))
        %Definiere Jacobian
        d11Phi = @(z) h*a11*R.dF(t0+b1*h,x0+h*a11*z(1)+h*a12*z(2)) - 1;
        d12Phi = @(z) h*a12*R.dF(t0+b1*h,x0+h*a11*z(1)+h*a12*z(2));
        d22Phi = @(z) h*a22*R.dF(t0+b1*h,x0+h*a11*z(1)+h*a12*z(2)) - 1;
        d21Phi = @(z) h*a21*R.dF(t0+b1*h,x0+h*a11*z(1)+h*a12*z(2));
        dPhi = @(z)[d11Phi(z),d12Phi(z);...
                d21Phi(z),d22Phi(z)];
        K = newton(Phi,dPhi,K0);
    else
        K = In.zerosolver(Phi,K0);   
    end
    x1 = x0 + h*c1*K(1) + h*c2*K(2);
	%Speichern in der Outputvariablen
	L.x(:,it+1) = x1;
end
L.name = 'impl. Gauss s=4';
end