%matlab-Script zu Aufgabe 3, b) des ersten Uebungsblatts
clear all

%Wahl des Fixpunktverfahrens
NEWTON = true;

%Variable a und rechte Seite R.F festlegen
a = 0.2;
R.F  = @(t,y) 2*t*(y+a)^2;
R.dF = @(t,y) 4*t*(y+a);

%Gitter fuer den Fehlerplot: Schrittweiten I.h
D.grid = [0.05,0.02,0.01,0.005,0.001];

%Schleife ueber verschiedene Schrittweiten
for k=1:length(D.grid)
    
    %Input fuer numerische Routinen
    I.d = 1;
    I.xstart = 1;
    I.h = D.grid(k);
    I.grid = 0:I.h:0.5/(sqrt(1+a)); 
    if NEWTON
        I.newton = true;
    else
        maxIt = 10000;
        eps = 1.e-4;
        nullit = @(phi,x0) zeroIterate(phi,x0,maxIt,eps);
        I.zerosolver = nullit;
    end

    %aufrufen der Routinen
    G4 = gauss4(R,I);
%     EE = exp_euler(R,I);
%     IE = imp_euler(R,I);
    IM = implMipu(R,I); ...imp_mittelpunkt(R,I);
%     RK = runge_kutta_4(R,I);
    
    %Speichern der exakten Loesing in einen Vektor
    y = @(t) (a+1)./(1-(a+1)*t.^2) - a;
    for n = 1:length(I.grid)
        y_data(n) = y(I.grid(n));
    end
    
    %groessten Fehler berechnen
    D.G4(k) = norm(y_data-G4.x,Inf);
%     D.EE(k) = norm(y_data-EE.x,Inf);
%     D.IE(k) = norm(y_data-IE.x,Inf);
    D.IM(k) = norm(y_data-IM.x,Inf);
%     D.RK(k) = norm(y_data-RK.x,Inf);
end

%Ausgabe der exakten und approximierten Loesungen
% plot(I.grid,y_data,EE.grid,EE.x,IE.grid,IE.x,IM.grid,IM.x,RK.grid,RK.x)

%Fit zur Bestimmung der Konvergenzordnung
D.plot(:,1) = polyfit(log(D.grid)/log(10),log(D.G4)/log(10),1)';
% D.plot(:,2) = polyfit(log(D.grid)/log(10),log(D.EE)/log(10),1)';
% D.plot(:,3) = polyfit(log(D.grid)/log(10),log(D.IE)/log(10),1)';
D.plot(:,4) = polyfit(log(D.grid)/log(10),log(D.IM)/log(10),1)';
% D.plot(:,5) = polyfit(log(D.grid)/log(10),log(D.RK)/log(10),1)';

%logarithmischer Plot der Fehler
plot(log(D.grid)/log(10),log(D.G4)/log(10),log(D.grid)/log(10),log(D.IM)/log(10)) ...log(D.grid)/log(10),log(D.EE)/log(10),log(D.grid)/log(10),log(D.IE)/log(10),log(D.grid)/log(10),log(D.IM)/log(10),log(D.grid)/log(10),log(D.RK)/log(10))

%Ausgabe der numerischen Konvergenzordnungen in erster Zeile von D.plot
D.plot