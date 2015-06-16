function k = Aufgabe3_blatt4(theta)
%matlab-Script zu Aufgabe 3 von Blatt 4: Bestimmung der
%numerischen Konvergenzordnung des theta-Verfahrensm für unterschiedl.
%Werte von theta
%input: theta, float (sinnvollerweise im Bereich zwische ca. 0 und 1, kann aber
%auch andere Werte annehmen
%output: numerische Konvergenzordnung
%authors: Alexander Blech, Christian Winkler
%emails: Alexander.Blech@uni.kn, Christian.Winkler@uni.kn

%Wahl des Fixpunktverfahrens
NEWTON = false;

%Variable a und rechte Seite R.F festlegen
a = 0.2;
R.F  = @(t,y) 2*t*(y+a)^2;
R.dF = @(t,y) 4*t*(y+a);

%Gitter fuer den Fehlerplot: Schrittweiten I.h
D.grid = [0.05,0.02,0.01,0.005,0.001,0.0005];

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
        maxIt = 100000;
        eps = 1.e-12;
        nullit = @(phi,x0) zeroIterate(phi,x0,maxIt,eps);
        I.zerosolver = nullit;
    end

    %aufrufen der Routinen
    I.theta = theta;
    L_0 = thetaVerfahren(R,I);
      
    
    %Speichern der exakten Loesing in einen Vektor
    y = @(t) (a+1)./(1-(a+1)*t.^2) - a;
    y_data = zeros(1,length(I.grid));
    for n = 1:length(I.grid)
        y_data(n) = y(I.grid(n));
    end
    
    %groessten Fehler berechnen
    D.L_0(k) = norm(y_data-L_0.x,Inf);
end

%Fit zur Bestimmung der Konvergenzordnung
D.plot(:,1) = polyfit(log(D.grid)/log(10),log(D.L_0)/log(10),1)';
%logarithmischer Plot der Fehler
plot(log(D.grid)/log(10),log(D.L_0)/log(10));

%Ausgabe der numerischen Konvergenzordnungen in erster Zeile von D.plot
fprintf('%s:Konvergenzordnung:%.3f \n',L_0.name, D.plot(1));
k =  D.plot(1);
end