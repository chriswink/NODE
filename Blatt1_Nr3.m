% matlab-Script zu Aufgabe 3, b) des ersten Übungsblatts

% definiere input
clear all


a = 1;
R.F  = @(t,y) 2*t*(y+a)^2;
R.dF = @(t,y) 4*t*(y+a);

D.grid = [0.1,0.05,0.01,0.005,0.001];

for k=1:length(D.grid)

    I.d = 1;
    I.xstart = 1;
    I.h = 0.01;
    I.grid = 0:I.h:0.5/(sqrt(1+a));
    I.zerosolver = @fsolve;

    y = @(t) (a+1)./(1-(a+1)*t.^2) - a;
    for n = 1:length(I.grid)
        y_data = y(I.grid(n));
    end
    plot(I.grid,y_data)
    
    G4 = gauss4(R,I);
%     EE = exp_euler(R,I);
%     IE = imp_euler(R,I);
    IM = imp_mittelpunkt(R,I);
%     RK = runge_kutta_4(R,I);
% 
    D.G4(k) = norm(y_data-G4.x,Inf);
%     D.EE(k) = norm(y_data-EE.x,Inf);
%     D.IE(k) = norm(y_data-IE.x,Inf);
    D.IM(k) = norm(y_data-IM.x,Inf);
%     D.RK(k) = norm(y_data-RK.x,Inf);
end

plot(I.grid,y_data)...,G4.grid,G4.x)...,IM.grid,IM.x)...EE.grid,EE.x,IE.grid,IE.x,IM.grid,IM.x,RK.grid,RK.x)

D.plot(:,1) = polyfit(D.grid,D.G4,1)';
% D.plot(:,2) = polyfit(D.grid,D.EE,1)';
% D.plot(:,3) = polyfit(D.grid,D.IE,1)';
D.plot(:,4) = polyfit(D.grid,D.IM,1)';
% D.plot(:,5) = polyfit(D.grid,D.RK,1)';

D.plot