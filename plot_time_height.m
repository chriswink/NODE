function p = plot_time_height(L,bahn,solver,g,linestyle)
% description: Plot für die Achterbahn Lösungen, Zeit vs. Höhe
% 
% input:
% L ... struct mit Lösung aus unserem DGL-Solver
% bahn ... string mit Name der Bahn
% solver ... string mit Name des Solvers
% g ... funktion der bahn (function handle)
% linestyle... str mit linestyle specifier,z.B.'--', '.', '.-'

% output:
% plot
%
% author: Christian Winkler, Alex Blech. christian.winkler@uni.kn, alexander.blech@uni.kn

ort=g(L.x(1,:)); %(x,y) Paare: Ort in R^2 zu Gitterzeitpunkten
h = L.grid(2)-L.grid(1); %Schrittweite
ax=gca();
p = plot(ax,L.grid,ort(2,:),linestyle); %plotte Zeit auf x-Achse und Höhe auf y-Achse
xlabel('Time');
ylabel('Height');
leg = sprintf('Solver:%s Schrittweite:%.2f',solver,h);
legend(leg);
titel = sprintf('Bahn: %s',bahn);
title(titel);

%checken, wo zum ersten Ziellinie überquert wird
n = find(ort(1,:) >=1, 1);
fprintf('\n %s :Ziellinie überquert bei t=%.5f\n',solver,L.grid(n));
end