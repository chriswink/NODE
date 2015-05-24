function p = plot_time_height(L,bahn,solver,g)
% description: Plot für die Achterbahn Lösungen, Zeit vs. Höhe
% 
% input:
% L ... struct mit Lösung aus unserem DGL-Solver
% bahn ... string mit Name der Bahn
% solver ... string mit Name des Solvers
% g ... funktion der bahn (function handle)
%
% output:
% plot
%
% author: Christian Winkler, Alex Blech. christian.winkler@uni.kn, alexander.blech@uni.kn

ort=g(L.x(1,:)); %(x,y) Paare: Ort in R^2 zu Gitterzeitpunkten
ax=gca();
p = plot(ax,L.grid,ort(2,:),'p'); %plotte Zeit auf x-Achse und Höhe auf y-Achse
xlabel('Time');
ylabel('Height');
leg = sprintf('Bahn: %s, Solver:%s',bahn,solver);
legend(leg);

%checken, wo zum ersten Ziellinie überquert wird
n = find(ort(1,:) >=1, 1);
printf('\nZiellinie überquert bei t=%.2f\n',L.grid(n));
end