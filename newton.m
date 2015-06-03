function x = newton(f, DF, x0)
% description: Newton Verfahren zur Nullstellensuche
% 
% input:
% f ... function handle, f:R^n->R^n
% Df ... function handle, returns Jacobian (nxn Matrix)
% x0... Startwert aus R^n
% eps_rel ... Toleranz für relative Verbesserung der Approximation pro Iteration
% eps_abs ... Toleranz für absolute Verbesserung der Approximation pro Iteration
% maxIt ... maximale Anzahl der Iterationen
%
% output:
% x ... approximierte Nullstelle der Fkt.
%
%author: Christian Winkler
maxIt = 1000;
x=1;
i=1;
while (i <= maxIt && norm(x,Inf) >= 10^(-6) && norm(x-x0)/norm(x0) >= 10^(-6))
	%Iteration
	A = DF(x0);
	b = -f(x0);
	z = A \ b;
	
	x = x0 + z;
	
	x0 = x;
	i = i+1;
end
if i == maxIt
    error('Max. Iterationen  (%d) erreicht, keine Nullstelle approximiert.',maxIt);
end
end