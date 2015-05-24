function x = newton(f, Df, x0, eps_rel, eps_abs, maxIt)
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
i=0;
while i<maxIt
	%Iteration
	A = DF(x0);
	b = -f(x0);
	z = A \ b;
	
	x = x0 + z;
	
	%Fehler betrachten
	err_abs = max(abs(x));
	err_rel = norm(x-x0)/norm(x0);
	if err_abs <= eps_abs and err_rel <= eps_rel
		break;
	end
	
	x0 = x;
	i = i+1;
end