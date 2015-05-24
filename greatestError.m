function err = greatestError(sol,L)
%calculate greatest abs(sol(t_i)-L.x(t_i)), 
%sol is scalar function handle for exakt solution
%L is output from DGL solver


t = L.grid;
for i=1:1:length(t)
	x_ex(i) = sol(t(i));
end
%  plot(t,x_ex,t,L.x);
err = max(abs(L.x-x_ex));

end