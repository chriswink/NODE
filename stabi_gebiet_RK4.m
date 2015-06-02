%Plotte die Stabilitätsfunktion des klass. RK Verfahrens in der komplexen
%Ebene
%author: Christian Winkler. Alexander Blech
R = @(z) (1 + z + 0.5*z.^2 + 1/6*z.^3 + 1/24*z.^4);

N=100;

x=linspace(-4,2,N);
x = complex(x);
y = 1i * linspace(-3,3,N);

r = zeros(N,N);
for j=1:N
    for k=1:N
        r(j,k) = R(x(j)+y(k));
    end
end

%Für Aufgabe 2a):
g = @(z)(1+ z + 0.5.*z.^2 + 1/6 .* z.^3 + 1/24 .* z.^4);
x = linspace(0,2,100);
arg = x*i*sqrt(9.81/1);
plot(x,abs(g(arg)));

Z = abs(r);
image(real(x),imag(y),rot90(Z),'CDataMapping','scaled');

cbar = colorbar();
% title_handle = get(cbar,'Title');
ylabel(cbar,'Abs(R(z)');
set(gca, 'CLim', [0, 1]);
xlabel('Re(z)');
ylabel('Im(z)');
title('Stabilitätsgebiet RK4');