function L = imp_mittelpunkt(R,I)

% description: implizite Mittelpunktregel
% input: 
% output:
% author: Alexander Blech

L.x(:,1) = I.xstart;
M=100;

for n = 1:length(I.grid)-1
    
    y(1) = L.x(:,n);
    for k = 1:M
        y(k+1) = L.x(:,n) + I.h * R.F(I.grid(n)+0.5*I.h,0.5*(y(k)+L.x(:,n)));
    end
    
    L.x(:,n+1) = y(k+1);...L.x(:,n) + I.h * R.F(0.5*(I.grid(n+1)+I.grid(n)),0.5*(y(k+1)+L.x(:,n)));
end

L.grid = I.grid;