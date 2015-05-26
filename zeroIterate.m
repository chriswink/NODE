function x = zeroIterate(f,x0,maxIt,eps)
% description: Fixpunktiteration zur Nulsstellenbestimmung
% 
% input:
%
% output:
%
%
% author: Christian Winkler, Alex Blech. christian.winkler@uni.kn, alexander.blech@uni.kn
foundRoot=false;
for i = 1:maxIt
    x = f(x0);
    if norm(f(x)) < eps
        foundRoot = true;
        break;
    else
        x0 = x;
    end
end
if foundRoot == false
    error('Keine Nullstelle gefunden nach %d Iter.: norm(f(x))=%lf', maxIt, norm(f(x)));
end
end