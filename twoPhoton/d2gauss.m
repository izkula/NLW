function h = d2gauss(m,n,sig)
    theta = 0;
    std1 = sig;
    std2 = sig;
    n2 = m;
    n1 = m;
    r=[cos(theta) -sin(theta);
        sin(theta) cos(theta)];
    for i = 1 : n2
        for j = 1 : n1
            u = r * [j-(n1+1)/2 i-(n2+1)/2]';
            h(i,j) = gauss(u(1),std1)*gauss(u(2),std2);
        end
    end
    h = h / sqrt(sum(sum(h.*h)));
end

% Function "gauss.m":
function y = gauss(x,std)
    y = exp(-x^2/(2*std^2)) / (std*sqrt(2*pi));
end