function rho = GetRho(rho,r,s )
incr = 2;
decr = 2;
delta = 10;

for i = 1:size(rho,1)
    if r(i) > delta*s(i)
        rho(i) = rho(i)*incr;
    elseif s(i) > delta*r(i)
        rho(i) = rho(i)/decr;
    else
        rho(i) = rho(i);
    end
end

end

