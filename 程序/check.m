AAW = zeros(1, 100);
AAW(1) = W;
CP = 299.6;
g = zeros(1, 200);

for i = 1:99
    g(2 * i - 1) = AAW(i) - W;
    g(2 * i) = -AAW(i);
    P_average = (sol_best(i) + sol_best(i + 1)) / 2;

    if P_average >= CP
        AAW(i + 1) = AAW(i) - (P_average - CP) * (t0(i + 1) - t0(i));
    elseif P_average < CP
        taoW = 546 * exp(-0.01 * (CP - P_average)) + 316;
        AAW(i + 1) = W - (W - AAW(i)) * exp(- (t0(i + 1) - t0(i)) / taoW);
    end

end

g(199) = AAW(100) - W;
g(200) = -AAW(100);
plot(x0, AAW, '-o')
