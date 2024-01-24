function PW = beta_expansion(N)
    n = log2(N);
    index2bin = zeros(N, n);
    for i = 0: N-1
        index2bin(i+1, :) = arrayfun(@(x) str2double(x), dec2bin(i, n));
    end
    W = zeros(1, n);
    beta = 2^(1/4);
    for i = 1: n
        W(i) = beta^(n-i);
    end
    PW = sum(index2bin .* W, 2);
end