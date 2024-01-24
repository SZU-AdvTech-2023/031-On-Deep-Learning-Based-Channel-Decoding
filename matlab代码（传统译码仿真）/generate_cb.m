function codebook = generate_cb(is_info, G)
    N = length(is_info);
    K = sum(is_info == 1);
    U = zeros(2^K, N);
    Info = zeros(2^K, K);
    for i = 0 : (2^K - 1)
        Info(i + 1, :) = arrayfun(@(x) str2double(x), dec2bin(i, K));
    end
    U(:, is_info == 1) = Info;
    codebook = mod(U * G, 2);
end