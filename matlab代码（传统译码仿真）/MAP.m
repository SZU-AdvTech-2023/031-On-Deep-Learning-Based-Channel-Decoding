function u = MAP(codebook, y, G)
    num = length(codebook);
    s = 1 - 2 * codebook;
    mse = zeros(num, 1);
    for i = 1: num
        mse(i) = sum((s(i, :) - y).^2);
    end
    [~, index] = min(mse);
    u = mod(codebook(index, :)  * G, 2);
end
