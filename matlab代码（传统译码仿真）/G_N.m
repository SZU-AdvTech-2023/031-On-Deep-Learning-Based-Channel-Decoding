function G = G_N(n)
    G = [1 0; 1 1]; % n=1,N=2的生成矩阵G2
    for i = 1:n-1
        % G_N迭代生成G_(2N)
        N = 2^i; % 当前生成矩阵的行列数
        G_up = [];
        G_down = [];
        for j = 1:N
            % 例如，N=2，G2=[c1,c2]
            % 则,G_up=[c1,0,c2,0],G_donw=[c1,c1,c2,c2]
            % 有 G=[G_up; G_down]=[c1,0,c2,0
            %                      c1,c1,c2,c2]
            G_up = [G_up, G(:, j), zeros(N, 1)];
            G_down = [G_down, G(:, j), G(:, j)];         
        end
        G = [G_up; G_down];
    end
end