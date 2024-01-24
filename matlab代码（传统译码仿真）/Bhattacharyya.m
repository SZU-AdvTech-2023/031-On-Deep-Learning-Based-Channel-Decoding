function Z = Bhattacharyya(N, ZW)
% 码长N，信道巴氏参数Z_init
    Z = [ZW];
    n = log2(N);
    for j = 0:n-1
        index = 2^j;
        % 模仿队列的形式，弹出队首，代入递归公式计算出两个，压入队列
        % 还可以不用列表[]实现，用列向量就行，参照于永润讲义中的代码
        for i = 1:index
            % 一生二，二生四，四生八，...
            % 不用bitreserve的实现
            Z1 = 2*Z(1)-(Z(1))^2;
            Z2 = (Z(1))^2;
            Z = [Z(2:end), Z1, Z2];

%             % 也可以这样实现
%             temp = Z(i);
%             Z(i) = 2*temp-temp^2;
%             Z(i+index) = temp^2;
%             % 最后还需要进行一次bitreserve
%             % Z = bit_reserve(Z);
        end
    end
    Z = Z';
end