function u_best = SCL_decoder(LLR, L, is_info)
% 信道LLR，列表大小L，is_info为1代表信息位
    N = length(LLR);
    n = log2(N);
    division = 2.^(0:n); % 用于将P和C数组分割的向量
    lazy_copy = zeros(n, L);
    P = zeros(N-1, L);
    C = zeros(2*N-1, 2*L); % 第l个译码器对应第2l-1和2l两列
    u = zeros(N, L);
    PM = zeros(L, 1);
    activepath = zeros(L, 1);
    activepath(1) = 1; % 开始时仅译码器1运行
    lazy_copy(:, 1) = 1;
    stack = 2:L;
    stack_size = L-1;

    for i = 1: N
        % 译码ui
        LLR_layer = get_LLR_layer(i);
        is_odd = mod(i, 2);

        % 对每一个译码器进行P矩阵的计算，注意用到Lazy Copy的地方
        for l_index = 1: L
            if activepath(l_index) == 0
                continue;
            end
            switch i
                case 1
                    index1 = division(n); % N/2
                    for j = 1: index1
                        P(index1-1+j, l_index) = f(LLR(2*j-1), LLR(2*j));
                    end
                    for j = n-1:-1:1
                        index1 = division(j);
                        index2 = division(j+1);
                        for k = index1: index2-1
                            P(k, l_index) = f(P(2*k, l_index), P(2*k+1, l_index));
                        end
                    end
                case N/2+1
                    index1 = division(n);
                    for j = 1: index1
                        x_tmp = C(index1-1+j, 2*l_index-1);
                        P(index1-1+j, l_index) = g(LLR(2*j-1), LLR(2*j), x_tmp);
                    end
                    for j = n-1:-1:1
                        index1 = division(j);
                        index2 = division(j+1);
                        for k = index1: index2-1
                            P(k, l_index) = f(P(2*k, l_index), P(2*k+1, l_index));
                        end
                    end
                otherwise
                    index1 = division(LLR_layer + 1);
                    index2 = division(LLR_layer + 2);
                    for j = index1: index2-1
                        P(j, l_index) = g(P(2*j, lazy_copy(LLR_layer+2, l_index)),...
                            P(2*j+1, lazy_copy(LLR_layer+2, l_index)), C(j, 2*l_index-1));
                    end       
                    for level = LLR_layer:-1:1
                        index1 = division(level);
                        index2 = division(level + 1);
                        for j = index1: index2-1
                            P(j, l_index) = f(P(2*j, l_index), P(2*j+1, l_index));
                        end
                    end
            end
        end

        if is_info(i) == 1
            % 信息位硬判决译码
            PM_pair = realmax * ones(2, L); % PM选小，初始值大避免被选中
            for l_index = 1: L
                if activepath(l_index) == 0
                    continue;
                end
                % 默认译码器选择ui=0，复制的选择ui=1，即第一行为译0的结果，第二行为译1的结果
                if P(1, l_index) >= 0
                    PM_pair(1, l_index) = PM(l_index);
                    PM_pair(2, l_index) = PM(l_index) + P(1, l_index); % 不符合硬判决，增加惩罚项
                else
                    PM_pair(1, l_index) = PM(l_index) - P(1, l_index); % 不符合硬判决，增加惩罚项
                    PM_pair(2, l_index) = PM(l_index);
                end
            end
            select_count = min(2*sum(activepath), L); % 要选择的译码路径个数
            [~, PM_index] = sort(PM_pair(:)); % 对PM_pair所有PM值进行升序，得到顺序PM_index
            compare = zeros(2, L);
            for j = 1: select_count
                row = 2 - mod(PM_index(j), 2);
                col = floor((PM_index(j) + 1) / 2);
                compare(row, col) = 1; % 对应的PM值较小，选择并保存下来
            end
           
            % 选择完译码路径后，每一列存在三种情况：

            % 一列都是0，如果译码器为激活的，那么译码器死亡，压入栈中
            for j = 1: L
                if (activepath(j) == 1) && (compare(1, j) == 0) && (compare(2, j) == 0)
                    % 译码器SC_j死亡
                    activepath(j) = 0;
                    stack_size = stack_size + 1;
                    stack(stack_size) = j;
                end
            end

            % 一列中仅有一个1，更新对应的u，C,PM即可
            % 一列中有两个1，需要激活栈顶中的译码器
            for l_index = 1: L
                if activepath(l_index) == 0
                    continue;
                end
                path_state = compare(1, l_index)*2 + compare(2, l_index);
                switch path_state
                    case 1
                        % 列中仅一个1，且选中的是第二行，u=1
                        u(i, l_index) = 1;
                        C(1, 2*l_index-is_odd) = 1;
                        PM(l_index) = PM_pair(2, l_index);
                    case 2
                        % 列中仅一个1，且选中的是第一行，u=0
                        u(i, l_index) = 0;
                        C(1, 2*l_index-is_odd) = 0;
                        PM(l_index) = PM_pair(1, l_index);
                    case 3
                        % 一列中有两个1，需要激活栈顶中的译码器，来复制当前译码器l_index
                        index = stack(stack_size);
                        stack_size = stack_size - 1;
                        activepath(index) = 1; % 激活新的译码器
                        lazy_copy(:, index) = lazy_copy(:, l_index); % 新激活的译码器index复制当前译码器l_index
                        u(:, index) = u(:, l_index); % 译码结果复制
                        u(i, l_index) = 0; % 当前译为0
                        u(i, index) = 1; % 新复制的译为1
                        C(1, 2*l_index-is_odd) = 0;
                        C(1, 2*index-is_odd) = 1;
                        PM(l_index) = PM_pair(1, l_index);
                        PM(index) = PM_pair(2, l_index);
                end
            end
        else
            % 冻结比特译码
            for l_index = 1: L
                if activepath(l_index) == 0
                    continue;
                end
                if P(1, l_index) < 0
                    PM(l_index) = PM(l_index) - P(1, l_index);
                end
                if is_odd == 1
                    C(1, 2*l_index-1) = 0;
                else
                    C(1, 2*l_index) = 0;
                end
            end
        end
        
        % 更新矩阵C
        for l_index = 1: L
            if activepath(l_index) == 0
                continue;
            end
            if is_odd == 0
                % i为偶数，需要进行位计算
                bit_layer = get_bit_layer(i);
                for j = 1:bit_layer
                    index1 = division(j);
                    index2 = division(j+1);
                    for k = index1: index2-1
                        C(2*k, 2*l_index) = mod(C(k, 2*lazy_copy(j, l_index)-1)...
                                                  + C(k, 2*l_index), 2);
                        C(2*k+1, 2*l_index) = C(k, 2*l_index);
                    end
                end
                index1 = division(bit_layer + 1);
                index2 = division(bit_layer + 2);
                for j = index1: index2-1
                    C(2*j, 2*l_index-1) = mod(C(j, 2*lazy_copy(bit_layer+1, l_index)-1)...
                                                + C(j, 2*l_index), 2);
                    C(2*j+1, 2*l_index-1) = C(j, 2*l_index);
                end
            end
        end
        % 更新lazy_copy
        if i < N
            for j = 1: get_LLR_layer(i+1)+1
                for l_index = 1: L
                    if activepath(l_index)
                        lazy_copy(j, l_index) = l_index;
                    end
                end
            end
        end
    end

    [~, min_index] = min(PM);
    u_best = u(:, min_index)';

end

function layer = get_LLR_layer(i)
% 对于ui而言，计算i-1对应的二进制数的低位中连续0的个数
% 即，该树节点能够连续向右上方走的最大步数，也是之后进行f运算的次数
    x = i-1;
    layer = 0;
    if x ~= 0
        while (mod(x, 2) == 0)
            x = x/2;
            layer = layer + 1;
        end
    end
end

function layer = get_bit_layer(i)
% 对于ui而言，计算i-1对应的二进制数的低位中连续1的个数再减去1
% 低位中连续1的个数代表beta返回计算的层数
% 其中前i-1层是将运算结果放在第二列（都是作为右子树），最后一层是放在第一列（作为左子树的返回）
% 这规律是容易理解的，1代表是向右，自然是作为右子树，假设最后一个1（低位往高位数）作为右子树的返回，那应该还有个1，矛盾
    x = floor((i-1)/2);
    layer = 0;
    while (mod(x, 2) == 1)
        x = floor(x/2);
        layer = layer + 1;
    end
end

function LLR = f(x, y)
    LLR = sign(x) * sign(y) * min(abs(x), abs(y));
end

function LLR = g(x, y, u)
    LLR = (1-2*u)*x + y;
end
