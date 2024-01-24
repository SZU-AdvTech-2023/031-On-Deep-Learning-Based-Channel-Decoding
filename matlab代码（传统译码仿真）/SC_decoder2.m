function u = SC_decoder2(LLR, is_info)
% 输入：LLR为一个N×1的向量，存储信道W的LLR，is_info也是N×1的向量，1代表信息位
    N = length(LLR);
    n = log2(N);
    division = 2.^(0:n); % 用于将P和C数组分割的向量
    u = zeros(1, N); % 译码结果，行向量
    P = zeros(N-1, 1); % P用来存储LLR的中间计算结果    
    C = zeros(2*N-1, 2); % C用来存储比特值的中间计算结果，第一列为左子树的beta值，第二列为右子树的beta值
    
    % 流程：更新P，译码ui（ui需添入C(1,1 or 2))，更新C（当i为偶数时进行）
    for i = 1: N
        % 译码ui
        switch i
            case 1
                % 译码u1
                index1 = division(n); % 就是N/2
                for j = 1: index1
                    % 信道LLR,相邻的两个作一次f运算
                    P(index1-1+j) = f(LLR(2*j-1), LLR(2*j));
                end
                for j = n-1:-1:1
                    index1 = division(j);
                    index2 = division(j+1);
                    for k = index1: index2-1
                        P(k) = f(P(2*k), P(2*k+1));
                    end
                end
                % P(1:N-1)更新完毕

            case N/2+1
                % 根节点的左子树全部译码完毕，接下来得译码右子树
                % 和u1基本一致，不同之处在于：
                % u1的P是用f运算，这里是借助左子树得到的C数组用g运算
                index1 = division(n);
                for j = 1: index1
                    % 信道LLR,相邻的两个作一次g运算
                    P(index1-1+j) = g(LLR(2*j-1), LLR(2*j), C(index1-1+j, 1));
                end
                for j = n-1:-1:1
                    index1 = division(j);
                    index2 = division(j+1);
                    for k = index1: index2-1
                        P(k) = f(P(2*k), P(2*k+1));
                    end
                end

            otherwise
                LLR_layer = get_LLR_layer(i);
                % 一直往右上方走，走到的终点的为当前需要进行g运算的右子树的根，这句话好好体会一下
                % layer=0,一次g运算，P(1)；layer=1，2次g运算,P(2:3)；layer=2，4次g运算,P(4:7)
                index1 = division(LLR_layer + 1);
                index2 = division(LLR_layer + 2);
                for j = index1: index2-1
                    P(j) = g(P(2*j), P(2*j+1), C(j,1));
                end
                % ui为某棵子树的最左侧的叶节点，该子树的根为叶节点ui一直向右上方走的终点
                % 该子树树的高度即为LLR_layer，计算ui的LLR还需要进行LLR_layer层f运算，每一层又有若干次f运算
                for level = LLR_layer:-1:1
                    index1 = division(level);
                    index2 = division(level + 1);
                    for j = index1: index2-1
                        P(j) = f(P(2*j), P(2*j+1));
                    end
                end
        end
        is_odd = mod(i, 2); % 1代表i是奇数
        if is_info(i) == 0
            % 冻结比特
            u(i) = 0;
        else
            % 信息位，借助P(1)硬判决ui
            u_LLR = P(1);
            if u_LLR >= 0
                u(i) = 0;
            else
                u(i) = 1;
            end
        end
        % odd=1,C(1,1); odd=0,C(1,2)
        C(1, 2-is_odd) = u(i); % 奇在第一列（左），偶在第二列（右）
        if is_odd == 0
            % （i为偶数时）译码出右节点后，需要进行位计算
            bit_layer = get_bit_layer(i);
            for j = 1:bit_layer
                index1 = division(j);
                index2 = division(j+1);
                for k = index1: index2-1
                    % 重点中的重点！！！！ 
                    % f和g的运算使用的是相邻的两个LLR或者P作为参数，从第一个开始，每两个为一个极化结构的上下LLR值
                    % 每一个2×2的极化结构的上下LLR都是相邻的，上下比特值也是相邻的，下面的代码中C的放置就是相邻的
                    % 这里C放置的位置和网上的代码都不同，网上的代码是没做bitreverse的，我默认都是在G中做了的，所有C放的位置也不同
                    C(2*k, 2) = mod(C(k, 1)+C(k, 2), 2);
                    C(2*k+1, 2) = C(k, 2);
                end
            end
            index1 = division(bit_layer + 1);
            index2 = division(bit_layer + 2);
            for j = index1: index2-1
                % C放的行进行了更改，G进行了bitreverse，每一对极化的LLR和比特值C都是相邻放的！
                C(2*j, 1) = mod(C(j, 1)+C(j, 2), 2);
                C(2*j+1, 1) = C(j, 2);
            end
        end
    end
end

function layer = get_LLR_layer(i)
% 对于ui而言，计算i-1对应的二进制数的低位中连续0的个数
% 即，该树节点能够连续向右上方走的最大步数，也是之后进行f运算的次数
    x = i-1;
    layer = 0;
    while (mod(x, 2) == 0)
        x = x/2;
        layer = layer + 1;
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

function L = f(x, y)
    L = sign(x) * sign(y) * min(abs(x), abs(y));
end

function L = g(x, y, u)
    L = (1-2*u)*x + y;
end