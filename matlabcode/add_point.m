% ��һ��������
if flag_he(i) % �����һ��������Ѿ������˵�
    k1 = flag_he(i);
else          % �����һ������ϻ�û�����µ�
    k1 = idxp;
    flag_he(i) = idxp;
    % ������ҲҪ������Ӧ�ĵ�
    name = hedge_oppo_name(a, b);
    flag_he(name) = idxp;
    idxp = idxp + 1;
end

% �ڶ���������
if flag_he(i + nf) % ����ڶ���������Ѿ������˵�
    k2 = flag_he(i + nf);
else          % ����ڶ�������ϻ�û�����µ�
    k2 = idxp;
    flag_he(i + nf) = idxp;
    % ������ҲҪ������Ӧ�ĵ�
    name = hedge_oppo_name(b, c);
    flag_he(name) = idxp;
    idxp = idxp + 1;
end

% ������������
if flag_he(i + 2*nf) % �����һ��������Ѿ������˵�
    k3 = flag_he(i + 2*nf);
else          % �����һ������ϻ�û�����µ�
    k3 = idxp;
    flag_he(i + 2*nf) = idxp;
    % ������ҲҪ������Ӧ�ĵ�
    name = hedge_oppo_name(c, a);
    flag_he(name) = idxp;
    idxp = idxp + 1;
end