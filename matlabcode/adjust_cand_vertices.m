% 1.�ų����̫С����
% 2.���롮�����㡯

f_choosed  = ones(1, nf)*9; % 1.�ų����̫С����
f_choosed(si_small) = 0; % ������ candidate_point ����
rem_fname = find(f_choosed == 9);

% temp  = 1:20:length(rem_fname); % 2.���롮�����㡯
temp  = 1:3:length(rem_fname); % 2.���롮�����㡯
leap = rem_fname(temp); % leap ѡȡ����
f_choosed(leap) = 1; % Ӧ���� candidate_point �����

add_idx1 = zeros(1,nf); % ��old_face_name1���벻�ò������棬����
add_idx2 = zeros(1,nf); % ��old_face_name2�������һ������棬����
t1 = 1; t2 = 1;
for i = 1:nf
    if f_choosed(i) == 0
        old_face_name1 = old_face_name1(old_face_name1 ~= i);
        old_face_name2 = old_face_name2(old_face_name2 ~= i);
        old_face_name3 = old_face_name3(old_face_name3 ~= i);
        add_idx1(t1) = i; t1 = t1 + 1;
    end
    
    if (f_choosed(i) == 1) && any(old_face_name1 == i)
        old_face_name1 = old_face_name1(old_face_name1 ~= i);
        add_idx2(t2) = i; t2 = t2 + 1;
    end 
end
add_idx1 = add_idx1(1:t1-1); add_idx2 = add_idx2(1:t2-1);
old_face_name1 = [old_face_name1; add_idx1'];
old_face_name2 = [old_face_name2; add_idx2'];

k1 = length(old_face_name1);
k2 = length(old_face_name2);
k3 = length(old_face_name3);

% �鿴�������Ƿ񹹳���ȫ����
old_face_name = [old_face_name1;old_face_name2;old_face_name3];
old_face_name = sort(old_face_name);
% if isempty(find(old_face_name - (1:nf)',1))
%     disp('������������������е��档')
% else
%     warning('������������');
%     return;
% end