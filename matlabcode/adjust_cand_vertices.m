% 1.排除面积太小的面
% 2.插入‘美化点’

f_choosed  = ones(1, nf)*9; % 1.排除面积太小的面
f_choosed(si_small) = 0; % 不该有 candidate_point 的面
rem_fname = find(f_choosed == 9);

% temp  = 1:20:length(rem_fname); % 2.插入‘美化点’
temp  = 1:3:length(rem_fname); % 2.插入‘美化点’
leap = rem_fname(temp); % leap 选取的面
f_choosed(leap) = 1; % 应该有 candidate_point 点的面

add_idx1 = zeros(1,nf); % 向old_face_name1加入不该插入点的面，索引
add_idx2 = zeros(1,nf); % 向old_face_name2加入插入一个点的面，索引
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

% 查看三类面是否构成了全部面
old_face_name = [old_face_name1;old_face_name2;old_face_name3];
old_face_name = sort(old_face_name);
% if isempty(find(old_face_name - (1:nf)',1))
%     disp('三类面包含了网格所有的面。')
% else
%     warning('三类面分配出错。');
%     return;
% end