function h = hmat(svmat,usrpos)
%输入参数：
% svmat：在用户定义的直角坐标系卫星坐标矩阵
% usrpos：在用户定义的直角坐标系中用户的位置
% 输出参数：
% h：GPS定位的方向余弦矩阵

N = max(size(svmat));
if N < 4,
    error('insufficient number of satellites')
else
    tmppos = usrpos;
    [m,n] = size(tmppos);
    if m > n, tmppos = tmppos'; end,
    h = ones(N,4);
    % linear expansion in user position, the precision depends on the user coordinates
    for i = 1:N,
        tmpvec = tmppos - svmat(i,:);
        h(i,1:3) = tmpvec./norm(tmpvec);
    end,
end,
