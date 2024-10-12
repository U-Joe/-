% 读取图像
img_y = imread('ty.bmp');
img = rgb2gray(img_y);
bw_b = imbinarize(img);  % 二值化,将背景和前景分离
bw = imcomplement(bw_b);
bw_fill = imfill(bw, 'holes');  % 填充孔洞
% 计算面积和周长
edge_img = edge(bw, 'canny');
count = sum (sum (edge_img==1));
stats = regionprops( edge_img, 'all');
area = stats.FilledArea;
perimeter = stats.PerimeterOld;

% 计算致密性(区域致密性通常指的是形状的紧密程度，可以通过面积与周长的平方比值来计算)
compactness = perimeter^2 / area;

%归一化致密度
gui_compactness = 1-(4*pi)/compactness;

%计算图像的二阶矩和三阶矩
[row,col]=size(bw);
m00=0;m10=0;m01=0;
m20=0; m02=0;m11=0;
m30=0;m12=0;m21=0;m03=0;
for i=1:row
    for j=1:col
        m00=m00+bw(i,j);
        m10=m10+i*bw(i,j);
        m01=m01+j*bw(i,j);
        m20=m20+i^2*bw(i,j);
        m02=m02+j^2*bw(i,j);
        m11=m11+j*i*bw(i,j);
        m30=m30+i^3*bw(i,j);
        m12=m12+i*j^2*bw(i,j);
        m21=m21+i^2*j*bw(i,j);
        m03=m03+j^3*bw(i,j);
    end
end

i_j =m10/m00;
j_j = m01/m00;
u01=0;u10=0;u00=m00;
u11=m11-((m10*m01)/m00);u20=m20-i_j*m10;
u02=m02-j_j*m01;u12=m12-2*j_j*m11-i_j*m02+2*j_j^2*m10;
u21=m21-2*i_j*m11-j_j*m20+2*i_j^2*m01;
u30=m30-3*i_j*m20+2*i_j^2*m10;
u03=m03-3*j_j*m02+2*j_j^2*m01;

%形状的主方向,x 轴与椭圆长轴（该椭圆与区域具有相同的二阶矩）之间的角度

orientation = 0.5*atan((2*u11)/(u20-u02))

%偏心度是指形状的离心率，是椭圆焦距与其长轴的比值：
eccentricity = ((u20 - u02)^2 + 4 * (u11^2))/(u20+u02)^2;  % 偏心度

stats = regionprops( edge_img, 'all');
centroid = stats.Centroid;
[rows, cols] = find(edge_img);
for i=1:count 
    distances(i) = sqrt((rows(i) - centroid(2)).^2 + (cols(1) - centroid(1)).^2); % 从形心到边界各点的距离
end
for i=1:count 
    g_distances(i) = distances(i)/ max(distances);
end

mmm = 0;
for i=1:count 
    mmm = mmm+g_distances(i);
end
mmm = (1/count)*mmm;
uuu2 = 0;
for i=1:count 
    uuu2 = uuu2+(g_distances(i)-mmm)^2;
end
uuu2 = (1/count)*uuu2;
uuu4 = 0;
for i=1:count 
    uuu4 = uuu4+(g_distances(i)-mmm)^4;
end
uuu4 = (1/count)*uuu4;
f1=uuu2^(1/2)/mmm;
f2=uuu4^(1/4)/mmm;
f21=f2-f1;

% 计算径向距离测度（例如，平均径向距离）
mean_radial_distance = mean(distances);

%不变空间矩是图像的一些不变量，在平移、旋转和尺度变化时保持不变。可以通过Hu矩来计算：
% 归一化中心矩
n20 = u20 / u00^2;
n02 = u02 / u00^2;
n11 = u11 / u00^2;
n30 = u30 / u00^2.5;
n12 = u12 / u00^2.5;
n21 = u21 / u00^2.5;
n03 = u03 / u00^2.5;

% 计算Hu不变矩
hu_moments = calculateHuMoments(n20, n02, n11, n30, n03, n12, n21);

 % 输出结果
fprintf('面积: %f\n', area);
fprintf('周长: %f\n', perimeter);
fprintf('区域致密度: %f\n', compactness);
fprintf('归一化致密度: %f\n', gui_compactness);
fprintf('主方向: %f\n', orientation);
fprintf('偏心度: %f\n', eccentricity);
fprintf('径向距离方差: %f\n', f21);
fprintf('Hu不变矩: %s\n\n', mat2str(hu_moments));

% 计算Hu不变矩的函数
function hu = calculateHuMoments(n20, n02, n11, n30, n03, n12, n21)
    % 计算七个Hu不变矩
    hu(1) = n20 + n02;
    hu(2) = (n20 - n02)^2 + 4 * n11^2;
    hu(3) = (n30 - 3*n12)^2 + (3*n21 - n03)^2;
    hu(4) = (n30 + n12)^2 + (n21 + n03)^2;
    hu(5) = (n30 - 3*n12) * (n30 + n12) * ((n30 + n12)^2 - 3 * (n21 + n03)^2) + ...
            (3*n21 - n03) * (n21 + n03) * (3 * (n30 + n12)^2 - (n21 + n03)^2);
    hu(6) = (n20 - n02) * ((n30 + n12)^2 - (n21 + n03)^2) + 4 * n11 * (n30 + n12) * (n21 + n03);
    hu(7) = (3 * n21 - n03) * (n30 + n12) * ((n30 + n12)^2 - 3 * (n21 + n03)^2) - ...
            (n30 - 3 * n12) * (n21 + n03) * (3 * (n30 + n12)^2 - (n21 + n03)^2);
end