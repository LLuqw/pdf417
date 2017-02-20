clear all;
close all;

I=imread('test3.png');
%I=rgb2gray(I);%test1、2、3不用加上这句，test4、5需要
[Row,Col] = size(I);
res_img = mat2gray(I);%归一化灰度值
%计算层数--------------------------------------------------------------------------------------------
%Robert算子提取边缘.保存到new_res_img
new_res_img=res_img;%为保留图像的边缘一个像素
robertsNum=0; %经roberts算子计算得到的每个像素的值
robertThreshold=0.2; %设定阈值
for j=1:Row-1 %进行边界提取
    for k=2:Col-1
        robertsNum = abs(res_img(j,k)-res_img(j+1,k+1)) + abs(res_img(j+1,k)-res_img(j,k+1));
        if(robertsNum > robertThreshold)
            new_res_img(j,k)=255;
        else
            new_res_img(j,k)=0;
        end
    end
end
figure,imshow(new_res_img);
title('roberts算子的处理结果')

% 水平投影，去掉竖直方向的线
hori = zeros(Row,Col);
for i=1:Row
    for j=2:Col-1
        if(new_res_img(i,j) == 255 && (new_res_img(i,j-1) == 255 && new_res_img(i,j+1) == 255))%如果左右都是白点
            hori(i,j) = 1;
        end
    end
end
figure;
imshow(hori);title('水平投影');

%统计每一行的白点，并找出峰值,保存在Line数组
for i=1:Row
    k(i)=0;
    for j=1:Col
        if(hori(i,j) == 1)
            k(i)=k(i)+1;
        end
    end
end
figure;
plot(k);title('峰值');
num_top=0;%峰值个数
for i=2:Row-1
    if(k(i-1) < k(i)&&k(i) > k(i+1))
        num_top=num_top+1;
        Line(num_top)=i;
    end
end
num_top=num_top+1;
Line(num_top)=Row;%加上下面最后一行
%layer数组保存每一层的中心行数
for i=1:num_top
    if(i==1)
        layer(i)=round(Line(i)/2);
    else
        layer(i)=round((Line(i-1)+Line(i))/2);
    end
end
figure;
imshow(I);
for i=1:num_top
    %line([0,Col],[Line(i),Line(i)]);%画出峰值线，也就是层与层的分隔
    line([0,Col],[layer(i),layer(i)]);%画出层的中心线
end
%提取码字
% 垂直投影，去掉水平方向的线
veri = zeros(Row,Col);
for i=2:Row-1
    for j=1:Col
        if(new_res_img(i,j) == 255 && (new_res_img(i-1,j) == 255 && new_res_img(i+1,j) == 255))%如果左右都是白点
            veri(i,j) = 1;
        end
    end
end
figure;
imshow(veri);title('垂直投影');
%计算距离，保存在code
code=zeros(num_top,Col);
for i=1:num_top
    mark=1;
    for j=2:Col
        if(veri(layer(i),j) == 1)
            code(i,j)=j-mark;
            mark=j;
        end
    end
end
%去掉code中为0的数据，结果保存在ccode
for i=1:num_top
    pos=1;
    for j=2:Col
        if(code(i,j) ~= 0)
            ccode(i,pos)=code(i,j);
            pos=pos+1;
        end
    end
end
aunit=round(sum(ccode(:,1))/num_top/8);
ccode=round(ccode./aunit);
[ccode_row,ccode_col]=size(ccode);
ccode=ccode(1:num_top,17:ccode_col-16);%去掉起始跟终止符
%读取PDF417符号转化表
load 'symcodes.mat' -ascii;
[symcodes_row,symcodes_col]=size(symcodes);
%查找码字，用temp保存查表数据，decode保存查表后的结果
for i=1:num_top
    for j=0:(ccode_col-16-17)/8
        %计算转换后的值
        temp(i,j+1)=0;
        for l=1:8
            temp(i,j+1)=temp(i,j+1)+ccode(i,j*8+l)*10^(8-l);
        end
        %转换成码字
        for p=1:symcodes_col
            if(temp(i,j+1)==symcodes(mod(i-1,3)+1,p))
                decode(i,j+1)=p-1;
                break;
            end
        end
    end
end
%将二维数组转换为一维数组来处理
[m , n ]= size(decode);
decodes = zeros(m * n);
k = 1;
for i = 1: m
   for j =1 : n
      decodes(k) = decode(i ,j);
      k= k+1;
   end
end
str=GetCode(decodes);



