 clear all;
close all;

I=imread('lv3 - hard.jpg');
I=rgb2gray(I);
[Row, Col] = size(I);%图像尺寸
figure;
imshow(I);
title('原始图像');
%OSTU二值化算法
[Count,x] = imhist(I);%统计各灰度级像素个数count,对应像素级x;

Count = Count./(Row*Col);%每个灰度级在图像中占的比例
L = 256; %0-255个灰度值
for i=1:L
    if(Count(i) ~= 0) %找到左边非0的灰度值
        st=i;
        break;
    end
end
for i=L:-1:1
    if(Count(i) ~= 0) %找到右边非0的灰度值
        nd=i;
        break;
    end
end
Count=Count(st:nd);%去掉两边不存在的灰度级
x=x(st:nd);
wo=zeros(L,1);
for t=1:nd+1-st
    wo(t)=sum(Count(1:t));%计算前t个像素的累计概率
    uo(t)=sum(x(1:t).*Count(1:t))/wo(t);%前t个像素的像素期望
end
for t=1:nd-st
    g(t)=wo(t)*(uo(nd+1-st)-uo(t))^2./(1-wo(t));
end%找到最佳阈值
[MAX,pos]=max(g);
MAX=x(pos);

res_img=zeros(Row,Col);
for i=1:Row
    for j=1:Col
        if(I(i,j)>=MAX) res_img(i,j)=0;
        else res_img(i,j)=255;
        end
    end
end%根据阈值二值化
figure;
imshow(res_img);
title('二值化');



% %线扫描
% leftline=zeros(Row,2);
% k=1;
% for i=1:Row
%     for j=1:Col
%         if(res_img(i,j)==0)
%             if(k==1 || (j-leftline(k-1,2)<10))
%                 leftline(k,1)=i;
%                 leftline(k,2)=j;
%                 k=k+1;
%                 break;
%             end
%         end
%     end
% end
% k;
% figure;
% plot(leftline(1:k,2));
% title('线扫描');

%膨胀函数 imdilate
%腐蚀函数 imerode
res_img_tmp = res_img;%开闭运算后的图
for i=1:1
    %开运算
    SE_small = strel('disk',1);%较小结构元素1
    res_img_tmp = corrosionMatrix(res_img_tmp,SE_small.getnhood);
    res_img_tmp = expandMatrix(res_img_tmp,SE_small.getnhood);
%      res_img_tmp = imerode(res_img_tmp,SE_small);
%      res_img_tmp = imdilate(res_img_tmp,SE_small);
    %闭运算
    SE_big = strel('disk',12);%较大结构元素24
%     res_img_tmp = imdilate(res_img_tmp,SE_big);
%     res_img_tmp = imerode(res_img_tmp,SE_big);
    res_img_tmp = expandMatrix(res_img_tmp,SE_big.getnhood);
    res_img_tmp = corrosionMatrix(res_img_tmp,SE_big.getnhood);
    
end
figure;
imshow(res_img_tmp);
title('开闭运算后');

%计算连通域数量
[L,num]=bwlabel(res_img_tmp);
%统计连通域面积
s=zeros(num,1);
for i=1:Row
    for j=1:Col
        if(L(i,j)~=0)
            s(L(i,j))=s(L(i,j))+1;
        end
    end
end
%找出最大得连通域
max_area=-1;
for i=1:num
    if(s(i)>max_area)
        max_area = s(i);
    end
end
%去掉面积小于最大连通域一半的连通域
for i=1:Row
    for j=1:Col
        if(L(i,j)~=0 && s(L(i,j)) < max_area/2)
            res_img_tmp(i,j)=0;
        end
        if(res_img_tmp(i,j)==0)
            res_img(i,j)=0;
        end
    end
end
figure;
imshow(res_img);
title('去掉小面积连通域');
res_img_kuang = zeros(Row,Col);%二维码的边框
SE_small = strel('disk',1);
res_img_tmp2 = expandMatrix(res_img_tmp,SE_small.getnhood);%较小的结构元素膨胀后的图
for i=1:Row
    for j=1:Col
        res_img_kuang(i,j)=res_img_tmp2(i,j)-res_img_tmp(i,j);
        %去除噪声
        if(res_img_tmp(i,j)==0) 
            res_img(i,j)=0;
        end
    end
end
%去除噪声
struct_ele=strel('disk',1);
res_img = expandMatrix(res_img,struct_ele.getnhood);
res_img = corrosionMatrix(res_img,struct_ele.getnhood);
figure;
imshow(res_img_kuang);
title('边框');

%极坐标系--------------------------------------------------------
H=zeros(2*round((Row^2+Col^2)^0.5)+1,180);
for i=1:Row
    for j=1:Col
        if(res_img_kuang(i,j) == 255)
            for n=1:180
                c = round(j*cos(n/180*pi)+i*sin(n/180*pi))+round((Row^2+Col^2)^0.5);
                H(c,n) = H(c,n) + 1;
            end
        end
    end
end
lines2 = zeros(4,2);
k=1;
H_pre=H(:,1:90);
H_nex=H(:,91:180);
for i=1:2
    [lines2(k,1),lines2(k,2)]=matr_max(H_pre);
%     [lines2(k,1),lines2(k,2)]=find(H_pre == max(max(H_pre)));
    for j=lines2(k,1)-10:lines2(k,1)+10
        for p=lines2(k,2)-7:lines2(k,2)+7
            H_pre(j,p)=0;
        end
    end
    k=k+1;
    [lines2(k,1),lines2(k,2)]=matr_max(H_nex);
%     [lines2(k,1),lines2(k,2)]=find(H_nex == max(max(H_nex)));
    for j=lines2(k,1)-16:lines2(k,1)+16
        for p=lines2(k,2)-5:lines2(k,2)+5
            H_nex(j,p)=0;
        end
    end
    k=k+1;
end
lines2(2,2)=lines2(2,2)+90;
lines2(4,2)=lines2(4,2)+90;
lines2
%找出四条线的斜率跟截距
k1=-cot(lines2(1,2)/180*pi);b1=(lines2(1,1)-round((Row^2+Col^2)^0.5))/sin(lines2(1,2)/180*pi);
k2=-cot(lines2(2,2)/180*pi);b2=(lines2(2,1)-round((Row^2+Col^2)^0.5))/sin(lines2(2,2)/180*pi);
k3=-cot(lines2(3,2)/180*pi);b3=(lines2(3,1)-round((Row^2+Col^2)^0.5))/sin(lines2(3,2)/180*pi);
k4=-cot(lines2(4,2)/180*pi);b4=(lines2(4,1)-round((Row^2+Col^2)^0.5))/sin(lines2(4,2)/180*pi);
%四个交点
[x1(1),y1(1)]=find_node(k1,b1,k2,b2);
[x1(2),y1(2)]=find_node(k1,b1,k4,b4);
[x1(3),y1(3)]=find_node(k3,b3,k2,b2);
[x1(4),y1(4)]=find_node(k3,b3,k4,b4);

%对点进行排序
for i=1:3
    for j=i+1:4
        if(x1(i)^2+y1(i)^2 > x1(j)^2+y1(j)^2)
            tmp_x=x1(i);
            tmp_y=y1(i);
            x1(i)=x1(j);
            y1(i)=y1(j);
            x1(j)=tmp_x;
            y1(j)=tmp_y;
        end
    end
end

%二维码两条边
L12 = ((x1(1)-x1(2))^2+(y1(1)-y1(2))^2)^0.5;
L13 = ((x1(1)-x1(3))^2+(y1(1)-y1(3))^2)^0.5;
if(L12 < L13)
    xielv=(y1(3)-y1(1))/(x1(3)-x1(1));
else
    xielv=(y1(2)-y1(1))/(x1(2)-x1(1));
end
line([0,Col],[(lines2(1,1)-round((Row^2+Col^2)^0.5)-0*cos(lines2(1,2)/180*pi))/sin(lines2(1,2)/180*pi),(lines2(1,1)-round((Row^2+Col^2)^0.5)-Col*cos(lines2(1,2)/180*pi))/sin(lines2(1,2)/180*pi)]);
line([0,Col],[(lines2(2,1)-round((Row^2+Col^2)^0.5)-0*cos(lines2(2,2)/180*pi))/sin(lines2(2,2)/180*pi),(lines2(2,1)-round((Row^2+Col^2)^0.5)-Col*cos(lines2(2,2)/180*pi))/sin(lines2(2,2)/180*pi)]);
line([(lines2(3,1)-round((Row^2+Col^2)^0.5)-0*sin(lines2(3,2)/180*pi))/cos(lines2(3,2)/180*pi),(lines2(3,1)-round((Row^2+Col^2)^0.5)-Row*sin(lines2(3,2)/180*pi))/cos(lines2(3,2)/180*pi)],[0,Row]);
line([(lines2(4,1)-round((Row^2+Col^2)^0.5)-0*sin(lines2(4,2)/180*pi))/cos(lines2(4,2)/180*pi),(lines2(4,1)-round((Row^2+Col^2)^0.5)-Row*sin(lines2(4,2)/180*pi))/cos(lines2(4,2)/180*pi)],[0,Row]);
hold on;
plot(x1(1),y1(1),'Marker','o','Color','red');
plot(x1(2),y1(2),'Marker','o','Color','blue');
plot(x1(3),y1(3),'Marker','o','Color','green');
plot(x1(4),y1(4),'Marker','o','Color','yellow');
hold off;

% %旋转变换
% % n=round(k/2);
% % for i=1:n
% %     tanarray(i)=(leftline(i+n,2)-leftline(i,2))/(leftline(i+n,1)-leftline(i,1));
% % end %算出所有斜率
% % tanarray=sort(tanarray);
% % ang=atan(0-tanarray(round(n/2)))%找出中间的斜率，并换算成角度(弧度制)
% % tmp=imrotate(res_img,ang*180/pi,'bilinear');
% % figure;
% % imshow(tmp)
% ang=atan(xielv);
% new_row=round(Col*sin(ang)+Row*cos(ang));
% new_col=round(Col*cos(ang)+Row*sin(ang));
% dst=zeros(new_row,new_col);
% for i=1:Row
%     for j=1:Col
%         X=(i-Row/2)*cos(ang)-(j-Col/2)*sin(ang)+(new_row/2);
%         Y=(i-Row/2)*sin(ang)+(j-Col/2)*cos(ang)+(new_col/2);
%         X=abs(round(X))+1;
%         Y=abs(round(Y))+1;
%         if(j==x1(1)&&i==y1(1)) loca(1,1)=X;loca(1,2)=Y;
%         elseif(j==x1(4)&&i==y1(4)) loca(2,1)=X;loca(2,2)=Y;
%         end
%         dst(X,Y)=res_img(i,j);
%     end
% end%图象旋转，以图像中心为旋转点
% figure;
% imshow(dst);
% title('旋转变换后');
% hold on;
% plot(loca(1,2),loca(1,1),'Marker','o','Color','red');
% plot(loca(2,2),loca(2,1),'Marker','o','Color','red');
% hold off;
% 
% %仿射变换---------------------------- 用于lv0，1，2
% %原三组点
% in_points = [x1(1),y1(1) 1;x1(2) y1(2) 1;x1(3) y1(3) 1]
% flag=0;%1则结果图像需要逆时针转90度
% %对应三组点
% if(abs(y1(1)-y1(2)) < abs(x1(1) -x1(2)))
%     out_points=[50,50,1;50+L12,50,1;50,50+L13,1];
%     if(L13 > L12)
%         flag=1;
%     end
% elseif(abs(y1(1)-y1(2)) > abs(x1(1) -x1(2)))
%     out_points=[50,50,1;50,50+L12,1;50+L13,50,1];
%     if(L12 > L13)
%         flag=1;
%     end
% end
% %变换矩阵 /右除 AB=C A=C/B;\左除 AB=C  B=A\C
% tform=in_points\out_points
% kk=1;
% for i=1:Row
%     for j=1:Col
%         if(res_img(i,j) == 255)
%             mmt(kk,:)=[j,i,1];
%             kk=kk+1;
%         end
%     end
% end
% kk
% mmt=mmt*tform;
% %结果图像
% for i=1:kk-1
%     dst(round(mmt(i,2)+50),round(mmt(i,1))+50)=255;
% end

%双线性变换----------------------------------------------------用于lv3
in_points=[x1(1) x1(2) x1(3) x1(4); y1(1) y1(2) y1(3) y1(4);x1(1)*y1(1) x1(2)*y1(2) x1(3)*y1(3) x1(4)*y1(4);1 1 1 1];
flag=0;
if(abs(y1(1)-y1(2)) < abs(x1(2) -x1(1)))
    out_points=[50 50 50+L13 50+L13;50 50+L12 50 50+L12];
    if(L13 > L12)
        flag=1;
    end
elseif(abs(y1(1)-y1(2)) > abs(x1(2) -x1(1)))
    out_points=[50 50+L12 50 50+L12;50 50+L13 50 50+L13];
    if(L12 > L13)
        flag=1;
    end
end
%变换矩阵
doublelinear=out_points/in_points
kk=1;
for i=1:Row
    for j=1:Col
        if(res_img(i,j) == 255)
            mt(1,kk)=j;mt(2,kk)=i;mt(3,kk)=i*j;mt(4,kk)=1;
            kk=kk+1;
        end
    end
end
mt=doublelinear*mt;
for i=1:kk-1
    dst(round(mt(1,i))+50,round(mt(2,i))+50)=255;
end
[new_row,new_col]=size(dst);
%转为黑色二维码，白色背景
for i=1:new_row
    for j=1:new_col
        if(dst(i,j)==255) dst(i,j)=0;
        else dst(i,j)=255;
        end
    end
end


figure;
imshow(dst);
title('插值前');
[new_row,new_col]=size(dst);
%邻近插值法
for i=2:new_row-1
    for j=2:new_col-1
        if(dst(i,j)~=dst(i-1,j)&& dst(i,j)~=dst(i+1,j)&& dst(i,j)~=dst(i,j-1)&& dst(i,j)~=dst(i,j+1))
            dst(i,j)=dst(i,j-1);
        end
    end
end
figure;
imshow(dst);
title('邻近插值');

% if(loca(2,1)>loca(1,1)&&loca(2,2)>loca(1,2)) mmin=loca(2,1)-loca(1,1);mmax=loca(2,2)-loca(1,2);
% else mmin=loca(1,1)-loca(2,1);mmax=loca(2,2)-loca(1,2);
% end
% for i=1:mmin;
%     for j=1:mmax;
%         if(loca(2,1)>loca(1,1)&&loca(2,2)>loca(1,2)) res_dst(i,j)=dst(i+loca(1,1),j+loca(1,2));
%         else res_dst(i,j)=dst(loca(2,1)+i,j+loca(1,2));
%         end
%     end
% end
% figure;
% imshow(res_dst);
% title('裁剪后');
if(flag==1||(flag==0&&L13 < L12))
   mmin=L13;mmax=L12;
else
    mmax=L13;mmin=L12;
end
for i=1:mmin
    for j=1:mmax
        res_dst(i,j)=dst(49+i+50,49+j+50);
    end
end
if(flag == 1)
    res_dst=imrotate(res_dst,90,'bilinear');
end
figure;
imshow(res_dst);
title('裁剪后');