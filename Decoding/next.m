clear all;
close all;

I=imread('test3.png');
%I=rgb2gray(I);%test1��2��3���ü�����䣬test4��5��Ҫ
[Row,Col] = size(I);
res_img = mat2gray(I);%��һ���Ҷ�ֵ
%�������--------------------------------------------------------------------------------------------
%Robert������ȡ��Ե.���浽new_res_img
new_res_img=res_img;%Ϊ����ͼ��ı�Եһ������
robertsNum=0; %��roberts���Ӽ���õ���ÿ�����ص�ֵ
robertThreshold=0.2; %�趨��ֵ
for j=1:Row-1 %���б߽���ȡ
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
title('roberts���ӵĴ�����')

% ˮƽͶӰ��ȥ����ֱ�������
hori = zeros(Row,Col);
for i=1:Row
    for j=2:Col-1
        if(new_res_img(i,j) == 255 && (new_res_img(i,j-1) == 255 && new_res_img(i,j+1) == 255))%������Ҷ��ǰ׵�
            hori(i,j) = 1;
        end
    end
end
figure;
imshow(hori);title('ˮƽͶӰ');

%ͳ��ÿһ�еİ׵㣬���ҳ���ֵ,������Line����
for i=1:Row
    k(i)=0;
    for j=1:Col
        if(hori(i,j) == 1)
            k(i)=k(i)+1;
        end
    end
end
figure;
plot(k);title('��ֵ');
num_top=0;%��ֵ����
for i=2:Row-1
    if(k(i-1) < k(i)&&k(i) > k(i+1))
        num_top=num_top+1;
        Line(num_top)=i;
    end
end
num_top=num_top+1;
Line(num_top)=Row;%�����������һ��
%layer���鱣��ÿһ�����������
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
    %line([0,Col],[Line(i),Line(i)]);%������ֵ�ߣ�Ҳ���ǲ����ķָ�
    line([0,Col],[layer(i),layer(i)]);%�������������
end
%��ȡ����
% ��ֱͶӰ��ȥ��ˮƽ�������
veri = zeros(Row,Col);
for i=2:Row-1
    for j=1:Col
        if(new_res_img(i,j) == 255 && (new_res_img(i-1,j) == 255 && new_res_img(i+1,j) == 255))%������Ҷ��ǰ׵�
            veri(i,j) = 1;
        end
    end
end
figure;
imshow(veri);title('��ֱͶӰ');
%������룬������code
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
%ȥ��code��Ϊ0�����ݣ����������ccode
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
ccode=ccode(1:num_top,17:ccode_col-16);%ȥ����ʼ����ֹ��
%��ȡPDF417����ת����
load 'symcodes.mat' -ascii;
[symcodes_row,symcodes_col]=size(symcodes);
%�������֣���temp���������ݣ�decode�������Ľ��
for i=1:num_top
    for j=0:(ccode_col-16-17)/8
        %����ת�����ֵ
        temp(i,j+1)=0;
        for l=1:8
            temp(i,j+1)=temp(i,j+1)+ccode(i,j*8+l)*10^(8-l);
        end
        %ת��������
        for p=1:symcodes_col
            if(temp(i,j+1)==symcodes(mod(i-1,3)+1,p))
                decode(i,j+1)=p-1;
                break;
            end
        end
    end
end
%����ά����ת��Ϊһά����������
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



