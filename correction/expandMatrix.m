function [result_matrix]=expandMatrix(res_img,struct_matrix)
[img_row,img_col]=size(res_img);
[struct_row,struct_col]=size(struct_matrix);
%定义结构元素的中心位置
mid_row=(struct_row+1)/2;
mid_col=(struct_col+1)/2;
%定义扩展矩阵的长和宽
expand_row=img_row+(struct_row-1);
expand_col=img_col+(struct_col-1);
expand_matrix=zeros(expand_row,expand_col);
%定义临时矩阵存储数据,并初始化为0
tmp_matrix=zeros(expand_row,expand_col);
%将图像像素赋值给扩展矩阵
for row=mid_row:expand_row-mid_row+1
   for col=mid_col:expand_col-mid_col+1
      expand_matrix(row,col)=res_img(row-mid_row+1,col-mid_col+1); 
   end
end
%移动结构元素中点位置，结构元素与图像像素进行或操作
for row = mid_row : expand_row-mid_row+1
    for col = mid_col : expand_col-mid_col+1
        %只有当该点像素为白色时，才会进行或操作
       if(expand_matrix(row,col)==255)
           %判断结构元素的中点是否也为白色，是的话进行膨胀
           if(struct_matrix(mid_row,mid_col)==1)
              for inner_row=1:struct_row
                 for inner_col=1:struct_col
                     %获取到临时矩阵所对应的行列
                     change_row=row-mid_row+inner_row;
                     change_col=col-mid_col+inner_col;
                     %将结构体中白色像素赋值给临时矩阵
                     if(struct_matrix(inner_row,inner_col)==1)
                        tmp_matrix(change_row,change_col)=255;
                     end
                 end
              end
           end
       end
    end
end


%去除扩展后的边界
result_matrix=zeros(img_row,img_col);
for row=1:img_row
   for col=1:img_col
      result_matrix(row,col)=tmp_matrix(row+mid_row-1,col+mid_col-1); 
   end
end
end