function [result_matrix]=corrosionMatrix(res_img,struct_matrix)
%获取原图像和结构元素的行列
[img_row,img_col]=size(res_img);
[struct_row,struct_col]=size(struct_matrix);
%获取结构元素的中心位置
mid_row=(struct_row+1)/2;
mid_col=(struct_col+1)/2;
%扩展后的行列
corrosion_row=img_row+(struct_row-1);
corrosion_col=img_col+(struct_col-1);
corrosion_matrix=zeros(corrosion_row,corrosion_col);
tmp_matrix=zeros(corrosion_row,corrosion_col);

for row=mid_row:corrosion_row-mid_row+1
   for col=mid_col:corrosion_col-mid_col+1
      corrosion_matrix(row,col)=res_img(row-mid_row+1,col-mid_col+1); 
   end
end

for row=mid_row:corrosion_row-mid_row+1
   for col=mid_col:corrosion_col-mid_col+1
       %如果该点为白色像素才进行处理
      if(corrosion_matrix(row,col)==255) 
          %定义一个标识，判断是否需要进行腐蚀
            flag=1;
            for inner_row=1:struct_row
               for inner_col=1:struct_col
                  change_row=row-mid_row+inner_row;
                  change_col=col-mid_col+inner_col;
                  %检查当结构元素中为白色像素时，原图所对应位置也必须为白色像素，否则设置flag为0，跳出循环
                  if(struct_matrix(inner_row,inner_col)==1)
                     if(corrosion_matrix(change_row,change_col)==0)
                        flag=0;
                        break;
                     end
                  end
               end
               if(flag==0)
                  break; 
               end
            end
            %需要腐蚀，将该点赋值为白色像素
            if(flag==1)
                tmp_matrix(row,col)=255;
            end
      end
   end
end

result_matrix=zeros(img_row,img_col);
for row=1:img_row
   for col=1:img_col
      result_matrix(row,col)=tmp_matrix(row+mid_row-1,col+mid_col-1); 
   end
end

end