function [result_matrix]=corrosionMatrix(res_img,struct_matrix)
%��ȡԭͼ��ͽṹԪ�ص�����
[img_row,img_col]=size(res_img);
[struct_row,struct_col]=size(struct_matrix);
%��ȡ�ṹԪ�ص�����λ��
mid_row=(struct_row+1)/2;
mid_col=(struct_col+1)/2;
%��չ�������
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
       %����õ�Ϊ��ɫ���زŽ��д���
      if(corrosion_matrix(row,col)==255) 
          %����һ����ʶ���ж��Ƿ���Ҫ���и�ʴ
            flag=1;
            for inner_row=1:struct_row
               for inner_col=1:struct_col
                  change_row=row-mid_row+inner_row;
                  change_col=col-mid_col+inner_col;
                  %��鵱�ṹԪ����Ϊ��ɫ����ʱ��ԭͼ����Ӧλ��Ҳ����Ϊ��ɫ���أ���������flagΪ0������ѭ��
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
            %��Ҫ��ʴ�����õ㸳ֵΪ��ɫ����
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