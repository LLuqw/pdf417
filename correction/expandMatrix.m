function [result_matrix]=expandMatrix(res_img,struct_matrix)
[img_row,img_col]=size(res_img);
[struct_row,struct_col]=size(struct_matrix);
%����ṹԪ�ص�����λ��
mid_row=(struct_row+1)/2;
mid_col=(struct_col+1)/2;
%������չ����ĳ��Ϳ�
expand_row=img_row+(struct_row-1);
expand_col=img_col+(struct_col-1);
expand_matrix=zeros(expand_row,expand_col);
%������ʱ����洢����,����ʼ��Ϊ0
tmp_matrix=zeros(expand_row,expand_col);
%��ͼ�����ظ�ֵ����չ����
for row=mid_row:expand_row-mid_row+1
   for col=mid_col:expand_col-mid_col+1
      expand_matrix(row,col)=res_img(row-mid_row+1,col-mid_col+1); 
   end
end
%�ƶ��ṹԪ���е�λ�ã��ṹԪ����ͼ�����ؽ��л����
for row = mid_row : expand_row-mid_row+1
    for col = mid_col : expand_col-mid_col+1
        %ֻ�е��õ�����Ϊ��ɫʱ���Ż���л����
       if(expand_matrix(row,col)==255)
           %�жϽṹԪ�ص��е��Ƿ�ҲΪ��ɫ���ǵĻ���������
           if(struct_matrix(mid_row,mid_col)==1)
              for inner_row=1:struct_row
                 for inner_col=1:struct_col
                     %��ȡ����ʱ��������Ӧ������
                     change_row=row-mid_row+inner_row;
                     change_col=col-mid_col+inner_col;
                     %���ṹ���а�ɫ���ظ�ֵ����ʱ����
                     if(struct_matrix(inner_row,inner_col)==1)
                        tmp_matrix(change_row,change_col)=255;
                     end
                 end
              end
           end
       end
    end
end


%ȥ����չ��ı߽�
result_matrix=zeros(img_row,img_col);
for row=1:img_row
   for col=1:img_col
      result_matrix(row,col)=tmp_matrix(row+mid_row-1,col+mid_col-1); 
   end
end
end