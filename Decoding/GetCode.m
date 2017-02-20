function [ str ] = GetCode(decodes)
%定义相关表格
tc_uc=[65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,32,201,202,204];%大写字母模式
tc_lc=[97,98,99,100,101,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,119,120,121,122,32,205,202,204];%小写字母模式
tc_mi=[48,49,50,51,52,53,54,55,56,57,38,13,09,44,58,35,45,46,36,47,43,37,42,61,94,203,32,201,200,204];%混合型模式
tc_do=[59,60,62,64,91,92,93,95,96,126,33,13,09,44,58,10,45,46,36,47,34,124,42,40,41,63,123,125,39,200];%标点型模式

%待解码码字长度
codelen=decodes(1);  
 % 2->数字型、3->字节型、11->文本模式大写型、12->文本模式小写型、13->文本模式混合型、14->文本模式标点型
mode=11; 
% 用于转移模式时记录模式值,第一个值为1（0）表示当前是（否）转移模式，第二个值表示要返回的模式值
premode=[0,0];
% 用于记录文本模式时的高低位数据，第一个值表示高位数据，第二个值表示低位数据
tcbyte=zeros(1,2);
% 用于记录字节模式和数字模式的缓存序列，其中第一个值表示序列是否有效及何种模式 0->无效、1->数字模式、2->字节模式，第二个值表示序列起始位置
valueindex=[0,0];
str=''; % 存放译码结果
for i=2:codelen
    if decodes(i) >= 900
        if decodes(i) == 900
            mode = 11;   % 11->文本压缩模型
        elseif decodes(i)==902
            mode = 2;    % 2->数字压缩型
        elseif (decodes(i)==901) || (decodes(i)==924) || (decodes(i)==913)
            mode = 3;    % 3->字节压缩型
        end
        
        %数字压缩模式解码
        if valueindex(1) == 1
            caches = decodes(valueindex(2):(i-1)); %数字压缩模式下的待解码码字
            len = size(caches,2);
            for k=1:15:len
                curvalues=zeros(1,15); %用来存储15个码字
                if (k+14)<len %当前组大于15个码字
                    curvalues = caches(k:k+14);
                else          %当前组不足15个码字
                    curvalues(k+15-len:15) = caches(k:len);
                end
                
                %开始解码
                longnum=0;
                for x=1:15
                    %使用符号位来处理
                    longnum = longnum+curvalues(x)*sym(9^(15-x))*(100^(15-x)); 
                end
                tempstr=char(longnum); %转换成字符
                str = strcat(str, tempstr(2:end)); %连接到str中
                
            end
        end 
        valueindex = [0,0];   
    else
        %文本模式处理
        if mode > 10
            % 分离高位字符和低位字符
            tcbyte(1)=floor(decodes(i)/30);  % 高位
            tcbyte(2)=decodes(i) - tcbyte(1)*30; % 低位
            for j=1:2
                if mode==11  % 大写字母
                    if premode(1) == 1
                        mode = premode(2);
                        premode(1) = 0;
                    end
                    if tc_uc(tcbyte(j)+1) == 201 % 锁定为小写字母型子模式
                        mode = 12;
                    elseif tc_uc(tcbyte(j)+1) == 202 % 锁定为混合型子模式
                        mode=13;
                    elseif tc_uc(tcbyte(j)+1) == 204 % 转移到标点模式
                        premode(1)=1;
                        premode(2)=mode;
                        mode=14;
                    else
                        str = strcat(str,char(tc_uc(tcbyte(j)+1))); 
                    end
                elseif mode == 12 % 小写字母
                    if premode(1) == 1
                        mode = premode(2);
                        premode(1)=0;
                    end
                    if tc_lc(tcbyte(j)+1) == 205 % 转移到大写字母型子模式
                        premode(1) = 1;
                        premode(2) = mode;
                        mode = 11;
                    elseif tc_lc(tcbyte(j)+1) == 202 % 锁定为混合型子模式
                        mode = 13;
                    elseif tc_lc(tcbyte(j)+1) == 204  % 转移到标点模式
                        premode(1) = 1;
                        premode(2) = mode;
                        mode = 14;
                    else
                        str = strcat(str,char(tc_lc(tcbyte(j)+1))); % 添加到解码结果字符串中
                    end
                elseif mode == 13 % 混合型子模式
                    if premode(1 )== 1
                        mode = premode(2);
                        premode(1) = 0;
                    end
                    if tc_mi(tcbyte(j)+1) == 200 % 锁定为大写字母模式
                        mode = 11;
                    elseif tc_mi(tcbyte(j)+1) == 201 % 转移到小写字母型子模式
                        mode = 12;
                    elseif tc_mi(tcbyte(j)+1) == 203 % 锁定为标点模式
                        mode = 14;
                    elseif tc_mi(tcbyte(j)+1) == 204 % 转移到标点模式
                        premode(1) = 1;
                        premode(2) = mode;
                        mode = 14;
                    else
                        str = strcat(str,char(tc_mi(tcbyte(j)+1)));  % 添加到解码结果字符串中
                    end
                elseif mode == 14 %标点
                    if premode(1)==1
                        mode = premode(2);
                        premode(1) = 0;
                    end
                    if tc_do(tcbyte(j)+1) == 200  % 锁定为大写字母模式
                        mode=11;
                    else
                        str = strcat(str,char(tc_do(tcbyte(j)+1)));  % 添加到解码结果字符串中
                    end
                end
            end
            
        elseif mode == 2
            if valueindex(1) == 0
                valueindex(1) = 1; % 数字模式
                valueindex(2) = i;
            end
        elseif mode == 3   
            %ignore
        end
    end
end
disp(['解码后的结果：',str])
end

