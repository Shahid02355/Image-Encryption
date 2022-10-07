
clc;  close all;   clear all;

theta = 0;
z1 = 2; z2 = 3; z3 = 4;
x(1)  = cos(theta*pi*(z1^1))^2;
y(1) = cos(theta*pi*(z2^1))^2;
z(1) = cos(theta*pi*(z3^1))^2;

for ii = 2:1:7000000
   x(ii) = cos(z1*acos(sqrt(x(ii-1))))^2;
   y(ii) = cos(z2*acos(sqrt(y(ii-1))))^2;
   z(ii) = cos(z3*acos(sqrt(z(ii-1))))^2;
end
x=ceil(mod((x*1000000000000000),512));
 y=ceil(mod((y*100000000000000),512));
 z=ceil(mod((z*100000000000000),512));
image=imread('download.jpg');
R=image(:,:,1);
G=image(:,:,2);
B=image(:,:,3);
rgbR=R;
[rowR,colR]=size(rgbR);
n=800;
p=900;
q=1000;
for j=1:1:rowR
    k(j)=x(j+n);
    l(j)=y(j+p);
end

% for i=1:1:row
%     k(i)=x(i+n);
% end
% for j=1:1:col
%     l(j)=y(j+p);
% end

for j=1:1:colR*rowR
    m1R(j)=z(j+q);
end
mR=uint8(255 * mat2gray(m1R));


for i=1:1:rowR
    for j=1:1:colR
        if(mod(k(i),2)==0)
            if((j+k(i))<=colR)   %shift right of row
               sh_rowG(i,j+k(i))=rgbR(i,j);
               row_shift_evenG(i,j)=j+k(i);
            else
               sh_rowG(i,(j+k(i)-colR))=rgbR(i,j); 
               row_shift_evenG(i,j)=(j+k(i)-colR);
            end
        else
            if((j-k(i))>=1)       %shift left of row
                sh_rowG(i,j-k(i))=rgbR(i,j);
                row_shift_oddB(i,j)=j-k(i);
            else
                sh_rowG(i,(colR+(j-k(i))))=rgbR(i,j);
                row_shift_oddB(i,j)=colR+j-k(i);
            end
        end
    end
end
   

for j=1:1:colR
    for i=1:1:rowR
        if(mod(l(j),2)==0)
           if((i-l(j))>=1)          %shift up of column
                sh_colG(i-l(j),j)=sh_rowG(i,j);
                col_shift_evenG(i,j)=i-l(j);
            else
                sh_colG((rowR+i-l(j)),j)=sh_rowG(i,j);
                col_shift_evenG(i,j)=rowR+i-l(j);
            end
        else
           if((i+l(j))<=rowR)           %shift down of column
               sh_colG(i+l(j),j)=sh_rowG(i,j);
               col_shift_oddB(i,j)=i+l(j);
               else
               sh_colG((i+l(j)-rowR),j)=sh_rowG(i,j); 
               col_shift_oddB(i,j)=(i+l(j)-rowR);
            end
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%XOR IMAGE%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
total_lengthG=rowR*colR;
column_imageG=reshape(sh_colG,1,total_lengthG);
for i=1:1:total_lengthG
xorr1G(1,i)=bitxor(column_imageG(i),mR(i));
end
enc_R=reshape(xorr1G,rowR,colR);
%Green Encryption
rgbG=G;
[rowG,colG]=size(rgbG);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%INITIALIZE THE VALUE OF ROTATION%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for j=1:1:rowG
    k(j)=x(j+n);
    l(j)=y(j+p);
end

% for i=1:1:row
%     k(i)=x(i+n);
% end
% for j=1:1:col
%     l(j)=y(j+p);
% end

for j=1:1:colG*rowG
    m1G(j)=z(j+100);
end
mG=uint8(255 * mat2gray(m1G));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% IF k is even right shift row else left shift row%%%%%%%%%%%%%%%%%
% If l is even shift up column else down shift column%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%ROTATION OPERATION%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i=1:1:rowG
    for j=1:1:colG
        if(mod(k(i),2)==0)
            if((j+k(i))<=colG)   %shift right of row
               sh_rowG(i,j+k(i))=rgbG(i,j);
               row_shift_evenG(i,j)=j+k(i);
            else
               sh_rowG(i,(j+k(i)-colG))=rgbG(i,j); 
               row_shift_evenG(i,j)=(j+k(i)-colG);
            end
        else
            if((j-k(i))>=1)       %shift left of row
                sh_rowG(i,j-k(i))=rgbG(i,j);
                row_shift_oddB(i,j)=j-k(i);
            else
                sh_rowG(i,(colG+(j-k(i))))=rgbG(i,j);
                row_shift_oddB(i,j)=colG+j-k(i);
            end
        end
    end
end
   

for j=1:1:colG
    for i=1:1:rowG
        if(mod(l(j),2)==0)
           if((i-l(j))>=1)          %shift up of column
                sh_colG(i-l(j),j)=sh_rowG(i,j);
                col_shift_evenG(i,j)=i-l(j);
            else
                sh_colG((rowG+i-l(j)),j)=sh_rowG(i,j);
                col_shift_evenG(i,j)=rowG+i-l(j);
            end
        else
           if((i+l(j))<=rowG)           %shift down of column
               sh_colG(i+l(j),j)=sh_rowG(i,j);
               col_shift_oddB(i,j)=i+l(j);
               else
               sh_colG((i+l(j)-rowG),j)=sh_rowG(i,j); 
               col_shift_oddB(i,j)=(i+l(j)-rowG);
            end
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%XOR IMAGE%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
total_lengthG=rowG*colG;
column_imageG=reshape(sh_colG,1,total_lengthG);
for i=1:1:total_lengthG
xorr1G(1,i)=bitxor(column_imageG(i),mG(i));
end
enc_G=reshape(xorr1G,rowG,colG);
%Blue colure
rgbB=B;
[rowB,colB]=size(rgbB);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%INITIALIZE THE VALUE OF ROTATION%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for j=1:1:rowB
    k(j)=x(j+n);
    l(j)=y(j+p);
end

% for i=1:1:row
%     k(i)=x(i+n);
% end
% for j=1:1:col
%     l(j)=y(j+p);
% end

for j=1:1:colB*rowB
    m1B(j)=z(j+700);
end
mB=uint8(255 * mat2gray(m1B));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% IF k is even right shift row else left shift row%%%%%%%%%%%%%%%%%
% If l is even shift up column else down shift column%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%ROTATION OPERATION%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i=1:1:rowB
    for j=1:1:colB
        if(mod(k(i),2)==0)
            if((j+k(i))<=colB)   %shift right of row
               sh_rowB(i,j+k(i))=rgbB(i,j);
               row_shift_evenB(i,j)=j+k(i);
            else
               sh_rowB(i,(j+k(i)-colB))=rgbB(i,j); 
               row_shift_evenB(i,j)=(j+k(i)-colB);
            end
        else
            if((j-k(i))>=1)       %shift left of row
                sh_rowB(i,j-k(i))=rgbB(i,j);
                row_shift_oddB(i,j)=j-k(i);
            else
                sh_rowB(i,(colB+(j-k(i))))=rgbB(i,j);
                row_shift_oddB(i,j)=colB+j-k(i);
            end
        end
    end
end
   

for j=1:1:colB
    for i=1:1:rowB
        if(mod(l(j),2)==0)
           if((i-l(j))>=1)          %shift up of column
                sh_colB(i-l(j),j)=sh_rowB(i,j);
                col_shift_evenB(i,j)=i-l(j);
            else
                sh_colB((rowB+i-l(j)),j)=sh_rowB(i,j);
                col_shift_evenB(i,j)=rowB+i-l(j);
            end
        else
           if((i+l(j))<=rowB)           %shift down of column
               sh_colB(i+l(j),j)=sh_rowB(i,j);
               col_shift_oddB(i,j)=i+l(j);
               else
               sh_colB((i+l(j)-rowB),j)=sh_rowB(i,j); 
               col_shift_oddB(i,j)=(i+l(j)-rowB);
            end
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%XOR IMAGE%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
total_lengthB=rowB*colB;
column_imageB=reshape(sh_colB,1,total_lengthB);
for i=1:1:total_lengthB
xorr1(1,i)=bitxor(column_imageB(i),mB(i));
end
enc_B=reshape(xorr1,rowB,colB);
encrypted=cat(3,enc_R,enc_G,enc_B);
imshow(encrypted)
imwrite(encrypted,'lena encrypted2.png')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%Decryption XOR%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

