clc;  close all;   clear all;
image=imread('download.jpg');
image_hist_RGB_3d(image)
image=imresize(image,[512,382]);
R=image(:,:,1);
G=image(:,:,2);
B=image(:,:,3);
n=40;
p=80;
q=180;

rgbR=R;
[rowR,colR]=size(rgbR);
theta = 90;
z1 = 2.23; z2 = 2.56; z3 = 2.567; z4=2.654; z5=2.543; z6=2.986; z7=2.999; z8=2.543; z9=2.56879;
x(1)  = sin(theta*pi*(z1^2))^2;
y(1) = sin(theta*pi*(z2))^2;
z(1) = sin(theta*pi*(z3))^2;
a(1)  = sin(theta*pi*(z4))^2;
b(1) = sin(theta*pi*(z5))^2;
c(1) = sin(theta*pi*(z6))^2;
d(1)  = sin(theta*pi*(z7))^2;
e(1) = sin(theta*pi*(z8))^2;
f(1) = sin(theta*pi*(z9))^2;

for ii = 2:1:80000
    x(ii) = sin(z1*asin(sqrt(x(ii-1))))^2;
   y(ii) = sin(z2*asin(sqrt(y(ii-1))))^2;
   z(ii) = sin(z3*asin(sqrt(z(ii-1))))^2;
   a(ii) = sin(z4*asin(sqrt(a(ii-1))))^2;
   b(ii) = sin(z5*asin(sqrt(b(ii-1))))^2;
   c(ii) = sin(z6*asin(sqrt(c(ii-1))))^2;
    d(ii) = sin(z7*asin(sqrt(d(ii-1))))^2;
   e(ii) = sin(z8*asin(sqrt(e(ii-1))))^2;
   f(ii) = sin(z9*asin(sqrt(f(ii-1))))^2;
end
x=ceil(mod((x*1000000000000000),256));
 y=ceil(mod((y*100000000000000),256));
 z=ceil(mod((z*100000000000000),256));
 a=ceil(mod((a*1000000000000000),256));
 b=ceil(mod((b*100000000000000),256));
 c=ceil(mod((c*100000000000000),256));
 d=ceil(mod((d*1000000000000000),256));
 e=ceil(mod((e*100000000000000),256));
 f=ceil(mod((f*100000000000000),256));
 for j1=1:1:rowR
    k1(j1)=x(j1+n);
    l1(j1)=y(j1+p);
end

% for i=1:1:row
%     k(i)=x(i+n);
% end
% for j=1:1:col
%     l(j)=y(j+p);
% end

for j1=1:1:colR*rowR
    m1R(j1)=z(j1+q);
end
mR=uint8(255 * mat2gray(m1R));


for i=1:1:rowR
    for j1=1:1:colR
        if(mod(k1(i),2)==0)
            if((j1+k1(i))<=colR)   %shift right of row
               sh_rowG(i,j1+k1(i))=rgbR(i,j1);
               row_shift_evenG(i,j1)=j1+k1(i);
            else
               sh_rowG(i,(j1+k1(i)-colR))=rgbR(i,j1); 
               row_shift_evenG(i,j1)=(j1+k1(i)-colR);
            end
        else
            if((j1-k1(i))>=1)       %shift left of row
                sh_rowG(i,j1-k1(i))=rgbR(i,j1);
                row_shift_oddB(i,j1)=j1-k1(i);
            else
                sh_rowG(i,(colR+(j1-k1(i))))=rgbR(i,j1);
                row_shift_oddB(i,j1)=colR+j1-k1(i);
            end
        end
    end
end
   

for j1=1:1:colR
    for i=1:1:rowR
        if(mod(l1(j1),2)==0)
           if((i-l1(j1))>=1)          %shift up of column
                sh_colG(i-l1(j1),j1)=sh_rowG(i,j1);
                col_shift_evenG(i,j1)=i-l1(j1);
            else
                sh_colG((rowR+i-l1(j1)),j1)=sh_rowG(i,j1);
                col_shift_evenG(i,j1)=rowR+i-l1(j1);
            end
        else
           if((i+l1(j1))<=rowR)           %shift down of column
               sh_colG(i+l1(j1),j1)=sh_rowG(i,j1);
               col_shift_oddB(i,j1)=i+l1(j1);
               else
               sh_colG((i+l1(j1)-rowR),j1)=sh_rowG(i,j1); 
               col_shift_oddB(i,j1)=(i+l1(j1)-rowR);
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
%Green encryption
rgbG=G;
[rowG,colG]=size(rgbG);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%INITIALIZE THE VALUE OF ROTATION%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for j2=1:1:rowG
    k2(j2)=a(j2+n);
    l2(j2)=b(j2+p);
end

% for i=1:1:row
%     k(i)=x(i+n);
% end
% for j=1:1:col
%     l(j)=y(j+p);
% end

for j2=1:1:colG*rowG
    m1G(j2)=c(j2+100);
end
mG=uint8(255 * mat2gray(m1G));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% IF k is even right shift row else left shift row%%%%%%%%%%%%%%%%%
% If l is even shift up column else down shift column%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%ROTATION OPERATION%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i=1:1:rowG
    for j2=1:1:colG
        if(mod(k2(i),2)==0)
            if((j2+k2(i))<=colG)   %shift right of row
               sh_rowG(i,j2+k2(i))=rgbG(i,j2);
               row_shift_evenG(i,j2)=j2+k2(i);
            else
               sh_rowG(i,(j2+k2(i)-colG))=rgbG(i,j2); 
               row_shift_evenG(i,j2)=(j2+k2(i)-colG);
            end
        else
            if((j2-k2(i))>=1)       %shift left of row
                sh_rowG(i,j2-k2(i))=rgbG(i,j2);
                row_shift_oddB(i,j2)=j2-k2(i);
            else
                sh_rowG(i,(colG+(j2-k2(i))))=rgbG(i,j2);
                row_shift_oddB(i,j2)=colG+j2-k2(i);
            end
        end
    end
end
   

for j2=1:1:colG
    for i=1:1:rowG
        if(mod(l2(j2),2)==0)
           if((i-l2(j2))>=1)          %shift up of column
                sh_colG(i-l2(j2),j2)=sh_rowG(i,j2);
                col_shift_evenG(i,j2)=i-l2(j2);
            else
                sh_colG((rowG+i-l2(j2)),j2)=sh_rowG(i,j2);
                col_shift_evenG(i,j2)=rowG+i-l2(j2);
            end
        else
           if((i+l2(j2))<=rowG)           %shift down of column
               sh_colG(i+l2(j2),j2)=sh_rowG(i,j2);
               col_shift_oddB(i,j2)=i+l2(j2);
               else
               sh_colG((i+l2(j2)-rowG),j2)=sh_rowG(i,j2); 
               col_shift_oddB(i,j2)=(i+l2(j2)-rowG);
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
%Blue encryption
rgbB=B;
[rowB,colB]=size(rgbB);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%INITIALIZE THE VALUE OF ROTATION%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for j3=1:1:rowB
    k3(j3)=d(j3+n);
    l3(j3)=e(j3+p);
end

% for i=1:1:row
%     k(i)=x(i+n);
% end
% for j=1:1:col
%     l(j)=y(j+p);
% end

for j3=1:1:colB*rowB
    m1B(j3)=f(j3+700);
end
mB=uint8(255 * mat2gray(m1B));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% IF k is even right shift row else left shift row%%%%%%%%%%%%%%%%%
% If l is even shift up column else down shift column%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%ROTATION OPERATION%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i=1:1:rowB
    for j3=1:1:colB
        if(mod(k3(i),2)==0)
            if((j3+k3(i))<=colB)   %shift right of row
               sh_rowB(i,j3+k3(i))=rgbB(i,j3);
               row_shift_evenB(i,j3)=j3+k3(i);
            else
               sh_rowB(i,(j3+k3(i)-colB))=rgbB(i,j3); 
               row_shift_evenB(i,j3)=(j3+k3(i)-colB);
            end
        else
            if((j3-k3(i))>=1)       %shift left of row
                sh_rowB(i,j3-k3(i))=rgbB(i,j3);
                row_shift_oddB(i,j3)=j3-k3(i);
            else
                sh_rowB(i,(colB+(j3-k3(i))))=rgbB(i,j3);
                row_shift_oddB(i,j3)=colB+j3-k3(i);
            end
        end
    end
end
   

for j3=1:1:colB
    for i=1:1:rowB
        if(mod(l3(j3),2)==0)
           if((i-l3(j3))>=1)          %shift up of column
                sh_colB(i-l3(j3),j3)=sh_rowB(i,j3);
                col_shift_evenB(i,j3)=i-l3(j3);
            else
                sh_colB((rowB+i-l3(j3)),j3)=sh_rowB(i,j3);
                col_shift_evenB(i,j3)=rowB+i-l3(j3);
            end
        else
           if((i+l3(j3))<=rowB)           %shift down of column
               sh_colB(i+l3(j3),j3)=sh_rowB(i,j3);
               col_shift_oddB(i,j3)=i+l3(j3);
               else
               sh_colB((i+l3(j3)-rowB),j3)=sh_rowB(i,j3); 
               col_shift_oddB(i,j3)=(i+l3(j3)-rowB);
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
image_encryption=xorr1;

enc_B=reshape(xorr1,rowB,colB);
encrypted=cat(3,enc_R,enc_G,enc_B);
imshow(encrypted)
imwrite(encrypted,'encrypted.png')