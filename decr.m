function decryption=ima_ge(A)
image=A;
R=image(:,:,1);
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
image=A;
R=image(:,:,1);
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
[rowR,colR]=size(rgbR);
total_length=rowR*colR;
column_image=reshape(rgbR,1,total_length);
%Red encryption
for i=1:1:total_length
        xorr1R(1,i)=bitxor(column_image(i),mR(i));
end
shuffledR=reshape(xorr1R,rowR,colR);
 for j=1:1:colR
    for i=1:1:rowR
        if(mod(l(j),2)==0)
           if((i+l(j))<=rowR)           %shift down of column
               sh_col(i+l(j),j)=shuffledR(i,j);
               col_shift_even(i,j)=i+l(j);
               else
               sh_col((i+l(j)-rowR),j)=shuffledR(i,j); 
               col_shift_even(i,j)=(i+l(j)-rowR);
            end
        else
           if((i-l(j))>=1)          %shift up of column
                sh_col(i-l(j),j)=shuffledR(i,j);
                col_shift_odd(i,j)=i-l(j);
            else
                sh_col((rowR+i-l(j)),j)=shuffledR(i,j);
                col_shift_odd(i,j)=rowR+i-l(j);
            end
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


for i=1:1:rowR
    for j=1:1:colR
        if(mod(k(i),2)==0)
            if((j-k(i))>=1)       %shift left of row
                sh_rowR(i,j-k(i))=sh_col(i,j);
                row_shift_even(i,j)=j-k(i);
            else
                sh_rowR(i,(colR+j-k(i)))=sh_col(i,j);
                row_shift_even(i,j)=colR+j-k(i);
            end
        else
            if((j+k(i))<=colR)   %shift right of row
               sh_rowR(i,j+k(i))=sh_col(i,j);
               row_shift_even(i,j)=j+k(i);
            else
               sh_rowR(i,(j+k(i)-colR))=sh_col(i,j); 
               row_shift_even(i,j)=(j+k(i)-colR);
            end
        end
    end
end
