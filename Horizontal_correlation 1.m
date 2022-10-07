
Image = imread('azan enc.jpg');      % Select any image to check Horizontal correlation
x = Image(:,1:end-1);             % All rows and columns 1 through 255 
y = Image(:,2:end);               % All rows and columns 2 through 256
randIndex = randperm(numel(x));   % A random permutation of the integer from 1 to numel(x)
randIndex = randIndex(1:1500);    % 1000 random Horizontal Pairs
xRand = x(randIndex);             % 1000 random values from x
yRand = y(randIndex);             % The corresponding 1000 values from y

scatter(xRand,yRand,5,'filled')
corr2(xRand,yRand)                % numerical value of correlation coefficient