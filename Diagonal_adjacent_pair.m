
Image = imread('azan1.jpg');     % Select any image to check Diagonal correlation
x = Image(1:end-1,1:end-1);      % All but the last row and column
y = Image(2:end,2:end);          % All but the first row and column
randIndex = randperm(numel(x));  % A random permutation of the integer from 1 to numel(x) 
randIndex = randIndex(1:1500);   % 1000 random Diagonal Pairs
xRand = x(randIndex);            % 1000 random values from x
yRand = y(randIndex);            % The corresponding 1000 values from y

scatter(xRand,yRand,5,'filled')       

corr2(xRand,yRand)              % numerical value of correlation coefficient
% axis ([0 270 0 270])