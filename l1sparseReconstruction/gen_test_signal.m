function [y,t] = gen_test_signal(type)
clear all;
close all;
% square wave signal
t=0:.005:1;
y0=square(2*pi*0.98*t);
y=y0;
% gaussian blur
% sigma=1;
% size=5;
% x = linspace(-size / 2, size / 2, size);
% gaussFilter = exp(-x .^ 2 / (2 * sigma ^ 2));
% gaussFilter = gaussFilter / sum (gaussFilter); % normalize
% y= filter (gaussFilter,1, y);
% y= conv (y, gaussFilter, 'same');
% add noise
y=awgn(y,15,'measured');
% plot
plot(t,y,'.r-');
hold on;
plot(t,y0,'.g-');
hold off;
axis([-0.2 1.2 -2 2]);
end

