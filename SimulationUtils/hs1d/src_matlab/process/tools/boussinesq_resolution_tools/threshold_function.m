function OUT=threshold_function(x)
% sigmoid function centered on 1 so it is concave
%     OUT=1-2./(1+exp(-1000000*(x-1)));
%     OUT=1-2./(1+exp(-1000*(x-1)));
%     OUT=1-exp(-100000*(1-x));
    OUT=1-exp(-1000*(1-x));
%     OUT=1*ones(size(x));
%     OUT=1-1./(1+exp(-1000*(x-0.99)));   
% OUT=x<1;
%     OUT=(1-1./(1+exp(-600*(x-0.97))));
%     OUT=(1-1./(1+exp(-70*(x-0.8))));
end