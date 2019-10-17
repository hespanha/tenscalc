function y=relu(x)
    y=x;
    y(y<0)=0;
end