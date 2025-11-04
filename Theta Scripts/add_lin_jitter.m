function [jitter_x,jitter_y]=add_lin_jitter(x,y,jitter)

jitter_x=x;
jitter_y=y;

if numel(x)==numel(y)
    for i=1:length(x)
        jitter_x(i)=(rand(1)*jitter*(x(i)))+(x(i));
        jitter_y(i)=(rand(1)*jitter*(y(i)))+(y(i));
    end
end