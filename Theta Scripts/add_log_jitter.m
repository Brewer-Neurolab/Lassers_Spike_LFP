function [jitter_x,jitter_y]=add_log_jitter(x,y,jitter)

jitter_x=x;
jitter_y=y;

if numel(x)==numel(y)
    for i=1:length(x)
        jitter_x(i)=10^((rand(1)*jitter*log10(x(i)))+log10(x(i)));
        jitter_y(i)=10^((rand(1)*jitter*log10(y(i)))+log10(y(i)));
    end
end