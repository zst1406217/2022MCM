t0=zeros(1,100);
x0=zeros(1,100);
t0(1)=0;
dx_t0=22050/99;
size_x=size(x);
j0_=zeros(1,100);
for i=95:99
    j0=0;
    for j=1:size_x(1)-1
        if (x(j)<i*dx_t0)&&(x(j+1)>=i*dx_t0)
            j0=j;
            break;
        end
    end
    j0_(i+1)=j0
    a1_=x(j0)
    a2_=x(j0+1)
    b1_=t(j0)
    b2_=t(j0+1)
    t0(i+1)=(b2_-b1_)*(i*dx_t0-a1_)/(a2_-a1_)+b1_
    x0(i+1)=i*dx_t0;
end