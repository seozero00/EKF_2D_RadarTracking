%% partial_h/partial_X : Jacobian of measurement eq
function H = fn_JH(xr,yr,x_)
x = x_(1);
y = x_(4);
r = sqrt((x-xr)^2+(y-yr)^2);
temp = cos(atan2(y-yr,x-xr))^2;
H = [
    (x-xr)/r 0 0 (y-yr)/r 0 0;
    -temp*((y-yr)/(x-xr)^2) 0 0 temp/(x-xr) 0 0
    ];