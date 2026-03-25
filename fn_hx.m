%% h(X):Nonlinear measurement eqn
function h = fn_hx(xr,yr,x_)
x = x_(1);
y = x_(4);
h = [
    sqrt((x-xr)^2 + (y-yr)^2)
    atan2(y-yr,x-xr)
    ];