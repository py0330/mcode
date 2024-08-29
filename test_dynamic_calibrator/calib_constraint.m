function [c,ceq] = calib_constraint(x)

c = [];

idx = 0;

mass = x(idx + 1);
mcx  = x(idx + 2);
mcy  = x(idx + 3);
mcz  = x(idx + 4);
Ixx  = x(idx + 5);
Iyy  = x(idx + 6);
Izz  = x(idx + 7);
Ixy  = x(idx + 8);
Ixz  = x(idx + 9);
Iyz  = x(idx + 10);

c = [   c
        ...(mcy*mcy + mcz*mcz) - Ixx*mass;
        ...(mcx*mcx + mcz*mcz) - Iyy*mass;
        (mcx*mcx + mcy*mcy) - Izz*mass;
    ];

idx = 10;
mass = x(idx + 1);
mcx  = x(idx + 2);
mcy  = x(idx + 3);
mcz  = x(idx + 4);
Ixx  = x(idx + 5);
Iyy  = x(idx + 6);
Izz  = x(idx + 7);
Ixy  = x(idx + 8);
Ixz  = x(idx + 9);
Iyz  = x(idx + 10);

c = [   c
        ...(mcy*mcy + mcz*mcz) - Ixx*mass;
        ...(mcx*mcx + mcz*mcz) - Iyy*mass;
        (mcx*mcx + mcy*mcy) - Izz*mass;
    ];

ceq = [];

end