// 4-point, 3rd-order Hermite (z-form)
float z = x - 1/2.0;
float even1 = y[-1]+y[2], odd1 = y[-1]-y[2];
float even2 = y[0]+y[1], odd2 = y[0]-y[1];
float c0 = 9/16.0*even2 - 1/16.0*even1;
float c1 = 1/8.0*odd1 - 11/8.0*odd2;
float c2 = 1/4.0*(even1-even2);
float c3 = 3/2.0*odd2 - 1/2.0*odd1;
return ((c3*z+c2)*z+c1)*z+c0;
