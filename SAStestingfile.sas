

%memore(y = y1 y2, m = m1_1 m1_2 m2_1 m2_2, model = 1, data = parallelserial, save = serial, serial = 1);
%memore(y = y1 y2, m = m1_1 m1_2 m2_1 m2_2 m3_1 m3_2, model = 1, data = parallelserial, save = parallel);

%sobel(a=1.44, sea = 0.222, b = 0.402, seb = 0.070);
%cdfinvt(p=.05, df = 1000);
%dichot(modcount=5, dat = 1);


proc iml;
coeftest = {1,2};
print coeftest;
valuestest = {4,5};
print valuestest;
semattest = {1 2}//{4 5};
print semattest;
temp = .05;
%probres(coef = coeftest, values = valuestest, semat = semattest, df = 100);
quit;

%memore(y = m1 m2, w = dich1 dich2 dich3, model = 2, data = judd);

proc iml;
testdat = {1, 2, 3, 4, 5}||{1, 0, 0, 1, 0};
print testdat;
center = 1;
%centerd(centdat = testdat);
quit;
