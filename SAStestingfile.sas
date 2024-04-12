

%memore(y = y1 y2, m = m11 m12 m21 m22, model = 1, data = parallelserial, save = serial, serial = 1);
%memore(y = y1 y2, m = m11 m12 m21 m22 m31 m32, model = 1, data = parallelserial, save = parallel);

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
*testdat = {1, 2, 3, 4, 5}||{1, 0, 0, 1, 0};
testdat = {1, 0, 0, 1, 0};
print testdat;
center = 2;
%centerd(centdat = testdat);
quit;

%memore(y = y1 y2, w = m11 m12, model = 3, data = parallelserial);

proc iml;
moddat = {-1,1, 2, 3, 4, 5,600};
*testdat = {1, 0, 0, 1, 0};
N = nrow(moddat);
print moddat;
%JNprobe(coefone = 1, coefTwo = 2, seOne = 3, seTwo = 4, cov = -10, critt = 1.96, dfJN = 100);
quit;


%memore(y = y1 y2, m = m11 m12, w = m21, model = 4, center = 2, plot = 1, data = parallelserial);
