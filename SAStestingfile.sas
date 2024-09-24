

%memore(y = y1 y2, m = m11 m12 m21 m22, model = 1, data = parallelserial, save = serial, serial = 1);
%memore(y = y1 y2, m = m11 m12 m21 m22 m31 m32, model = 1, data = parallelserial, save = parallel);

proc iml;
aval = 1.44;
seaval = 0.222;
bval = 0.402;
sebval = 0.070;
%sobel(a=aval, sea = seaval, b = bval, seb = sebval);
print sobelres;
quit;
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

proc iml;
testdat = {-1,1, 2, 3, 4, 5,600,1, 2, 3, 4, 5,6001, 2, 3, 4, 5,6001, 2, 3, 4, 5,6001, 2, 3, 4, 5,6001, 2, 3, 4, 5,600};
N = nrow(moddat);
print moddat;
%JNprobe(coefone = 1, coefTwo = 2, seOne = 3, seTwo = 4, cov = -10, critt = 1.96, dfJN = 100);
quit;


%memore(y = y1 y2, m = m11 m12 m21 m22 m31 m32, serial = 1, model = 1, normal = 1, mc = 1, xmint = 0, data = parallelserial);

%memore(y = y1 y2, m = m11 m12 m31 m32, w = m21, normal = 1, mc = 1, model = 5, contrast = 1, samples = 1000, data = parallelserial);

%memore(y = y1 y2, m = m11 m12 m31 m32, w = m21, model = 4, plot = 1, jn=1, xmint = 0, samples = 1000, data = parallelserial, save = est3);
%memore(y = y1 y2, m = m11 m12 m31 m32, w = m21, model = 4, plot = 1, jn=1, xmint = 0, mc = 1, samples = 1000, data = parallelserial, save = est4);
%memore(y = y1 y2, m = m11 m12 m31 m32, w = m21, model = 7, plot = 1, jn=1, samples = 1000, data = parallelserial);
%memore(y = y1 y2, m = m11 m12 m31 m32, w = m21, model = 8, plot = 1, jn=1, samples = 1000, data = parallelserial);

%memore(y = Y1 Y2, m = M1_1 M1_2 M2_1 M2_2 M3_1 M3_2 M4_1 M4_2 M5_1 M5_2, model = 1, normal = 1, mc = 1, contrast = 1, data = work.FakeSerialData);

%memore(y = y1 y2, w = m21 m22, model = 3, plot = 1, jn = 1, data = parallelserial);

%memore(y = y1 y2, m = m11 m12, model = 1, data = parallelserial);
%memore(y = y1 y2, m = m11 m12, model = 1, normal = 1, mc = 1, data = parallelserial);
%memore(y = y1 y2, m = m11 m12 m21 m22, model = 1, normal = 1, mc = 1, data = parallelserial);

%memore(y = y1 y2, m = m11 m12 m31 m32, w = m21, model = 5, samples = 1000, data = parallelserial);

%memore(y = y1 y2, m = m11 m12 m31 m32, model = 1, normal = 1, samples = 1000, data = parallelserial);
%memore(y = y1 y2, m = m11 m12, model = 1, normal = 1, samples = 1000, data = parallelserial);

%memore(y = int_g int_I, w = order, model = 2, data = CompSci_WS);
%memore(w=order,y = int_G int_I, model = 2, plot = 1, data = CompSci_WS);
%memore(w=grppref order,y = int_G int_I, model = 2, 
        data = CompSci_WS);

%memore(y = int_g int_I, m = comm_g comm_i diff_g diff_i, model = 1, bc = 1, data = CompSci_WS, samples = 1000);
%memore(y = int_g int_I, m = comm_g comm_i diff_g diff_i, w = order, model = 4, data = CompSci_WS, samples = 20000);
%memore(y = int_g int_I, m = comm_g comm_i diff_g diff_i, w = grppref, model = 12, data = CompSci_WS);


%memore(y = int_g int_I, m = comm_g comm_i diff_g diff_i, w = order, model = 4, bc = 1, samples = 10000, data = CompSci_WS);
