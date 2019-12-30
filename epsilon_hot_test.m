function tests = epsilon_hot_test

tests = functiontests(localfunctions);

end

function test_epsilon_hot(testCase)

f=42e6;
amu=2;
Z=1;
B=1;
n=3.677e19;
T_eV = 100;
k_per = 10;
k_par = 100;

[actSolution_eps_cold,actSolution_sig_cold] = epsilon_cold(f,amu,Z,B,n);
[actSolution_eps,actSolution_sig] = epsilon_hot(f,amu,Z,B,n,T_eV,k_per,k_par);

expSolution_sig = complex(zeros(3,3));

S = 0 + 1i*1.1142493;
D = -0.20369692 + 1i*0;
P = 0 + 1i*1.0770133;

expSolution_sig(1,1) = S;
expSolution_sig(1,2) = D;
expSolution_sig(1,3) = 0;

expSolution_sig(2,1) = -D;
expSolution_sig(2,2) = S;
expSolution_sig(2,3) = 0;

expSolution_sig(3,1) = 0;
expSolution_sig(3,2) = 0;
expSolution_sig(3,3) = P;

verifyEqual(testCase,actSolution_sig,expSolution_sig,'RelTol',1e-4)

end