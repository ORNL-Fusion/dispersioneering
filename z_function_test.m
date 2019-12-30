function tests = z_function_test

tests = functiontests(localfunctions);

end

function test_z_function(testCase)

% Plot for sanity
%
% z_fun_arg_re   = linspace(-10,10,100);
% z_fun_arg_re_C = z_fun_arg_re * 0;
% 
% [Z,Zp] = z_function(z_fun_arg_re, z_fun_arg_re_C);
% 
% subplot(2,1,1)
% p=plot(z_fun_arg_re,real(Z));
% hold on
% plot(z_fun_arg_re,imag(Z));
% hold off
% 
% subplot(2,1,2)
% plot(z_fun_arg_re,real(Zp));
% hold on
% plot(z_fun_arg_re,imag(Zp));
% hold off



% Check against specfic inputs with the mathematica

arg_re = 0.8;
[act_Z,act_Zp] = z_function(arg_re,0);

exp_Z = -1.064203414112731 + 0.9346014875484057i;
exp_Zp = -0.2972745374196304-1.495362380077449i;

verifyEqual(testCase,act_Z,exp_Z,'RelTol',1e-8)
verifyEqual(testCase,act_Zp,exp_Zp,'RelTol',1e-8)



arg_re = 1.8;
[act_Z,act_Zp] = z_function(arg_re,0);

exp_Z = -0.6935455382297442+0.06941619668465928i;
exp_Zp = 0.4967639376270792-0.2498983080647735i;

verifyEqual(testCase,act_Z,exp_Z,'RelTol',1e-7)
verifyEqual(testCase,act_Zp,exp_Zp,'RelTol',1e-7)



arg_re = 25.0;
[act_Z,act_Zp] = z_function(arg_re,0);

exp_Z = -0.04003207710893372-0i;
exp_Zp = 0.001603855446686175+0i;

verifyEqual(testCase,act_Z,exp_Z,'RelTol',1e-7)
verifyEqual(testCase,act_Zp,exp_Zp,'RelTol',1e-7)

end
