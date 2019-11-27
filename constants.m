function [phys] = constants()

phys = struct();

phys.('eps0') = 8.85418782e-12;
phys.('u0') = 1.25663706e-6;
phys.('c') = 2.99792458e8;
phys.('e') = 1.60217646d-19;
phys.('amu') = 1.66053904d-27;
phys.('me') = 9.10938188d-31;
phys.('mH') = 1.00794*phys.('amu');
phys.('me_amu') = phys.('me') / phys.('amu');

end
