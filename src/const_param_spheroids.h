﻿//
//*** Physical & Mathematical Constants ***//
double const
    wavel_fix = 0.34,
    pi = 3.1415926535897932384626433832795,
	twopi = 6.283185307179586476925286766559,
	rad = pi/180.0,
	scalef_um3_um2 = 1000.0,
	tiny = 0.0001;
//
//*** Fixed Kernels Parameters ***//
int const
    path_len_max = 256,
	fname_len_max = 64;
char const
//  kernels location
    rootdir[path_len_max] = "./kernels_fix_bin/",
//
	fname_fixkern_ext[fname_len_max] = "Rkext1.fix.bin",
	fname_fixkern_f11[fname_len_max] = "Rkernel1.11.fix.bin",
	fname_fixkern_f22[fname_len_max] = "Rkernel1.22.fix.bin",
	fname_fixkern_f33[fname_len_max] = "Rkernel1.33.fix.bin",
	fname_fixkern_f44[fname_len_max] = "Rkernel1.44.fix.bin",
	fname_fixkern_f12[fname_len_max] = "Rkernel1.12.fix.bin",
	fname_fixkern_f34[fname_len_max] = "Rkernel1.34.fix.bin";
//
//*** Scattering Angle Grid ***//
int const
    nsca_fine = 721, // dtheta = 0.25o - from DLS
	nsca_fix = 181;
double const
	sca_fix[nsca_fix] =
					  {  0.0,   1.0,   2.0,   3.0,   4.0,   5.0,   6.0,   7.0,   8.0,   9.0,
						10.0,  11.0,  12.0,  13.0,  14.0,  15.0,  16.0,  17.0,  18.0,  19.0,
						20.0,  21.0,  22.0,  23.0,  24.0,  25.0,  26.0,  27.0,  28.0,  29.0,
						30.0,  31.0,  32.0,  33.0,  34.0,  35.0,  36.0,  37.0,  38.0,  39.0,
						40.0,  41.0,  42.0,  43.0,  44.0,  45.0,  46.0,  47.0,  48.0,  49.0,
						50.0,  51.0,  52.0,  53.0,  54.0,  55.0,  56.0,  57.0,  58.0,  59.0,
						60.0,  61.0,  62.0,  63.0,  64.0,  65.0,  66.0,  67.0,  68.0,  69.0,
						70.0,  71.0,  72.0,  73.0,  74.0,  75.0,  76.0,  77.0,  78.0,  79.0,
						80.0,  81.0,  82.0,  83.0,  84.0,  85.0,  86.0,  87.0,  88.0,  89.0,
						90.0,  91.0,  92.0,  93.0,  94.0,  95.0,  96.0,  97.0,  98.0,  99.0,
					   100.0, 101.0, 102.0, 103.0, 104.0, 105.0, 106.0, 107.0, 108.0, 109.0,
					   110.0, 111.0, 112.0, 113.0, 114.0, 115.0, 116.0, 117.0, 118.0, 119.0,
					   120.0, 121.0, 122.0, 123.0, 124.0, 125.0, 126.0, 127.0, 128.0, 129.0,
					   130.0, 131.0, 132.0, 133.0, 134.0, 135.0, 136.0, 137.0, 138.0, 139.0,
					   140.0, 141.0, 142.0, 143.0, 144.0, 145.0, 146.0, 147.0, 148.0, 149.0,
					   150.0, 151.0, 152.0, 153.0, 154.0, 155.0, 156.0, 157.0, 158.0, 159.0,
					   160.0, 161.0, 162.0, 163.0, 164.0, 165.0, 166.0, 167.0, 168.0, 169.0,
					   170.0, 171.0, 172.0, 173.0, 174.0, 175.0, 176.0, 177.0, 178.0, 179.0,
					   180.0};
//
//*** Radii Grids ***//
int const
	nrgrid_fix = 34;
double const // nodes are equdistant in log(r)-space. See 'rgrid' below.
	rmin_grid_fix = 0.3306652e-02,
	rmax_grid_fix = 0.2582305e+02,
	rgrid_fix[nrgrid_fix] = {0.3306652e-02, 0.4338570e-02, 0.5692523e-02, 0.7469008e-02, 0.9799887e-02,
                             0.1285817e-01, 0.1687087e-01, 0.2213582e-01, 0.2904382e-01, 0.3810762e-01,
	                         0.5000000e-01, 0.6560367e-01, 0.8607684e-01, 0.1129391e+00, 0.1481844e+00,
	                         0.1944289e+00, 0.2551050e+00, 0.3347165e+00, 0.4391726e+00, 0.5762267e+00,
	                         0.7560518e+00, 0.9919955e+00, 0.1301571e+01, 0.1707757e+01, 0.2240702e+01,
	                         0.2939966e+01, 0.3857452e+01, 0.5061260e+01, 0.6640745e+01, 0.8713145e+01,
	                         0.1143229e+02, 0.1500000e+02, 0.1968110e+02, 0.2582305e+02};
//
//*** Refractive Index Grids: m = n - i*k, k > 0 ***//
int const
    nrefre = 22,
	nrefim = 15;
double const
    refre_min = 1.291429, refre_max = 1.696429,
    refim_min = 0.0005,   refim_max = 0.5,
    refre_fix[nrefre] = {0.1291429E+01, 0.1310714E+01, 0.1330000E+01, 0.1349286E+01, 0.1368571E+01,
                         0.1387857E+01, 0.1407143E+01, 0.1426429E+01, 0.1445714E+01, 0.1465000E+01,
                         0.1484286E+01, 0.1503571E+01, 0.1522857E+01, 0.1542143E+01, 0.1561429E+01,
                         0.1580714E+01, 0.1600000E+01, 0.1619286E+01, 0.1638571E+01, 0.1657857E+01,
                         0.1677143E+01, 0.1696429E+01},
	refim_fix[nrefim] = {0.5000000E-03, 0.8189469E-03, 0.1341348E-02, 0.2196985E-02, 0.3598428E-02,
	                     0.5893843E-02, 0.9653489E-02, 0.1581139E-01, 0.2589737E-01, 0.4241714E-01,
	                     0.6947477E-01, 0.1137923E+00, 0.1863797E+00, 0.3052701E+00, 0.5000000E+00};
/*------------------------------------------------------------------------------------------------*/