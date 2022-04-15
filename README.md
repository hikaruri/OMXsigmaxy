# OMXsigmaxy

This is a post-processing code for computing anomalous Hall conductivity in [OpenMX](http://openmx-square.org/) by expanded Fukui-Hatsugai-Suzuki method [\[J. Phys. Soc. Jpn. **74**, 1674 (2005)\]](https://journals.jps.jp/doi/abs/10.1143/JPSJ.76.053702). This code will be included in OpenMX (Original Fukui-Hatsugai-Suzuki method was implemented in [OpenMX3.9](http://www.openmx-square.org/openmx_man3.9/node184.html)).
This code can compute the metallic system efficiently, so this is suitable for thermoelectric material design based on the high-throughput first-principles screening.

Please cite the following article:
- Hikaru Sawahata, Naoya Yamaguchi, Susumu Minami, and Fumiyuki Ishii, First-principles calculation of anomalous Hall and Nernst conductivity by local Berry phase, [arXiv:2204.05949](https://arxiv.org/abs/2204.05949) (2022).

# Installation (how to make)

Please see the [OpenMX manual](http://www.openmx-square.org/openmx_man3.9/node4.html).
If you success the installation and compiling of OpenMX, you have to use same compiler option by specifying CC, FC and LIB in makefile.
If you use [MateriApps LIVE!](http://cmsi.github.io/MateriAppsLive/) environment, you can use our code by specifying the compiler option as:
```
 CC = mpicc -O3 -fopenmp -I/usr/include
 LIB = -L/usr/lib -llapack -lblas -lgfortran -lmpi_mpifh
```
and please type as
```
$ make
```
on your terminal. An executable file 'sigmaxy' will be generated.

# How to excute

Please prepare '.scfout' file in OpenMX by turning on SCF calculation's option as
```shell script
HS.fileout   on
```
Please refere to [OpenMX manual](http://www.openmx-square.org/openmx_man3.9/node213.html) and [How to use calB.c](http://www.openmx-square.org/openmx_man3.9/node185.html).

And, you can execute 'sigmaxy' as
```shell script
($PATH)/sigmaxy ($SystemName).scfout > (Standard output file)
```
or
```shell script
($PATH)/sigmaxy ($SystemName).scfout < sigmaxy.in > (Standard output file)
```

After carrying out calculation, you can obtain '($SystemName).scfout.AHC.dat'.
This file stores the anomalous Hall conductivity on every chemical energy as:
```shell script
# Chemical Energy (eV)  sigma_yz   sigma_zx   sigma_xy
...
...
```
# Input parameter
Please specify the input parameters in 'sigmaxy.in' as:
```
Elow Ehigh Emesh
Mesh_Number
AHC_Unit ([0]S/cm or [1]e^2/h)
[0]2D or [1]3D?
```
For example, in Febcc case,
```
-1.0 1.0 400
200 200 200
0
1
```
this code compute a chemical potential dependence of 
\sigma_{xy} -0.5 eV to 0.5 eV
(energy mesh number is 400 and descritized integration mesh is 100\*100\*100).

When you choose a 2D option, sigmaxy compute AHC of k_3 surface (sigma_{xy}).
