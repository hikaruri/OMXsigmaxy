/**********************************************************************
  sigmaxy.c

  Code for calculating anomalous Hall conductivity.

  30/Jul/2019    Rename calB.c -> sigmaxy.c by Hikaru Sawahata
  04/Sep/2019    Modified for efficiency and stability by Naoya Yamaguchi
  13/Jun/2021    Tuned for acceleration by Naoya Yamaguchi
  12/Nov/2021    Modified for the maximum band indices by Naoya Yamaguchi
  21/Nov/2021    Modified for a memory leak by Naoya Yamaguchi
  22/Nov/2021    Modified for k-point indices by Naoya Yamaguchi
  11/Apr/2022    Modified for release by Hikaru Sawahata


  polB.c:

  code for calculating the electric polarization of bulk system using
  Berry's phase.

  Log of polB.c:

     30/Nov/2005   Released by Taisuke Ozaki
     27/Feb/2006   Modified by Fumiyuki Ishii
     28/July/2006  Modified for MPI by Fumiyuki Ishii
     19/Jan/2007   Modified by Taisuke Ozaki
***********************************************************************/

#include "f77func.h"
#include "lapack_prototypes.h"
#include "mpi.h"
#include "read_scfout.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <sys/times.h>
#include <sys/types.h>
#include <time.h>

#define Host_ID 0 /* ID of the host CPU in MPI */

#define printout 0 /* 0:off, 1:on */
#define PI 3.1415926535897932384626
#define measure_time 0
#define dste_flag 2
#define AU2Debye 2.54174776
#define AU2Mucm 5721.52891433 /* 100 e/bohr/bohr */
#define Bohr 0.529177

/* Added by N. Yamaguchi ***/
#define HartreeEV 27.2113845
/* ***/

#define e2parh 3.874046e-5      /* S       */
#define Boltzmankb 8.617333e-5  /* eV / K  */
#define ElemCharge 1.602177e-19 /* A * sec */
#define EVtoJ 1.602177e-19

#define ALLOCATE(TYPE, NAME, SIZE)                                             \
  TYPE *NAME = (TYPE *)malloc(sizeof(TYPE) * SIZE)
#define DEALLOCATE(NAME)                                                       \
  free(NAME);                                                                  \
  NAME = NULL
#define PRINTFF(...)                                                           \
  do {                                                                         \
    printf(__VA_ARGS__);                                                       \
    fflush(stdout);                                                            \
  } while (0)

struct timeval2 {
  long tv_sec;  /* second */
  long tv_usec; /* microsecond */
};

/* Added by N. Yamaguchi ***/
dcomplex ****expOLP;
dcomplex *****expOLP_HWF;
static void expOLP_Band(dcomplex ****expOLP, dcomplex **T, int *MP, double k1,
                        double k2, double k3, double dka, double dkb,
                        double dkc, char sign);
static void setexpOLP(double dka, double dkb, double dkc, int calcOrderMax,
                      dcomplex ****expOLP);
static dcomplex ****memoryAllocation_dcomplex(dcomplex ****in);
/* ***/

/* Modified by N. Yamaguchi ***/
static void Overlap_k1k2(int diag_flag, double k1[4], double k2[4],
                         int spinsize, int fsize, int fsize2, int fsize3,
                         int calcBandMin, int calcBandMax, int calcOrderMax,
                         int *MP, dcomplex ***Sop, dcomplex ***Wk1,
                         dcomplex ***Wk2, double **EigenVal1,
                         double **EigenVal2);
/* ***/

/* For Black box */
void Calc_OverlapSopk1234(double k1[4], double k2[4], double k3[4],
                          double k4[4], dcomplex ***Wk1, dcomplex ***Wk2,
                          dcomplex ***Wk3, dcomplex ***Wk4, double **EigenVal1,
                          double **EigenVal2, double **EigenVal3,
                          double **EigenVal4, dcomplex ***Sop12,
                          dcomplex ***Sop23, dcomplex ***Sop34,
                          dcomplex ***Sop41, int spinsize, int fsize,
                          int fsize2, int fsize3, int fsize4, int *MP);

void Calc_BerryCEps(double **BerryCproc, int ABloop, int Ecount,
                    double FermiEnergy, double ChemP, int fsize2,
                    double **EigenVal1, double **EigenVal2, double **EigenVal3,
                    double **EigenVal4, dcomplex ***Sop12, dcomplex ***Sop23,
                    dcomplex ***Sop34, dcomplex ***Sop41, dcomplex *work1,
                    dcomplex *work2, dcomplex Cdet[2], INTEGER *ipiv);

static void Eigen_HH(dcomplex **ac, double *ko, int n, int EVmax);
static void lapack_dstevx2(INTEGER N, double *D, double *E, double *W,
                           dcomplex **ev);
static void Overlap_Band(double ****OLP, dcomplex **S, int *MP, double k1,
                         double k2, double k3);
static void Hamiltonian_Band(double ****RH, dcomplex **H, int *MP, double k1,
                             double k2, double k3);
static void Hamiltonian_Band_NC(double *****RH, double *****IH, dcomplex **H,
                                int *MP, double k1, double k2, double k3);
static void determinant(int spin, int N, dcomplex **a, INTEGER *ipiv,
                        dcomplex *a1d, dcomplex *work, dcomplex Cdet[2]);

static void dtime(double *t);
static double tmp[4];

void LU_fact(int n, dcomplex **a);
dcomplex RCdiv(double x, dcomplex a);
dcomplex Csub(dcomplex a, dcomplex b);
dcomplex Cmul(dcomplex a, dcomplex b);
dcomplex Cdiv(dcomplex a, dcomplex b);

void Cross_Product(double a[4], double b[4], double c[4]);
double Dot_Product(double a[4], double b[4]);

int main(int argc, char *argv[]) {
  int fsize, fsize2, fsize3, fsize4;
  int spin, spinsize, diag_flag;
  int i, j, k, hos, po, wan, hog2, valnonmag;
  int n1, n2, n3, i1, i2, i3;
  int Nk[4], kloop[4][4];
  int pflag[4];
  int hog[1];
  int metal;
  double k1[4], k2[4], k3[4], k4[4];
  double Gabs[4], Parb[4];
  double Phidden[4], Phidden0;
  double tmpr, tmpi;
  double mulr[2], muli[2];
  double **kg;
  double sum, d;
  double norm;
  double pele, piony, pall;
  double gdd;
  double CellV;
  double Cell_Volume;
  double psi, sumpsi[2];
  double detr;
  double TZ;
  int ct_AN, h_AN;
  char *s_vec[20];
  int *MP;
  double **EigenVal1;
  double **EigenVal2;
  double **EigenVal3;
  double **EigenVal4;
  dcomplex ***Sop12;
  dcomplex ***Sop23;
  dcomplex ***Sop34;
  dcomplex ***Sop41;
  dcomplex ***Wk1;
  dcomplex ***Wk2;
  dcomplex ***Wk3;
  dcomplex ***Wk4;
  MPI_Comm comm1;
  int numprocs, myid, ID, ID1;
  double TStime, TEtime;

  dcomplex *work1;  /* for determinant */
  dcomplex *work2;  /* for determinant */
  dcomplex Cdet[2]; /* for determinant */
  INTEGER *ipiv;    /* for determinant */
  /* for MPI*/
  int *arpo;
  int AB_knum, S_knum, E_knum, num_ABloop0;
  int AB_knum0;
  int ik1, ik2, ik3, ABloop, ABloop0, abcount;
  double tmp4;

  /* MPI initialize */

  MPI_Status stat;
  MPI_Request request;
  MPI_Init(&argc, &argv);
  comm1 = MPI_COMM_WORLD;
  MPI_Comm_size(comm1, &numprocs);
  MPI_Comm_rank(comm1, &myid);

  dtime(&TStime);
  if (myid == Host_ID) {
    printf("\n*****************************************************************"
           "*\n");
    printf(
        "******************************************************************\n");
    printf(" sigmaxy:\n");
    printf(" code for calculating the Hall conductivity every chemical "
           "potential\n");
    printf(" by expanded Fukui-Hatsugai-Suzuki method\n");
    printf(" [J. Phys. Soc. Jpn. 74, 1674 (2005)].\n");
    printf(" Copyright (C), 2006-2022, Hikaru Sawahata, Naoya Yamaguchi,\n");
    printf(" Susumu Minami, Fumiyuki Ishii and Taisuke Ozaki\n");
    printf(" This is free software, and you are welcome to         \n");
    printf(" redistribute it under the constitution of the GNU-GPL.\n");
    printf(" Please cite the following article:\n");
    printf(" H. Sawahata, N. Yamaguchi, S. Minami and F. Ishii,\n");
    printf(" arXiv:2204.05949 (2021).\n");
    printf(
        "******************************************************************\n");
    printf(
        "******************************************************************\n");
  }

  /* ****************************************
   *          scan Mesh,BandNumber
   * ****************************************/
  int mesh1, mesh2, mesh3;
  int EnergyMesh;
  double Elow, Ehigh, FermiEnergy, Eunit;
  int Unitflag, Dimensionflag;

  if (myid == Host_ID) {
    printf("Energy Range (eV) : Elow Ehigh Emesh\n");
    scanf("%lf %lf %d", &Elow, &Ehigh, &EnergyMesh);
    printf("Mesh Number:");
    scanf("%d %d %d", &mesh1, &mesh2, &mesh3);
    printf("AHC Unit ([0]S/cm or [1]e^2/h):");
    scanf("%d", &Unitflag);
    printf("[0]2D or [1]3D?:");
    scanf("%d", &Dimensionflag);
  }

  MPI_Barrier(comm1);

  MPI_Bcast(&mesh1, 1, MPI_INT, 0, comm1);
  MPI_Bcast(&mesh2, 1, MPI_INT, 0, comm1);
  MPI_Bcast(&mesh3, 1, MPI_INT, 0, comm1);

  MPI_Bcast(&hog2, 1, MPI_INT, 0, comm1);

  MPI_Bcast(&Elow, 1, MPI_DOUBLE, 0, comm1);
  MPI_Bcast(&Ehigh, 1, MPI_DOUBLE, 0, comm1);
  MPI_Bcast(&EnergyMesh, 1, MPI_INT, 0, comm1);

  MPI_Bcast(&Dimensionflag, 1, MPI_INT, 0, comm1);

  Eunit = (Ehigh - Elow) / EnergyMesh;

  int Ecount;

  double *sigmaxyProc;
  double *sigmaxySum;
  double *sigmaxy[4];

  sigmaxyProc = (double *)malloc(sizeof(double) * EnergyMesh);
  sigmaxySum = (double *)malloc(sizeof(double) * EnergyMesh);

  for (i = 0; i < 4; i++) {
    sigmaxy[i] = (double *)malloc(sizeof(double) * EnergyMesh);
  }

  for (i = 0; i < EnergyMesh; i++) {
    sigmaxy[1][i] = 0.0;
    sigmaxy[2][i] = 0.0;
    sigmaxy[3][i] = 0.0;
  }

  double kx, ky;

  double Kspace;
  double A12, A23, A34, A41;
  /******************************************
              read a scfout file
  ******************************************/

  read_scfout(argv);

  if (SpinP_switch != 3) {
    PRINTFF("Error: \"sigmaxy\" is available only for non-collinear cases.\n");
    MPI_Finalize();
    return 0;
  }

  s_vec[0] = "Recursion";
  s_vec[1] = "Cluster";
  s_vec[2] = "Band";
  s_vec[3] = "NEGF";
  s_vec[4] = "DC";
  s_vec[5] = "GDC";
  s_vec[6] = "Cluster2";
  s_vec[7] = "Krylov";

  if (myid == Host_ID) {
    printf(" Previous eigenvalue solver = %s\n", s_vec[Solver - 1]);
    printf(" atomnum                    = %i\n", atomnum);
    printf(" ChemP                      = %15.12f (Hartree)\n", ChemP);
    printf(" E_Temp                     = %15.12f (K)\n", E_Temp);
    printf(" Total_SpinS                = %15.12f (K)\n", Total_SpinS);
  }

  s_vec[0] = "collinear spin-unpolarized";
  s_vec[1] = "collinear spin-polarized";
  s_vec[3] = "non-collinear";

  if (myid == Host_ID) {
    printf(" Spin treatment             = %s\n", s_vec[SpinP_switch]);
  }

  if (myid == Host_ID) {
    /* printf("\n fsize4=%d\n",fsize4); */
    printf("\n r-space primitive vector (Bohr)\n");
    printf("  tv1=%10.6f %10.6f %10.6f\n", tv[1][1], tv[1][2], tv[1][3]);
    printf("  tv2=%10.6f %10.6f %10.6f\n", tv[2][1], tv[2][2], tv[2][3]);
    printf("  tv3=%10.6f %10.6f %10.6f\n", tv[3][1], tv[3][2], tv[3][3]);
    printf(" k-space primitive vector (Bohr^-1)\n");
    printf("  rtv1=%10.6f %10.6f %10.6f\n", rtv[1][1], rtv[1][2], rtv[1][3]);
    printf("  rtv2=%10.6f %10.6f %10.6f\n", rtv[2][1], rtv[2][2], rtv[2][3]);
    printf("  rtv3=%10.6f %10.6f %10.6f\n\n", rtv[3][1], rtv[3][2], rtv[3][3]);
  }

  Cross_Product(tv[2], tv[3], tmp);
  CellV = Dot_Product(tv[1], tmp);
  Cell_Volume = fabs(CellV);

  if (myid == Host_ID) {
    printf("  Cell_Volume=%10.6f (Bohr^3)\n\n", Cell_Volume);
  }

  /******************************************
      find the size of the full matrix
  ******************************************/

  /* MP:
     a pointer which shows the starting number
     of basis orbitals associated with atom i
     in the full matrix */

  MP = (int *)malloc(sizeof(int) * (atomnum + 1));
  arpo = (int *)malloc(sizeof(int) * numprocs);

  fsize = 1;
  for (i = 1; i <= atomnum; i++) {
    MP[i] = fsize;
    fsize += Total_NumOrbs[i];
  }
  fsize--;

  /******************************************
               allocation arrays
  ******************************************/
  if (SpinP_switch == 0) {
    spinsize = 1;
    fsize2 = fsize;
    fsize3 = fsize + 2;
    fsize4 = Valence_Electrons / 2;
  } else if (SpinP_switch == 1) {
    spinsize = 2;
    fsize2 = fsize;
    fsize3 = fsize + 2;
    fsize4 = Valence_Electrons / 2 + abs(floor(Total_SpinS)) * 2 + 1;
  } else if (SpinP_switch == 3) {
    spinsize = 1;
    fsize2 = 2 * fsize;
    fsize3 = 2 * fsize + 2;
    fsize4 = Valence_Electrons;
  }

  /* Sopij:
     overlap matrix between one-particle wave functions
     calculated at two k-points, ki and kj */

  Sop12 = (dcomplex ***)malloc(sizeof(dcomplex **) * spinsize);
  for (spin = 0; spin < spinsize; spin++) {
    Sop12[spin] = (dcomplex **)malloc(sizeof(dcomplex *) * fsize3);
    for (i = 0; i < fsize3; i++) {
      Sop12[spin][i] = (dcomplex *)malloc(sizeof(dcomplex) * fsize3);
      for (j = 0; j < fsize3; j++) {
        Sop12[spin][i][j].r = 0.0;
        Sop12[spin][i][j].i = 0.0;
      }
    }
  }

  Sop23 = (dcomplex ***)malloc(sizeof(dcomplex **) * spinsize);
  for (spin = 0; spin < spinsize; spin++) {
    Sop23[spin] = (dcomplex **)malloc(sizeof(dcomplex *) * fsize3);
    for (i = 0; i < fsize3; i++) {
      Sop23[spin][i] = (dcomplex *)malloc(sizeof(dcomplex) * fsize3);
      for (j = 0; j < fsize3; j++) {
        Sop23[spin][i][j].r = 0.0;
        Sop23[spin][i][j].i = 0.0;
      }
    }
  }

  Sop34 = (dcomplex ***)malloc(sizeof(dcomplex **) * spinsize);
  for (spin = 0; spin < spinsize; spin++) {
    Sop34[spin] = (dcomplex **)malloc(sizeof(dcomplex *) * fsize3);
    for (i = 0; i < fsize3; i++) {
      Sop34[spin][i] = (dcomplex *)malloc(sizeof(dcomplex) * fsize3);
      for (j = 0; j < fsize3; j++) {
        Sop34[spin][i][j].r = 0.0;
        Sop34[spin][i][j].i = 0.0;
      }
    }
  }

  Sop41 = (dcomplex ***)malloc(sizeof(dcomplex **) * spinsize);
  for (spin = 0; spin < spinsize; spin++) {
    Sop41[spin] = (dcomplex **)malloc(sizeof(dcomplex *) * fsize3);
    for (i = 0; i < fsize3; i++) {
      Sop41[spin][i] = (dcomplex *)malloc(sizeof(dcomplex) * fsize3);
      for (j = 0; j < fsize3; j++) {
        Sop41[spin][i][j].r = 0.0;
        Sop41[spin][i][j].i = 0.0;
      }
    }
  }

  /* Wk1:
     wave functions at k1 */

  Wk1 = (dcomplex ***)malloc(sizeof(dcomplex **) * spinsize);
  for (spin = 0; spin < spinsize; spin++) {
    Wk1[spin] = (dcomplex **)malloc(sizeof(dcomplex *) * fsize3);
    for (i = 0; i < fsize3; i++) {
      Wk1[spin][i] = (dcomplex *)malloc(sizeof(dcomplex) * fsize3);
      for (j = 0; j < fsize3; j++) {
        Wk1[spin][i][j].r = 0.0;
        Wk1[spin][i][j].i = 0.0;
      }
    }
  }

  /* Wk2:
     wave functions at k2 */

  Wk2 = (dcomplex ***)malloc(sizeof(dcomplex **) * spinsize);
  for (spin = 0; spin < spinsize; spin++) {
    Wk2[spin] = (dcomplex **)malloc(sizeof(dcomplex *) * fsize3);
    for (i = 0; i < fsize3; i++) {
      Wk2[spin][i] = (dcomplex *)malloc(sizeof(dcomplex) * fsize3);
      for (j = 0; j < fsize3; j++) {
        Wk2[spin][i][j].r = 0.0;
        Wk2[spin][i][j].i = 0.0;
      }
    }
  }

  /* Wk3:
     wave functions at k3 */

  Wk3 = (dcomplex ***)malloc(sizeof(dcomplex **) * spinsize);
  for (spin = 0; spin < spinsize; spin++) {
    Wk3[spin] = (dcomplex **)malloc(sizeof(dcomplex *) * fsize3);
    for (i = 0; i < fsize3; i++) {
      Wk3[spin][i] = (dcomplex *)malloc(sizeof(dcomplex) * fsize3);
      for (j = 0; j < fsize3; j++) {
        Wk3[spin][i][j].r = 0.0;
        Wk3[spin][i][j].i = 0.0;
      }
    }
  }

  /* Wk4:
      wave functions at k4 */

  Wk4 = (dcomplex ***)malloc(sizeof(dcomplex **) * spinsize);
  for (spin = 0; spin < spinsize; spin++) {
    Wk4[spin] = (dcomplex **)malloc(sizeof(dcomplex *) * fsize3);
    for (i = 0; i < fsize3; i++) {
      Wk4[spin][i] = (dcomplex *)malloc(sizeof(dcomplex) * fsize3);
      for (j = 0; j < fsize3; j++) {
        Wk4[spin][i][j].r = 0.0;
        Wk4[spin][i][j].i = 0.0;
      }
    }
  }

  /* EigenVal1:
      eigenvalues at k1 */

  EigenVal1 = (double **)malloc(sizeof(double *) * spinsize);
  for (spin = 0; spin < spinsize; spin++) {
    EigenVal1[spin] = (double *)malloc(sizeof(double) * fsize3);
    for (j = 0; j < fsize3; j++)
      EigenVal1[spin][j] = 0.0;
  }

  /* EigenVal2:
     eigenvalues at k2 */

  EigenVal2 = (double **)malloc(sizeof(double *) * spinsize);
  for (spin = 0; spin < spinsize; spin++) {
    EigenVal2[spin] = (double *)malloc(sizeof(double) * fsize3);
    for (j = 0; j < fsize3; j++)
      EigenVal2[spin][j] = 0.0;
  }

  /* EigenVal3:
     eigenvalues at k3 */

  EigenVal3 = (double **)malloc(sizeof(double *) * spinsize);
  for (spin = 0; spin < spinsize; spin++) {
    EigenVal3[spin] = (double *)malloc(sizeof(double) * fsize3);
    for (j = 0; j < fsize3; j++)
      EigenVal3[spin][j] = 0.0;
  }

  /* EigenVal4:
   * eigenvalues at k4*/
  EigenVal4 = (double **)malloc(sizeof(double *) * spinsize);
  for (spin = 0; spin < spinsize; spin++) {
    EigenVal4[spin] = (double *)malloc(sizeof(double) * fsize3);
    for (j = 0; j < fsize3; j++)
      EigenVal4[spin][j] = 0.0;
  }

  /* for determinant */

  ipiv = (INTEGER *)malloc(sizeof(INTEGER) * fsize3);
  work1 = (dcomplex *)malloc(sizeof(dcomplex) * fsize3 * fsize3);
  work2 = (dcomplex *)malloc(sizeof(dcomplex) * fsize3);

  for (i = 1; i <= 3; i++) {
    Gabs[i] = sqrt(rtv[i][1] * rtv[i][1] + rtv[i][2] * rtv[i][2] +
                   rtv[i][3] * rtv[i][3]);
  }

  kloop[1][1] = 1;
  kloop[1][2] = 2;
  kloop[1][3] = 3;

  kloop[2][1] = 2;
  kloop[2][2] = 3;
  kloop[2][3] = 1;

  kloop[3][1] = 3;
  kloop[3][2] = 1;
  kloop[3][3] = 2;

  /* Computing sigmaxy in 3D system */
  for (int DirectionLoop = 1; DirectionLoop <= 3; DirectionLoop++) {
    /* only computing sigma_xy */
    if (!Dimensionflag) {
      if (DirectionLoop != 3)
        continue;
    }

    for (i = 0; i < EnergyMesh; i++) {
      sigmaxyProc[i] = 0.0;
      sigmaxySum[i] = 0.0;
    }

    if (myid == Host_ID)
      printf("Direction......... %d / 3\n", DirectionLoop);

    /* A direction of plane */

    // (1,0,0) -> b-c plane
    // (0,1,0) -> c-a plane
    // (0,0,1) -> a-b plane

    n1 = kloop[DirectionLoop][1];
    n2 = kloop[DirectionLoop][2];
    n3 = kloop[DirectionLoop][3];

    Cross_Product(rtv[n2], rtv[n3], tmp);
    Kspace = fabs(tmp[n1]);

    Nk[1] = mesh1;
    Nk[2] = mesh2;
    Nk[3] = mesh3;

    int ABNumprocs;

    int Summesh = Nk[n3] * Nk[n2];

    double ***BerryC;
    BerryC = (double ***)malloc(sizeof(double **) * EnergyMesh);
    for (k = 0; k < EnergyMesh; k++) {
      BerryC[k] = (double **)malloc(sizeof(double *) * Nk[n3]);
      for (i = 0; i < Nk[n3]; i++) {
        BerryC[k][i] = (double *)malloc(sizeof(double) * Nk[n2]);
        for (j = 0; j < Nk[n2]; j++) {
          BerryC[k][i][j] = 0.0;
        }
      }
    }

    /* for MPI */

    if (myid < Summesh % numprocs) {
      ABNumprocs = Summesh / numprocs + 1;
    } else {
      ABNumprocs = Summesh / numprocs;
    }

    MPI_Barrier(comm1);

    double **BerryCproc;

    BerryCproc = (double **)malloc(sizeof(double *) * EnergyMesh);
    for (i = 0; i < EnergyMesh; i++) {
      BerryCproc[i] = (double *)malloc(sizeof(double) * ABNumprocs);
      for (j = 0; j < ABNumprocs; j++) {
        BerryCproc[i][j] = 0.0;
      }
    }

    norm = sqrt(tv[n1][1] * tv[n1][1] + tv[n1][2] * tv[n1][2] +
                tv[n1][3] * tv[n1][3]);
    norm *= Bohr;

    kg = (double **)malloc(sizeof(double *) * 4);
    kg[0] = (double *)malloc(sizeof(double) * 1);
    for (i = 1; i <= 3; i++) {
      kg[i] = (double *)malloc(sizeof(double) * (Nk[i] + 1));
    }

    int direction;

    expOLP_HWF = (dcomplex *****)malloc(sizeof(dcomplex ****) * 4);
    for (direction = 0; direction < 4; direction++) {
      expOLP_HWF[direction] = memoryAllocation_dcomplex(expOLP_HWF[direction]);
    }

    /* make  k-point */

    for (i = 1; i <= 3; i++) {
      d = 1.0 / (double)Nk[i];
      for (j = 0; j <= Nk[i]; j++) {
        kg[i][j] = d * (double)j;
      }
    }

    double dk2, dk3, dk1;

    dk2 = 1.0 / (double)Nk[n2];
    dk3 = 1.0 / (double)Nk[n3];
    dk1 = 1.0 / (double)Nk[n1];

    /*****************************************************/
    /* Counting refinement plaquette                     */
    /*****************************************************/

    /* k3 direction's loop */
    for (i3 = 0; i3 < Nk[n1]; i3++) {
      /* Parallelization k1 k2 */
      ABloop = 0;

      if (Summesh % numprocs == 0) {
        AB_knum0 = myid * ABNumprocs;
      } else if (myid < Summesh % numprocs) {
        AB_knum0 = myid * ABNumprocs;
      } else if (myid == Summesh % numprocs) {
        AB_knum0 = (ABNumprocs + 1) * myid;
      } else {
        AB_knum0 = Summesh % numprocs + myid * ABNumprocs;
      }

      /* Parallelization */
      for (AB_knum = AB_knum0; AB_knum < AB_knum0 + ABNumprocs; AB_knum++) {

        /* Modified by N. Yamaguchi ***/
        if (n1 == 3) {
          i1 = AB_knum % mesh1;
          i2 = AB_knum / mesh1;
        } else if (n1 == 2) {
          i1 = AB_knum % mesh3;
          i2 = AB_knum / mesh3;
        } else if (n1 == 1) {
          i1 = AB_knum % mesh2;
          i2 = AB_knum / mesh2;
        }

        /* make k-point k1,k2,k3,k4 = kg[x,y,z][meshNumber] */
        k1[n2] = kg[n2][i1];
        k1[n3] = kg[n3][i2];
        k1[n1] = kg[n1][i3];

        k2[n2] = kg[n2][i1 + 1];
        k2[n3] = kg[n3][i2];
        k2[n1] = kg[n1][i3];

        k3[n2] = kg[n2][i1 + 1];
        k3[n3] = kg[n3][i2 + 1];
        k3[n1] = kg[n1][i3];

        k4[n2] = kg[n2][i1];
        k4[n3] = kg[n3][i2 + 1];
        k4[n1] = kg[n1][i3];

        /* Modified by N. Yamaguchi */
        if (n1 == 3) {
          setexpOLP(1.0 / Nk[1], 0.0, 0.0, 1, expOLP_HWF[0]);
          setexpOLP(0.0, 1.0 / Nk[2], 0.0, 1, expOLP_HWF[1]);
          setexpOLP(-1.0 / Nk[1], 0.0, 0.0, 1, expOLP_HWF[2]);
          setexpOLP(0.0, -1.0 / Nk[2], 0.0, 1, expOLP_HWF[3]);
        } else if (n1 == 2) {
          setexpOLP(0.0, 0.0, 1.0 / Nk[3], 1, expOLP_HWF[0]);
          setexpOLP(1.0 / Nk[1], 0.0, 0.0, 1, expOLP_HWF[1]);
          setexpOLP(0.0, 0.0, -1.0 / Nk[3], 1, expOLP_HWF[2]);
          setexpOLP(-1.0 / Nk[1], 0.0, 0.0, 1, expOLP_HWF[3]);
        } else {
          setexpOLP(0.0, 1.0 / Nk[2], 0.0, 1, expOLP_HWF[0]);
          setexpOLP(0.0, 0.0, 1.0 / Nk[3], 1, expOLP_HWF[1]);
          setexpOLP(0.0, -1.0 / Nk[2], 0.0, 1, expOLP_HWF[2]);
          setexpOLP(0.0, 0.0, -1.0 / Nk[3], 1, expOLP_HWF[3]);
        }

        if (myid == Host_ID)
          printf("Calculate Sop......... %d / %d\n",
                 ABloop + 1 + i3 * ABNumprocs, ABNumprocs * Nk[n1]);

        /* Added by N. Yamaguchi ***/
        EigenVal4[0][0] = Ehigh / HartreeEV + ChemP;
        /* ***/

        /* Computing S12, S23, S34, S41 */
        Calc_OverlapSopk1234(k1, k2, k3, k4, Wk1, Wk2, Wk3, Wk4, EigenVal1,
                             EigenVal2, EigenVal3, EigenVal4, Sop12, Sop23,
                             Sop34, Sop41, spinsize, fsize, fsize2, fsize3,
                             fsize4, MP);

        /* Added by N. Yamaguchi ***/
        fsize4 = MP[0];
        /* ***/

        FermiEnergy = Elow;
        for (Ecount = 0; Ecount < EnergyMesh; Ecount++) {
          Calc_BerryCEps(BerryCproc, ABloop, Ecount, FermiEnergy, ChemP, fsize2,
                         EigenVal1, EigenVal2, EigenVal3, EigenVal4, Sop12,
                         Sop23, Sop34, Sop41, work1, work2, Cdet, ipiv);
          sigmaxyProc[Ecount] += 0.5 * BerryCproc[Ecount][ABloop] / PI;
          FermiEnergy += Eunit;
        } /* Energy Loop */
        ABloop++;
      } /* AB_knum */
    }   /*  i3   */

    MPI_Barrier(comm1);

    /* Sum chernNum,Gather BerryCurvature */

    MPI_Reduce(sigmaxyProc, sigmaxySum, EnergyMesh, MPI_DOUBLE, MPI_SUM,
               Host_ID, comm1);

    /*************************************************/
    /* Output Results                                */
    /*************************************************/
    MPI_Barrier(comm1);

    for (i = 0; i < EnergyMesh; i++) {
      sigmaxySum[i] /= (double)Nk[n1];
      sigmaxySum[i] = sigmaxySum[i] * e2parh / norm * 1.0e+8; /* S / cm */
      sigmaxy[DirectionLoop][i] = sigmaxySum[i];
    }

    for (i = 0; i < EnergyMesh; i++) {
      free(BerryCproc[i]);
    }
    free(BerryCproc);

    for (i = 0; i < EnergyMesh; i++) {
      for (j = 0; j < Nk[n3]; j++) {
        free(BerryC[i][j]);
      }
      free(BerryC[i]);
    }
    free(BerryC);

    /* Added by N. Yamaguchi ***/
    for (direction = 0; direction < 4; direction++) {
      for (ct_AN = 1; ct_AN <= atomnum; ct_AN++) {
        int TNO1 = Total_NumOrbs[ct_AN];
        for (h_AN = 0; h_AN <= FNAN[ct_AN]; h_AN++) {
          for (i = 0; i < TNO1; i++) {
            free(expOLP_HWF[direction][ct_AN][h_AN][i]);
          }
          free(expOLP_HWF[direction][ct_AN][h_AN]);
        }
        free(expOLP_HWF[direction][ct_AN]);
      }
      free(expOLP_HWF[direction]);
    }
    free(expOLP_HWF);
    expOLP_HWF = NULL;
    MPI_Barrier(comm1);
  } // DirectionLoop

  double BerryCVec[4][4];
  Cross_Product(rtv[2], rtv[3], BerryCVec[1]);
  Cross_Product(rtv[3], rtv[1], BerryCVec[2]);
  Cross_Product(rtv[1], rtv[2], BerryCVec[3]);

  double normBerryCVec;
  for (j = 1; j <= 3; j++) {
    normBerryCVec = Dot_Product(BerryCVec[j], BerryCVec[j]);
    for (i = 1; i <= 3; i++)
      BerryCVec[j][i] /= sqrt(normBerryCVec);
  }

  if (myid == Host_ID) {
    printf(" primitive vector (Bohr^-1)\n");
    printf("  BerryCV1=%10.6f %10.6f %10.6f\n", BerryCVec[1][1],
           BerryCVec[1][2], BerryCVec[1][3]);
    printf("  BerryCV2=%10.6f %10.6f %10.6f\n", BerryCVec[2][1],
           BerryCVec[2][2], BerryCVec[2][3]);
    printf("  BerryCV3=%10.6f %10.6f %10.6f\n", BerryCVec[3][1],
           BerryCVec[3][2], BerryCVec[3][3]);
  }

  if (myid == Host_ID) {
    /* Output Sigmaxy */
    FILE *Cfp;
    char CdataFILE[100];
    sprintf(CdataFILE, "%s.AHC.dat", argv[1]);
    Cfp = fopen(CdataFILE, "w");
    fprintf(Cfp, "# (mesh1,mesh2,mesh3) = (%d, %d, %d)\n", mesh1, mesh2, mesh3);
    fprintf(Cfp, "# EnergyMesh: %d\n", EnergyMesh);
    fprintf(Cfp, "# Chemical Energy (eV)  sigma_yz   sigma_zx   sigma_xy\n");
    FermiEnergy = Elow;
    for (Ecount = 0; Ecount < EnergyMesh; Ecount++) {
      /* S / cm */
      double Fx, Fy, Fz;
      Fx = sigmaxy[1][Ecount] * BerryCVec[1][1] +
           sigmaxy[2][Ecount] * BerryCVec[2][1] +
           sigmaxy[3][Ecount] * BerryCVec[3][1];
      Fy = sigmaxy[1][Ecount] * BerryCVec[1][2] +
           sigmaxy[2][Ecount] * BerryCVec[2][2] +
           sigmaxy[3][Ecount] * BerryCVec[3][2];
      Fz = sigmaxy[1][Ecount] * BerryCVec[1][3] +
           sigmaxy[2][Ecount] * BerryCVec[2][3] +
           sigmaxy[3][Ecount] * BerryCVec[3][3];
      if (Unitflag) {
        double norm[4];
        for (i = 1; i <= 3; i++) {
          norm[i] = sqrt(tv[i][1] * tv[i][1] + tv[i][2] * tv[i][2] +
                         tv[i][3] * tv[i][3]);
          norm[i] *= Bohr;
        }
        fprintf(Cfp, "%10.6e %10.6e %10.6e %10.6e\n", FermiEnergy,
                Fx * norm[1] / e2parh * 1.0e-8, Fy * norm[2] / e2parh * 1.0e-8,
                Fz * norm[3] / e2parh * 1.0e-8);
      } else {
        fprintf(Cfp, "%10.6e %10.6e %10.6e %10.6e\n", FermiEnergy, Fx, Fy, Fz);
      }
      fflush(stdout);
      FermiEnergy += Eunit;
    }
    fclose(Cfp);
  } // if myid

  /******************************************
                  free arrays
  ******************************************/

  free(MP);

  free(sigmaxyProc);
  free(sigmaxySum);

  for (spin = 0; spin < spinsize; spin++) {
    for (i = 0; i < fsize3; i++) {
      free(Sop12[spin][i]);
    }
    free(Sop12[spin]);
  }
  free(Sop12);
  for (spin = 0; spin < spinsize; spin++) {
    for (i = 0; i < fsize3; i++) {
      free(Sop23[spin][i]);
    }
    free(Sop23[spin]);
  }
  free(Sop23);
  for (spin = 0; spin < spinsize; spin++) {
    for (i = 0; i < fsize3; i++) {
      free(Sop34[spin][i]);
    }
    free(Sop34[spin]);
  }
  free(Sop34);
  for (spin = 0; spin < spinsize; spin++) {
    for (i = 0; i < fsize3; i++) {
      free(Sop41[spin][i]);
    }
    free(Sop41[spin]);
  }
  free(Sop41);

  for (spin = 0; spin < spinsize; spin++) {
    for (i = 0; i < fsize3; i++) {
      free(Wk1[spin][i]);
    }
    free(Wk1[spin]);
  }
  free(Wk1);

  for (spin = 0; spin < spinsize; spin++) {
    for (i = 0; i < fsize3; i++) {
      free(Wk2[spin][i]);
    }
    free(Wk2[spin]);
  }
  free(Wk2);

  for (spin = 0; spin < spinsize; spin++) {
    for (i = 0; i < fsize3; i++) {
      free(Wk3[spin][i]);
    }
    free(Wk3[spin]);
  }
  free(Wk3);

  for (spin = 0; spin < spinsize; spin++) {
    for (i = 0; i < fsize3; i++) {
      free(Wk4[spin][i]);
    }
    free(Wk4[spin]);
  }
  free(Wk4);

  for (spin = 0; spin < spinsize; spin++) {
    free(EigenVal1[spin]);
  }
  free(EigenVal1);

  for (spin = 0; spin < spinsize; spin++) {
    free(EigenVal2[spin]);
  }
  free(EigenVal2);

  for (spin = 0; spin < spinsize; spin++) {
    free(EigenVal3[spin]);
  }
  free(EigenVal3);

  for (spin = 0; spin < spinsize; spin++) {
    free(EigenVal4[spin]);
  }
  free(EigenVal4);

  free(ipiv);
  free(work1);
  free(work2);

  for (i = 0; i <= 3; i++) {
    free(kg[i]);
  }
  free(kg);

  /* print message */

  MPI_Barrier(comm1);

  dtime(&TEtime);

  printf(" \nElapsed time = %lf (s) for myid=%3d\n", TEtime - TStime, myid);
  fflush(stdout);
  /*  */
  MPI_Barrier(comm1);

  printf(" \nThe calculation was finished normally in myid=%2d.\n", myid);
  fflush(stdout);
  /*  */

  /* MPI_Finalize */

  MPI_Finalize();

  /* return */

  return 0;
}

void Calc_OverlapSopk1234(double k1[4], double k2[4], double k3[4],
                          double k4[4], dcomplex ***Wk1, dcomplex ***Wk2,
                          dcomplex ***Wk3, dcomplex ***Wk4, double **EigenVal1,
                          double **EigenVal2, double **EigenVal3,
                          double **EigenVal4, dcomplex ***Sop12,
                          dcomplex ***Sop23, dcomplex ***Sop34,
                          dcomplex ***Sop41, int spinsize, int fsize,
                          int fsize2, int fsize3, int fsize4, int *MP) {
  /****************************************************
    calculate UUUU
       <u1|u2>  -> diag_flag = 0 diagonalize u1,u2
       <u2|u3>  -> diag_flag = 2 diagonalize only u3
       <u3|u4>  -> diag_flag = 2 diagonalize only u4
       <u4|u1>  -> diag_flag = 3 no diagonalize u4,u1
  *****************************************************/

  /* Modified by N. Yamaguchi ***/
  double upperBound = EigenVal4[0][0];
  int fsize4Modified = fsize4;
  Overlap_k1k2(0, k1, k2, spinsize, fsize, fsize2, fsize3, 0, -1, 1, MP, Sop12,
               Wk1, Wk2, EigenVal1, EigenVal2);
  Overlap_k1k2(2, k2, k3, spinsize, fsize, fsize2, fsize3, 0, -1, 1, MP, Sop23,
               Wk2, Wk3, EigenVal2, EigenVal3);
  Overlap_k1k2(2, k3, k4, spinsize, fsize, fsize2, fsize3, 0, -1, 1, MP, Sop34,
               Wk3, Wk4, EigenVal3, EigenVal4);
  int k;
  for (k = fsize4Modified; k <= fsize2; k++) {
    if (upperBound < EigenVal1[0][k]) {
      break;
    }
  }
  fsize4Modified = k;
  for (k = fsize4Modified; k <= fsize2; k++) {
    if (upperBound < EigenVal2[0][k]) {
      break;
    }
  }
  fsize4Modified = k;
  for (k = fsize4Modified; k <= fsize2; k++) {
    if (upperBound < EigenVal3[0][k]) {
      break;
    }
  }
  fsize4Modified = k;
  for (k = fsize4Modified; k <= fsize2; k++) {
    if (upperBound < EigenVal4[0][k]) {
      break;
    }
  }
  if (k > fsize2) {
    printf("The `Ehigh` is too high for this case.\n");
    MPI_Abort(MPI_COMM_WORLD, 1);
  }
  fsize4Modified = k - 1;
  expOLP = expOLP_HWF[0];
  Overlap_k1k2(3, k1, k2, spinsize, fsize, fsize2, fsize3, 1, fsize4Modified, 1,
               MP, Sop12, Wk1, Wk2, EigenVal1, EigenVal2);
  expOLP = expOLP_HWF[1];
  Overlap_k1k2(3, k2, k3, spinsize, fsize, fsize2, fsize3, 1, fsize4Modified, 1,
               MP, Sop23, Wk2, Wk3, EigenVal2, EigenVal3);
  expOLP = expOLP_HWF[2];
  Overlap_k1k2(3, k3, k4, spinsize, fsize, fsize2, fsize3, 1, fsize4Modified, 1,
               MP, Sop34, Wk3, Wk4, EigenVal3, EigenVal4);
  expOLP = expOLP_HWF[3];
  Overlap_k1k2(3, k4, k1, spinsize, fsize, fsize2, fsize3, 1, fsize4Modified, 1,
               MP, Sop41, Wk4, Wk1, EigenVal4, EigenVal1);
  MP[0] = fsize4Modified;
  /* ***/
}

void Calc_BerryCEps(double **BerryCproc, int ABloop, int Ecount,
                    double FermiEnergy, double ChemP, int fsize2,
                    double **EigenVal1, double **EigenVal2, double **EigenVal3,
                    double **EigenVal4, dcomplex ***Sop12, dcomplex ***Sop23,
                    dcomplex ***Sop34, dcomplex ***Sop41, dcomplex *work1,
                    dcomplex *work2, dcomplex Cdet[2], INTEGER *ipiv) {
  int RefinementFlag;
  int k;
  int Bandcount[5];
  int hog2;
  dcomplex UUUU;
  dcomplex U12, U23, U34, U41;

  for (k = 1; k <= 4; k++) {
    Bandcount[k] = 0;
  }

  for (k = 1; k <= fsize2; k++) {
    if (FermiEnergy > HartreeEV * (EigenVal1[0][k] - ChemP))
      Bandcount[1] += 1;
    if (FermiEnergy > HartreeEV * (EigenVal2[0][k] - ChemP))
      Bandcount[2] += 1;
    if (FermiEnergy > HartreeEV * (EigenVal3[0][k] - ChemP))
      Bandcount[3] += 1;
    if (FermiEnergy > HartreeEV * (EigenVal4[0][k] - ChemP))
      Bandcount[4] += 1;
  }

  /* Coumputing F(eps) by average method */
  if (Bandcount[1] == Bandcount[2] && Bandcount[2] == Bandcount[3] &&
      Bandcount[3] == Bandcount[4] && Bandcount[1] != 0) {
    hog2 = Bandcount[1];
    determinant(0, hog2, Sop12[0], ipiv, work1, work2, Cdet);
    U12.r = Cdet[0].r;
    U12.i = Cdet[0].i;
    determinant(0, hog2, Sop23[0], ipiv, work1, work2, Cdet);
    U23.r = Cdet[0].r;
    U23.i = Cdet[0].i;
    determinant(0, hog2, Sop34[0], ipiv, work1, work2, Cdet);
    U34.r = Cdet[0].r;
    U34.i = Cdet[0].i;
    determinant(0, hog2, Sop41[0], ipiv, work1, work2, Cdet);
    U41.r = Cdet[0].r;
    U41.i = Cdet[0].i;
    UUUU = Cmul(Cmul(U12, U23), Cmul(U34, U41));
    BerryCproc[Ecount][ABloop] = atan2(UUUU.i, UUUU.r);
  } else {
    /* Calculation Average */
    BerryCproc[Ecount][ABloop] = 0.0;
    for (k = 1; k <= 4; k++) {
      if (Bandcount[k] != 0) {
        hog2 = Bandcount[k];
        determinant(0, hog2, Sop12[0], ipiv, work1, work2, Cdet);
        U12.r = Cdet[0].r;
        U12.i = Cdet[0].i;
        determinant(0, hog2, Sop23[0], ipiv, work1, work2, Cdet);
        U23.r = Cdet[0].r;
        U23.i = Cdet[0].i;
        determinant(0, hog2, Sop34[0], ipiv, work1, work2, Cdet);
        U34.r = Cdet[0].r;
        U34.i = Cdet[0].i;
        determinant(0, hog2, Sop41[0], ipiv, work1, work2, Cdet);
        U41.r = Cdet[0].r;
        U41.i = Cdet[0].i;
        UUUU = Cmul(Cmul(U12, U23), Cmul(U34, U41));
        BerryCproc[Ecount][ABloop] += atan2(UUUU.i, UUUU.r);
      } else {
        BerryCproc[Ecount][ABloop] += 0.0;
      }
    }
    /* Average */
    BerryCproc[Ecount][ABloop] = 0.25 * BerryCproc[Ecount][ABloop];
  }
}

void Overlap_Band(double ****OLP, dcomplex **S, int *MP, double k1, double k2,
                  double k3) {
  static int i, j, wanA, wanB, tnoA, tnoB, Anum, Bnum, NUM, GA_AN, LB_AN, GB_AN;
  static int l1, l2, l3, Rn, n2;
  static double **S1, **S2;
  static double kRn, si, co, s;

  Anum = 1;
  for (i = 1; i <= atomnum; i++) {
    MP[i] = Anum;
    Anum += Total_NumOrbs[i];
  }
  NUM = Anum - 1;

  /****************************************************
                       Allocation
  ****************************************************/

  n2 = NUM + 2;

  S1 = (double **)malloc(sizeof(double *) * n2);
  for (i = 0; i < n2; i++) {
    S1[i] = (double *)malloc(sizeof(double) * n2);
  }

  S2 = (double **)malloc(sizeof(double *) * n2);
  for (i = 0; i < n2; i++) {
    S2[i] = (double *)malloc(sizeof(double) * n2);
  }

  /****************************************************
                       set overlap
  ****************************************************/

  S[0][0].r = NUM;

  for (i = 1; i <= NUM; i++) {
    for (j = 1; j <= NUM; j++) {
      S1[i][j] = 0.0;
      S2[i][j] = 0.0;
    }
  }

  for (GA_AN = 1; GA_AN <= atomnum; GA_AN++) {
    tnoA = Total_NumOrbs[GA_AN];
    Anum = MP[GA_AN];

    for (LB_AN = 0; LB_AN <= FNAN[GA_AN]; LB_AN++) {
      GB_AN = natn[GA_AN][LB_AN];
      Rn = ncn[GA_AN][LB_AN];
      tnoB = Total_NumOrbs[GB_AN];

      l1 = atv_ijk[Rn][1];
      l2 = atv_ijk[Rn][2];
      l3 = atv_ijk[Rn][3];
      kRn = k1 * (double)l1 + k2 * (double)l2 + k3 * (double)l3;

      si = sin(2.0 * PI * kRn);
      co = cos(2.0 * PI * kRn);
      Bnum = MP[GB_AN];
      for (i = 0; i < tnoA; i++) {
        for (j = 0; j < tnoB; j++) {
          s = OLP[GA_AN][LB_AN][i][j];
          S1[Anum + i][Bnum + j] += s * co;
          S2[Anum + i][Bnum + j] += s * si;
        }
      }
    }
  }

  for (i = 1; i <= NUM; i++) {
    for (j = 1; j <= NUM; j++) {
      S[i][j].r = S1[i][j];
      S[i][j].i = S2[i][j];
    }
  }

  /****************************************************
                       free arrays
  ****************************************************/

  for (i = 0; i < n2; i++) {
    free(S1[i]);
    free(S2[i]);
  }
  free(S1);
  free(S2);
}

void Hamiltonian_Band(double ****RH, dcomplex **H, int *MP, double k1,
                      double k2, double k3) {
  static int i, j, wanA, wanB, tnoA, tnoB, Anum, Bnum, NUM, GA_AN, LB_AN, GB_AN;
  static int l1, l2, l3, Rn, n2;
  static double **H1, **H2;
  static double kRn, si, co, h;

  Anum = 1;
  for (i = 1; i <= atomnum; i++) {
    MP[i] = Anum;
    Anum += Total_NumOrbs[i];
  }
  NUM = Anum - 1;

  /****************************************************
                       Allocation
  ****************************************************/

  n2 = NUM + 2;

  H1 = (double **)malloc(sizeof(double *) * n2);
  for (i = 0; i < n2; i++) {
    H1[i] = (double *)malloc(sizeof(double) * n2);
  }

  H2 = (double **)malloc(sizeof(double *) * n2);
  for (i = 0; i < n2; i++) {
    H2[i] = (double *)malloc(sizeof(double) * n2);
  }

  /****************************************************
                    set Hamiltonian
  ****************************************************/

  H[0][0].r = 2.0 * NUM;
  for (i = 1; i <= NUM; i++) {
    for (j = 1; j <= NUM; j++) {
      H1[i][j] = 0.0;
      H2[i][j] = 0.0;
    }
  }

  for (GA_AN = 1; GA_AN <= atomnum; GA_AN++) {
    tnoA = Total_NumOrbs[GA_AN];
    Anum = MP[GA_AN];

    for (LB_AN = 0; LB_AN <= FNAN[GA_AN]; LB_AN++) {
      GB_AN = natn[GA_AN][LB_AN];
      Rn = ncn[GA_AN][LB_AN];
      tnoB = Total_NumOrbs[GB_AN];

      l1 = atv_ijk[Rn][1];
      l2 = atv_ijk[Rn][2];
      l3 = atv_ijk[Rn][3];
      kRn = k1 * (double)l1 + k2 * (double)l2 + k3 * (double)l3;

      si = sin(2.0 * PI * kRn);
      co = cos(2.0 * PI * kRn);
      Bnum = MP[GB_AN];
      for (i = 0; i < tnoA; i++) {
        for (j = 0; j < tnoB; j++) {
          h = RH[GA_AN][LB_AN][i][j];
          H1[Anum + i][Bnum + j] += h * co;
          H2[Anum + i][Bnum + j] += h * si;
        }
      }
    }
  }

  for (i = 1; i <= NUM; i++) {
    for (j = 1; j <= NUM; j++) {
      H[i][j].r = H1[i][j];
      H[i][j].i = H2[i][j];
    }
  }

  /****************************************************
                       free arrays
  ****************************************************/

  for (i = 0; i < n2; i++) {
    free(H1[i]);
    free(H2[i]);
  }
  free(H1);
  free(H2);
}

void Hamiltonian_Band_NC(double *****RH, double *****IH, dcomplex **H, int *MP,
                         double k1, double k2, double k3) {
  static int i, j, k, wanA, wanB, tnoA, tnoB, Anum, Bnum;
  static int NUM, GA_AN, LB_AN, GB_AN;
  static int l1, l2, l3, Rn, n2;
  static double **H11r, **H11i;
  static double **H22r, **H22i;
  static double **H12r, **H12i;
  static double kRn, si, co, h;

  /* set MP */

  Anum = 1;
  for (i = 1; i <= atomnum; i++) {
    MP[i] = Anum;
    Anum += Total_NumOrbs[i];
  }
  NUM = Anum - 1;
  n2 = NUM + 2;

  /*******************************************
   allocation of H11r, H11i,
                 H22r, H22i,
                 H12r, H12i
  *******************************************/

  H11r = (double **)malloc(sizeof(double *) * n2);
  for (i = 0; i < n2; i++) {
    H11r[i] = (double *)malloc(sizeof(double) * n2);
  }

  H11i = (double **)malloc(sizeof(double *) * n2);
  for (i = 0; i < n2; i++) {
    H11i[i] = (double *)malloc(sizeof(double) * n2);
  }

  H22r = (double **)malloc(sizeof(double *) * n2);
  for (i = 0; i < n2; i++) {
    H22r[i] = (double *)malloc(sizeof(double) * n2);
  }

  H22i = (double **)malloc(sizeof(double *) * n2);
  for (i = 0; i < n2; i++) {
    H22i[i] = (double *)malloc(sizeof(double) * n2);
  }

  H12r = (double **)malloc(sizeof(double *) * n2);
  for (i = 0; i < n2; i++) {
    H12r[i] = (double *)malloc(sizeof(double) * n2);
  }

  H12i = (double **)malloc(sizeof(double *) * n2);
  for (i = 0; i < n2; i++) {
    H12i[i] = (double *)malloc(sizeof(double) * n2);
  }

  /****************************************************
                    set Hamiltonian
  ****************************************************/

  H[0][0].r = 2.0 * NUM;
  for (i = 1; i <= NUM; i++) {
    for (j = 1; j <= NUM; j++) {
      H11r[i][j] = 0.0;
      H11i[i][j] = 0.0;
      H22r[i][j] = 0.0;
      H22i[i][j] = 0.0;
      H12r[i][j] = 0.0;
      H12i[i][j] = 0.0;
    }
  }

  for (GA_AN = 1; GA_AN <= atomnum; GA_AN++) {
    tnoA = Total_NumOrbs[GA_AN];
    Anum = MP[GA_AN];

    for (LB_AN = 0; LB_AN <= FNAN[GA_AN]; LB_AN++) {
      GB_AN = natn[GA_AN][LB_AN];
      Rn = ncn[GA_AN][LB_AN];
      tnoB = Total_NumOrbs[GB_AN];

      l1 = atv_ijk[Rn][1];
      l2 = atv_ijk[Rn][2];
      l3 = atv_ijk[Rn][3];
      kRn = k1 * (double)l1 + k2 * (double)l2 + k3 * (double)l3;

      si = sin(2.0 * PI * kRn);
      co = cos(2.0 * PI * kRn);
      Bnum = MP[GB_AN];

      for (i = 0; i < tnoA; i++) {
        for (j = 0; j < tnoB; j++) {

          H11r[Anum + i][Bnum + j] +=
              co * RH[0][GA_AN][LB_AN][i][j] - si * IH[0][GA_AN][LB_AN][i][j];
          H11i[Anum + i][Bnum + j] +=
              si * RH[0][GA_AN][LB_AN][i][j] + co * IH[0][GA_AN][LB_AN][i][j];
          H22r[Anum + i][Bnum + j] +=
              co * RH[1][GA_AN][LB_AN][i][j] - si * IH[1][GA_AN][LB_AN][i][j];
          H22i[Anum + i][Bnum + j] +=
              si * RH[1][GA_AN][LB_AN][i][j] + co * IH[1][GA_AN][LB_AN][i][j];
          H12r[Anum + i][Bnum + j] +=
              co * RH[2][GA_AN][LB_AN][i][j] -
              si * (RH[3][GA_AN][LB_AN][i][j] + IH[2][GA_AN][LB_AN][i][j]);
          H12i[Anum + i][Bnum + j] +=
              si * RH[2][GA_AN][LB_AN][i][j] +
              co * (RH[3][GA_AN][LB_AN][i][j] + IH[2][GA_AN][LB_AN][i][j]);
        }
      }
    }
  }

  /******************************************************
    the full complex matrix of H
  ******************************************************/

  for (i = 1; i <= NUM; i++) {
    for (j = 1; j <= NUM; j++) {
      H[i][j].r = H11r[i][j];
      H[i][j].i = H11i[i][j];
      H[i + NUM][j + NUM].r = H22r[i][j];
      H[i + NUM][j + NUM].i = H22i[i][j];
      H[i][j + NUM].r = H12r[i][j];
      H[i][j + NUM].i = H12i[i][j];
      H[j + NUM][i].r = H[i][j + NUM].r;
      H[j + NUM][i].i = -H[i][j + NUM].i;
    }
  }

  /****************************************************
                       free arrays
  ****************************************************/

  for (i = 0; i < n2; i++) {
    free(H11r[i]);
  }
  free(H11r);

  for (i = 0; i < n2; i++) {
    free(H11i[i]);
  }
  free(H11i);

  for (i = 0; i < n2; i++) {
    free(H22r[i]);
  }
  free(H22r);

  for (i = 0; i < n2; i++) {
    free(H22i[i]);
  }
  free(H22i);

  for (i = 0; i < n2; i++) {
    free(H12r[i]);
  }
  free(H12r);

  for (i = 0; i < n2; i++) {
    free(H12i[i]);
  }
  free(H12i);
}

void Eigen_HH(dcomplex **ac, double *ko, int n, int EVmax) {
  /**********************************************************************
    Eigen_HH:

    Eigen_HH.c is a subroutine to solve a standard eigenvalue problem
    with a Hermite complex matrix using Householder method and lapack's
    F77_NAME(dstegr,DSTEGR)() or dstedc_().

    Log of Eigen_HH.c:

       Dec/07/2004  Released by T.Ozaki

  ***********************************************************************/

  static double ABSTOL = 1.0e-13;

  static dcomplex **ad, *u, *b1, *p, *q, tmp1, tmp2, u1, u2, p1, ss;
  static double *D, *E, *uu, *alphar, *alphai, s1, s2, s3, r, sum, ar, ai, br,
      bi, e, a1, a2, a3, a4, a5, a6, b7, r1, r2, r3, x1, x2, xap, bb, bb1, ui,
      uj, uij;

  static int jj, jj1, jj2, k, ii, ll, i3, i2, j2, i, j, i1, j1, n1, n2, ik, jk,
      po1, nn, count;

  static double Stime, Etime;
  static double Stime1, Etime1;
  static double Stime2, Etime2;
  static double time1, time2;

  /****************************************************
    allocation of arrays:
  ****************************************************/

  n2 = n + 5;

  ad = (dcomplex **)malloc(sizeof(dcomplex *) * n2);
  for (i = 0; i < n2; i++) {
    ad[i] = (dcomplex *)malloc(sizeof(dcomplex) * n2);
  }

  b1 = (dcomplex *)malloc(sizeof(dcomplex) * n2);
  u = (dcomplex *)malloc(sizeof(dcomplex) * n2);
  uu = (double *)malloc(sizeof(double) * n2);
  p = (dcomplex *)malloc(sizeof(dcomplex) * n2);
  q = (dcomplex *)malloc(sizeof(dcomplex) * n2);

  D = (double *)malloc(sizeof(double) * n2);
  E = (double *)malloc(sizeof(double) * n2);

  alphar = (double *)malloc(sizeof(double) * n2);
  alphai = (double *)malloc(sizeof(double) * n2);

  for (i = 1; i <= (n + 2); i++) {
    uu[i] = 0.0;
  }

  if (measure_time == 1)
    printf("size n=%3d EVmax=%2d\n", n, EVmax);
  if (measure_time == 1)
    dtime(&Stime);

  /****************************************************
               Householder transformation
  ****************************************************/

  for (i = 1; i <= (n - 1); i++) {

    s1 = ac[i + 1][i].r * ac[i + 1][i].r + ac[i + 1][i].i * ac[i + 1][i].i;
    s2 = 0.0;

    u[i + 1].r = ac[i + 1][i].r;
    u[i + 1].i = ac[i + 1][i].i;

    for (i1 = i + 2; i1 <= n; i1++) {

      tmp1.r = ac[i1][i].r;
      tmp1.i = ac[i1][i].i;

      s2 += tmp1.r * tmp1.r + tmp1.i * tmp1.i;

      u[i1].r = tmp1.r;
      u[i1].i = tmp1.i;
    }

    s3 = fabs(s1 + s2);

    if (ABSTOL < (fabs(ac[i + 1][i].r) + fabs(ac[i + 1][i].i))) {
      if (ac[i + 1][i].r < 0.0)
        s3 = sqrt(s3);
      else
        s3 = -sqrt(s3);
    } else {
      s3 = sqrt(s3);
    }

    if (ABSTOL < fabs(s2) || i == (n - 1)) {

      ss.r = ac[i + 1][i].r;
      ss.i = ac[i + 1][i].i;

      ac[i + 1][i].r = s3;
      ac[i + 1][i].i = 0.0;
      ac[i][i + 1].r = s3;
      ac[i][i + 1].i = 0.0;

      u[i + 1].r = u[i + 1].r - s3;
      u[i + 1].i = u[i + 1].i;

      u1.r = s3 * s3 - ss.r * s3;
      u1.i = -ss.i * s3;
      u2.r = 2.0 * u1.r;
      u2.i = 2.0 * u1.i;

      e = u2.r / (u1.r * u1.r + u1.i * u1.i);
      ar = e * u1.r;
      ai = e * u1.i;

      /* store alpha */
      alphar[i] = ar;
      alphai[i] = ai;

      /* store u2 */
      uu[i] = u2.r;

      /* store the first component of u */
      b1[i].r = ss.r - s3;
      b1[i].i = ss.i;

      r = 0.0;
      for (i1 = i + 1; i1 <= n; i1++) {

        p1.r = 0.0;
        p1.i = 0.0;
        for (j = i + 1; j <= n; j++) {
          p1.r += ac[i1][j].r * u[j].r - ac[i1][j].i * u[j].i;
          p1.i += ac[i1][j].r * u[j].i + ac[i1][j].i * u[j].r;
        }
        p[i1].r = p1.r / u1.r;
        p[i1].i = p1.i / u1.r;

        r += u[i1].r * p[i1].r + u[i1].i * p[i1].i;
      }

      r = 0.5 * r / u2.r;

      br = ar * r;
      bi = -ai * r;

      for (i1 = i + 1; i1 <= n; i1++) {
        tmp1.r = 0.5 * (p[i1].r - (br * u[i1].r - bi * u[i1].i));
        tmp1.i = 0.5 * (p[i1].i - (br * u[i1].i + bi * u[i1].r));
        q[i1].r = ar * tmp1.r - ai * tmp1.i;
        q[i1].i = ar * tmp1.i + ai * tmp1.r;
      }

      for (i1 = i + 1; i1 <= n; i1++) {
        tmp1.r = u[i1].r;
        tmp1.i = u[i1].i;
        tmp2.r = q[i1].r;
        tmp2.i = q[i1].i;
        for (j1 = i + 1; j1 <= n; j1++) {
          ac[i1][j1].r -= (tmp1.r * q[j1].r + tmp1.i * q[j1].i +
                           tmp2.r * u[j1].r + tmp2.i * u[j1].i);
          ac[i1][j1].i -= (-tmp1.r * q[j1].i + tmp1.i * q[j1].r -
                           tmp2.r * u[j1].i + tmp2.i * u[j1].r);
        }
      }
    }
  }

  for (i = 1; i <= n; i++) {
    for (j = 1; j <= n; j++) {
      ad[i][j].r = ac[i][j].r;
      ad[i][j].i = ac[i][j].i;
    }
  }

  if (measure_time == 1) {
    dtime(&Etime);
    printf("T1   %15.12f\n", Etime - Stime);
  }

  /****************************************************
                     call a lapack routine
  ****************************************************/

  if (measure_time == 1)
    dtime(&Stime);

  for (i = 1; i <= n; i++) {
    D[i - 1] = ad[i][i].r;
    E[i - 1] = ad[i][i + 1].r;
  }

  /*
  if      (dste_flag==0) lapack_dstegr2(n,D,E,ko,ac);
  else if (dste_flag==1) lapack_dstedc2(n,D,E,ko,ac);
  else if (dste_flag==2) lapack_dstevx2(n,D,E,ko,ac);
  */

  lapack_dstevx2(n, D, E, ko, ac);

  if (measure_time == 1) {
    dtime(&Etime);
    printf("T2   %15.12f\n", Etime - Stime);
  }

  /****************************************************
    transformation of eigenvectors to original space
  ****************************************************/

  if (measure_time == 1)
    dtime(&Stime);

  /* ad stores u */
  for (i = 2; i <= n; i++) {
    ad[i - 1][i].r = b1[i - 1].r;
    ad[i - 1][i].i = -b1[i - 1].i;
    ad[i][i - 1].r = b1[i - 1].r;
    ad[i][i - 1].i = b1[i - 1].i;
  }

  for (k = 1; k <= EVmax; k++) {

    for (nn = 1; nn <= n - 1; nn++) {

      if ((1.0e-3 * ABSTOL) < fabs(uu[n - nn])) {

        tmp1.r = 0.0;
        tmp1.i = 0.0;

        for (i = n - nn + 1; i <= n; i++) {
          tmp1.r += ad[n - nn][i].r * ac[k][i].r - ad[n - nn][i].i * ac[k][i].i;
          tmp1.i += ad[n - nn][i].i * ac[k][i].r + ad[n - nn][i].r * ac[k][i].i;
        }

        ss.r = (alphar[n - nn] * tmp1.r - alphai[n - nn] * tmp1.i) / uu[n - nn];
        ss.i = (alphar[n - nn] * tmp1.i + alphai[n - nn] * tmp1.r) / uu[n - nn];

        for (i = n - nn + 1; i <= n; i++) {
          ac[k][i].r -= ss.r * ad[n - nn][i].r + ss.i * ad[n - nn][i].i;
          ac[k][i].i -= -ss.r * ad[n - nn][i].i + ss.i * ad[n - nn][i].r;
        }
      }
    }
  }

  if (measure_time == 1) {
    dtime(&Etime);
    printf("T4   %15.12f\n", Etime - Stime);
  }

  /****************************************************
                     normalization
  ****************************************************/

  if (measure_time == 1)
    dtime(&Stime);

  for (j = 1; j <= EVmax; j++) {
    sum = 0.0;
    for (i = 1; i <= n; i++) {
      sum += ac[j][i].r * ac[j][i].r + ac[j][i].i * ac[j][i].i;
    }
    sum = 1.0 / sqrt(sum);
    for (i = 1; i <= n; i++) {
      ac[j][i].r = ac[j][i].r * sum;
      ac[j][i].i = ac[j][i].i * sum;
    }
  }

  if (measure_time == 1) {
    dtime(&Etime);
    printf("T5   %15.12f\n", Etime - Stime);
  }

  /****************************************************
                     transpose ac
  ****************************************************/

  for (i = 1; i <= n; i++) {
    for (j = (i + 1); j <= n; j++) {
      tmp1 = ac[i][j];
      tmp2 = ac[j][i];
      ac[i][j] = tmp2;
      ac[j][i] = tmp1;
    }
  }

  /****************************************************
                  freeing of arrays:
  ****************************************************/

  for (i = 0; i < n2; i++) {
    free(ad[i]);
  }
  free(ad);

  free(b1);
  free(u);
  free(uu);
  free(p);
  free(q);
  free(D);
  free(E);
  free(alphar);
  free(alphai);
}

void lapack_dstevx2(INTEGER N, double *D, double *E, double *W, dcomplex **ev) {
  static int i, j;

  char *JOBZ = "V";
  char *RANGE = "A";

  double VL, VU; /* dummy */
  INTEGER IL, IU;
  double ABSTOL = 1.0e-12;
  INTEGER M;
  double *Z;
  INTEGER LDZ;
  double *WORK;
  INTEGER *IWORK;
  INTEGER *IFAIL;
  INTEGER INFO;

  M = N;
  LDZ = N;

  Z = (double *)malloc(sizeof(double) * LDZ * N);
  WORK = (double *)malloc(sizeof(double) * 5 * N);
  IWORK = (INTEGER *)malloc(sizeof(INTEGER) * 5 * N);
  IFAIL = (INTEGER *)malloc(sizeof(INTEGER) * N);

  F77_NAME(dstevx, DSTEVX)
  (JOBZ, RANGE, &N, D, E, &VL, &VU, &IL, &IU, &ABSTOL, &M, W, Z, &LDZ, WORK,
   IWORK, IFAIL, &INFO);

  /* store eigenvectors */

  for (i = 0; i < N; i++) {
    for (j = 0; j < N; j++) {
      ev[i + 1][j + 1].r = Z[i * N + j];
      ev[i + 1][j + 1].i = 0.0;
    }
  }

  /* shift ko by 1 */
  for (i = N; i >= 1; i--) {
    W[i] = W[i - 1];
  }

  if (INFO > 0) {
    printf("\n error in dstevx_, info=%d\n\n", INFO);
  }
  if (INFO < 0) {
    printf("info=%d in dstevx_\n", INFO);
    exit(0);
  }

  free(Z);
  free(WORK);
  free(IWORK);
  free(IFAIL);
}

void determinant(int spin, int N, dcomplex **a, INTEGER *ipiv, dcomplex *a1d,
                 dcomplex *work, dcomplex Cdet[2]) {
  /********************************************************************
   void determinant( int spin, int N, dcomplex **a, INTEGER *ipiv,
                     dcomplex *a1d,
                     dcomplex *work, dcomplex Cdet[2] )

   a routine for calculating the determinant of overlap matrix
   for occupied states at k1- and k2-points.

   int spin (input variable)

     collinear spin-unpolarized:    spin = 0
     collinear spin-polarized:      spin = 0 or 1
     non-collinear spin-polarized:  spin = 0

   int N (input variable)

     the number of occupied states

   dcomplex **a (input variable)

     overlap matrix whose size is [fize3][fsize3].

   INTEGER *ipiv (work space)

     work space for a lapack routine, zgetrf,
     whose size is fsize3.

   dcomplex *a1d (work space)

     work space for a lapack routine, zgetrf,
     whose size is fsize3*fsize3.

   dcomplex *work (work space)

     work space for a lapack routine, zgetrf,
     whose size is fsize3.

   dcomplex Cdet[2] (output variable)

     the determinant of the matrix 'a' for occupied states
     at k1- and k2-points. Cdet[0] and Cdet[1] correspond
     to those with spin index 0 and 1.
  ********************************************************************/

  static int i, j;
  static INTEGER lda, info, lwork;
  dcomplex Ctmp;

  lda = N;
  lwork = N;

  /****************************************************
      a -> a1d
  ****************************************************/

  for (i = 0; i < N; i++) {
    for (j = 0; j < N; j++) {
      /* Modified by N. Yamaguchi ***/
      a1d[j * N + i] = a[i][j];
      /* ***/
    }
  }

  /****************************************************
                call zgetrf_() in clapack
  ****************************************************/

  F77_NAME(zgetrf, ZGETRF)(&N, &N, a1d, &lda, ipiv, &info);
  if (info != 0) {
    printf("ERROR in zgetrf_(), info=%2d\n", info);
  }

  /*
  for (i=0;i<=N;i++) {
    ipiv[i] = i + 1;
  }
  LU_fact(N,a);

  for (i=0; i<N; i++){
    a1d[i*N+i] = a[i+1][i+1];
  }
  */

  /****************************************************
               a1d -> a
  ****************************************************/

  /*
  printf("Re \n");
  for (i=0;i<N;i++) {
    for (j=0;j<N;j++) {
      printf("%15.12f ",a1d[j*N+i].r);
    }
    printf("\n");
  }

  printf("Im\n");
  for (i=0;i<N;i++) {
    for (j=0;j<N;j++) {
      printf("%15.12f ",a1d[j*N+i].i);
    }
    printf("\n");
  }
  */

  /*
  for (i=0;i<=N;i++) {
    printf("i=%2d ipiv[i]=%3d\n",i,ipiv[i]);
  }
  */

  Cdet[spin].r = 1.0;
  Cdet[spin].i = 0.0;
  for (i = 0; i < N; i++) {

    Ctmp = Cdet[spin];
    if (ipiv[i] != (i + 1)) {
      Ctmp.r = -Ctmp.r;
      Ctmp.i = -Ctmp.i;
    }

    Cdet[spin].r = Ctmp.r * a1d[i * N + i].r - Ctmp.i * a1d[i * N + i].i;
    Cdet[spin].i = Ctmp.r * a1d[i * N + i].i + Ctmp.i * a1d[i * N + i].r;
  }
}

void dtime(double *t) {
  /* real time */
  struct timeval timev;
  gettimeofday(&timev, NULL);
  *t = timev.tv_sec + (double)timev.tv_usec * 1e-6;

  /* user time + system time */
  /*
  float tarray[2];
  clock_t times(), wall;
  struct tms tbuf;
  wall = times(&tbuf);
  tarray[0] = (float) (tbuf.tms_utime / (float)CLOCKS_PER_SEC);
  tarray[1] = (float) (tbuf.tms_stime / (float)CLOCKS_PER_SEC);
  *t = (double) (tarray[0]+tarray[1]);
  printf("dtime: %lf\n",*t);
  */
}

void LU_fact(int n, dcomplex **a) {
  /*    LU factorization */
  int i, j, k;
  dcomplex w, sum;

  for (k = 1; k <= n - 1; k++) {
    w = RCdiv(1.0, a[k][k]);
    for (i = k + 1; i <= n; i++) {
      a[i][k] = Cmul(w, a[i][k]);
      for (j = k + 1; j <= n; j++) {
        a[i][j] = Csub(a[i][j], Cmul(a[i][k], a[k][j]));
      }
    }
  }
}

dcomplex RCdiv(double x, dcomplex a) {
  dcomplex c;
  double xx, yy, w;
  xx = a.r;
  yy = a.i;
  w = xx * xx + yy * yy;
  c.r = x * a.r / w;
  c.i = -x * a.i / w;
  return c;
}

dcomplex Csub(dcomplex a, dcomplex b) {
  dcomplex c;
  c.r = a.r - b.r;
  c.i = a.i - b.i;
  return c;
}

dcomplex Cmul(dcomplex a, dcomplex b) {
  dcomplex c;
  c.r = a.r * b.r - a.i * b.i;
  c.i = a.i * b.r + a.r * b.i;
  return c;
}
dcomplex Cdiv(dcomplex a, dcomplex b) {
  dcomplex c;
  c.r = (a.r * b.r + a.i * b.i) / (b.r * b.r + b.i * b.i);
  c.i = (a.i * b.r - a.r * b.i) / (b.r * b.r + b.i * b.i);
  return c;
}

void Cross_Product(double a[4], double b[4], double c[4]) {
  c[1] = a[2] * b[3] - a[3] * b[2];
  c[2] = a[3] * b[1] - a[1] * b[3];
  c[3] = a[1] * b[2] - a[2] * b[1];
}

double Dot_Product(double a[4], double b[4]) {
  static double sum;
  sum = a[1] * b[1] + a[2] * b[2] + a[3] * b[3];
  return sum;
}

/* Added by N. Yamaguchi ***/
static void expOLP_Band(dcomplex ****expOLP, dcomplex **T, int *MP, double k1,
                        double k2, double k3, double dka, double dkb,
                        double dkc, char sign) {
  int i, j, wanA, wanB, tnoA, tnoB, Anum, Bnum, NUM, GA_AN, LB_AN, GB_AN;
  int l1, l2, l3, Rn, n2;
  double **T1, **T2;
  double kRn, si, co;
  dcomplex t;
  double si2, co2;

  Anum = 1;
  for (i = 1; i <= atomnum; i++) {
    MP[i] = Anum;
    Anum += Total_NumOrbs[i];
  }
  NUM = Anum - 1;

  /****************************************************
    Allocation
   ****************************************************/

  n2 = NUM + 2;

  T1 = (double **)malloc(sizeof(double *) * n2);
  for (i = 0; i < n2; i++) {
    T1[i] = (double *)malloc(sizeof(double) * n2);
  }

  T2 = (double **)malloc(sizeof(double *) * n2);
  for (i = 0; i < n2; i++) {
    T2[i] = (double *)malloc(sizeof(double) * n2);
  }

  /****************************************************
    set T
   ****************************************************/

  T[0][0].r = NUM;

  for (i = 1; i <= NUM; i++) {
    for (j = 1; j <= NUM; j++) {
      T1[i][j] = 0.0;
      T2[i][j] = 0.0;
    }
  }

  if (sign == '-') {
    for (GA_AN = 1; GA_AN <= atomnum; GA_AN++) {
      tnoA = Total_NumOrbs[GA_AN];
      Anum = MP[GA_AN];

      for (LB_AN = 0; LB_AN <= FNAN[GA_AN]; LB_AN++) {
        GB_AN = natn[GA_AN][LB_AN];
        Rn = ncn[GA_AN][LB_AN];
        tnoB = Total_NumOrbs[GB_AN];

        l1 = atv_ijk[Rn][1];
        l2 = atv_ijk[Rn][2];
        l3 = atv_ijk[Rn][3];
        kRn = k1 * (double)l1 + k2 * (double)l2 + k3 * (double)l3;

        si = sin(2.0 * PI * kRn);
        co = cos(2.0 * PI * kRn);
        Bnum = MP[GB_AN];
        kRn -= dka * l1 + dkb * l2 + dkc * l3;
        si2 = sin(2.0 * PI * kRn);
        co2 = cos(2.0 * PI * kRn);
        for (i = 0; i < tnoA; i++) {
          for (j = 0; j < tnoB; j++) {
            t = expOLP[GA_AN][LB_AN][i][j];
            T1[Anum + i][Bnum + j] += t.r * co - t.i * si;
            T2[Anum + i][Bnum + j] += t.i * co + t.r * si;
            T1[Bnum + j][Anum + i] += t.r * co2 + t.i * si2;
            T2[Bnum + j][Anum + i] += t.i * co2 - t.r * si2;
          }
        }
      }
    }
  } else if (sign == '+') {
    for (GA_AN = 1; GA_AN <= atomnum; GA_AN++) {
      tnoA = Total_NumOrbs[GA_AN];
      Anum = MP[GA_AN];

      for (LB_AN = 0; LB_AN <= FNAN[GA_AN]; LB_AN++) {
        GB_AN = natn[GA_AN][LB_AN];
        Rn = ncn[GA_AN][LB_AN];
        tnoB = Total_NumOrbs[GB_AN];

        l1 = atv_ijk[Rn][1];
        l2 = atv_ijk[Rn][2];
        l3 = atv_ijk[Rn][3];
        kRn = k1 * (double)l1 + k2 * (double)l2 + k3 * (double)l3;

        si = sin(2.0 * PI * kRn);
        co = cos(2.0 * PI * kRn);
        Bnum = MP[GB_AN];
        kRn += dka * l1 + dkb * l2 + dkc * l3;
        si2 = sin(2.0 * PI * kRn);
        co2 = cos(2.0 * PI * kRn);
        for (i = 0; i < tnoA; i++) {
          for (j = 0; j < tnoB; j++) {
            t = expOLP[GA_AN][LB_AN][i][j];
            T1[Anum + i][Bnum + j] += t.r * co + t.i * si;
            T2[Anum + i][Bnum + j] -= t.i * co - t.r * si;
            T1[Bnum + j][Anum + i] += t.r * co2 - t.i * si2;
            T2[Bnum + j][Anum + i] -= t.i * co2 + t.r * si2;
          }
        }
      }
    }
  } else if (sign == 'M') {
    for (GA_AN = 1; GA_AN <= atomnum; GA_AN++) {
      tnoA = Total_NumOrbs[GA_AN];
      Anum = MP[GA_AN];

      for (LB_AN = 0; LB_AN <= FNAN[GA_AN]; LB_AN++) {
        GB_AN = natn[GA_AN][LB_AN];
        Rn = ncn[GA_AN][LB_AN];
        tnoB = Total_NumOrbs[GB_AN];

        l1 = atv_ijk[Rn][1];
        l2 = atv_ijk[Rn][2];
        l3 = atv_ijk[Rn][3];
        kRn = k1 * (double)l1 + k2 * (double)l2 + k3 * (double)l3;

        si = sin(2.0 * PI * kRn);
        co = cos(2.0 * PI * kRn);
        Bnum = MP[GB_AN];
        for (i = 0; i < tnoA; i++) {
          for (j = 0; j < tnoB; j++) {
            t = expOLP[GA_AN][LB_AN][i][j];
            T1[Anum + i][Bnum + j] += t.r * co - t.i * si;
            T2[Anum + i][Bnum + j] += t.i * co + t.r * si;
          }
        }
      }
    }
  } else if (sign == 'P') {
    for (GA_AN = 1; GA_AN <= atomnum; GA_AN++) {
      tnoA = Total_NumOrbs[GA_AN];
      Anum = MP[GA_AN];

      for (LB_AN = 0; LB_AN <= FNAN[GA_AN]; LB_AN++) {
        GB_AN = natn[GA_AN][LB_AN];
        Rn = ncn[GA_AN][LB_AN];
        tnoB = Total_NumOrbs[GB_AN];

        l1 = atv_ijk[Rn][1];
        l2 = atv_ijk[Rn][2];
        l3 = atv_ijk[Rn][3];
        kRn = k1 * (double)l1 + k2 * (double)l2 + k3 * (double)l3;

        si = sin(2.0 * PI * kRn);
        co = cos(2.0 * PI * kRn);
        Bnum = MP[GB_AN];
        for (i = 0; i < tnoA; i++) {
          for (j = 0; j < tnoB; j++) {
            t = expOLP[GA_AN][LB_AN][i][j];
            T1[Anum + i][Bnum + j] += t.r * co + t.i * si;
            T2[Anum + i][Bnum + j] += -t.i * co + t.r * si;
          }
        }
      }
    }
  }

  if (sign == '+' || sign == '-') {
    for (i = 1; i <= NUM; i++) {
      for (j = 1; j <= NUM; j++) {
        T[i][j].r = 0.5 * T1[i][j];
        T[i][j].i = 0.5 * T2[i][j];
      }
    }
  } else {
    for (i = 1; i <= NUM; i++) {
      for (j = 1; j <= NUM; j++) {
        T[i][j].r = T1[i][j];
        T[i][j].i = T2[i][j];
      }
    }
  }

  /****************************************************
    free arrays
   ****************************************************/

  for (i = 0; i < n2; i++) {
    free(T1[i]);
    free(T2[i]);
  }
  free(T1);
  free(T2);
}
static void Overlap_k1k2(int diag_flag, double k1[4], double k2[4],
                         int spinsize, int fsize, int fsize2, int fsize3,
                         int calcBandMin, int calcBandMax, int calcOrderMax,
                         int *MP, dcomplex ***Sop, dcomplex ***Wk1,
                         dcomplex ***Wk2, double **EigenVal1,
                         double **EigenVal2) {
  int spin;
  int ik, i1, j1, i, j, l, k;
  int ct_AN, h_AN, mu1, mu2;
  int Anum, Bnum, tnoA, tnoB;
  int mn, jj1, ii1, m;
  int Rnh, Gh_AN, l1, l2, l3;
  int recalc[2];
  double kpoints[2][4];
  double *ko, *M1;
  double sumr, sumi;
  double tmp1r, tmp1i;
  double tmp2r, tmp2i;
  double tmp3r, tmp3i;
  double si, co, kRn, k1da, k1db, k1dc;
  double si1, co1, si2, co2;
  double tmp2r1, tmp2i1;
  double tmp2r2, tmp2i2;
  double dx, dy, dz;
  double dkx, dky, dkz;
  double dka, dkb, dkc;
  double k1x, k1y, k1z;
  dcomplex **S, **H, **C;
  double OLP_eigen_cut = 1.0e-12;
  dcomplex Ctmp1, Ctmp2;

  /* Added by N. Yamaguchi ***/
  dcomplex *matrix1, *matrix2, *matrix3, opt, *work, alpha = {1.0, 0.0}, beta;
  matrix1 = (dcomplex *)malloc(sizeof(dcomplex) * fsize2 * fsize2);
  matrix2 = (dcomplex *)malloc(sizeof(dcomplex) * fsize2 * fsize2);
  matrix3 = (dcomplex *)malloc(sizeof(dcomplex) * fsize2 * fsize2);
  int itype, info, lwork;
  dcomplex zero = {0.0, 0.0};
  double *rwork;
  rwork = (double *)malloc(sizeof(double) * (3 * fsize2 - 2));
  /* ***/

  /****************************************************
    allocation of arrays:
   ****************************************************/

  ko = (double *)malloc(sizeof(double) * fsize3);
  M1 = (double *)malloc(sizeof(double) * fsize3);

  S = (dcomplex **)malloc(sizeof(dcomplex *) * fsize3);
  for (i = 0; i < fsize3; i++) {
    S[i] = (dcomplex *)malloc(sizeof(dcomplex) * fsize3);
    for (j = 0; j < fsize3; j++) {
      S[i][j].r = 0.0;
      S[i][j].i = 0.0;
    }
  }

  H = (dcomplex **)malloc(sizeof(dcomplex *) * fsize3);
  for (i = 0; i < fsize3; i++) {
    H[i] = (dcomplex *)malloc(sizeof(dcomplex) * fsize3);
    for (j = 0; j < fsize3; j++) {
      H[i][j].r = 0.0;
      H[i][j].i = 0.0;
    }
  }

  C = (dcomplex **)malloc(sizeof(dcomplex *) * fsize3);
  for (i = 0; i < fsize3; i++) {
    C[i] = (dcomplex *)malloc(sizeof(dcomplex) * fsize3);
    for (j = 0; j < fsize3; j++) {
      C[i][j].r = 0.0;
      C[i][j].i = 0.0;
    }
  }

  /* set kpoints */

  kpoints[0][1] = k1[1];
  kpoints[0][2] = k1[2];
  kpoints[0][3] = k1[3];

  kpoints[1][1] = k2[1];
  kpoints[1][2] = k2[2];
  kpoints[1][3] = k2[3];

  k1x = k1[1] * rtv[1][1] + k1[2] * rtv[2][1] + k1[3] * rtv[3][1];
  k1y = k1[1] * rtv[1][2] + k1[2] * rtv[2][2] + k1[3] * rtv[3][2];
  k1z = k1[1] * rtv[1][3] + k1[2] * rtv[2][3] + k1[3] * rtv[3][3];

  dka = k2[1] - k1[1];
  dkb = k2[2] - k1[2];
  dkc = k2[3] - k1[3];

  dkx = dka * rtv[1][1] + dkb * rtv[2][1] + dkc * rtv[3][1];
  dky = dka * rtv[1][2] + dkb * rtv[2][2] + dkc * rtv[3][2];
  dkz = dka * rtv[1][3] + dkb * rtv[2][3] + dkc * rtv[3][3];

  if (diag_flag == 0) {
    recalc[0] = 1;
    recalc[1] = 1;
  } else if (diag_flag == 1) {
    recalc[0] = 1;
    recalc[1] = 0;
  } else if (diag_flag == 2) {
    recalc[0] = 0;
    recalc[1] = 1;
  } else if (diag_flag == 3) {
    recalc[0] = 0;
    recalc[1] = 0;
  }

  /****************************************************
    diagonalize Bloch matrix at k-points, k1 and k2
   ****************************************************/

  /* Modified by N. Yamaguchi ***/
  itype = 1;
  lwork = -1;
  zhegv_(&itype, "V", "U", &fsize2, matrix1, &fsize2, matrix2, &fsize2, NULL,
         &opt, &lwork, rwork, &info);
  lwork = opt.r;
  work = (dcomplex *)malloc(sizeof(dcomplex) * lwork);
  for (ik = 0; ik < 2; ik++) {

    if (recalc[ik]) {

      Overlap_Band(OLP, S, MP, kpoints[ik][1], kpoints[ik][2], kpoints[ik][3]);
      for (spin = 0; spin < spinsize; spin++) {
        if (SpinP_switch == 3) {
          Hamiltonian_Band_NC(Hks, iHks, H, MP, kpoints[ik][1], kpoints[ik][2],
                              kpoints[ik][3]);
        } else {
          Hamiltonian_Band(Hks[spin], H, MP, kpoints[ik][1], kpoints[ik][2],
                           kpoints[ik][3]);
        }
        for (i1 = 1; i1 <= fsize2; i1++) {
          for (j1 = 1; j1 <= fsize2; j1++) {
            matrix1[(i1 - 1) + (j1 - 1) * fsize2] = H[i1][j1];
          }
        }
        if (SpinP_switch != 3) {
          for (i1 = 1; i1 <= fsize2; i1++) {
            for (j1 = 1; j1 <= fsize2; j1++) {
              matrix2[(i1 - 1) + (j1 - 1) * fsize2] = S[i1][j1];
            }
          }
        } else {
          for (i1 = 1; i1 <= fsize; i1++) {
            for (j1 = 1; j1 <= fsize; j1++) {
              matrix2[(i1 - 1) + (j1 - 1) * fsize2] = S[i1][j1];
              matrix2[(i1 - 1) + (j1 - 1 + fsize) * fsize2] = zero;
              matrix2[(i1 - 1 + fsize) + (j1 - 1) * fsize2] = zero;
              matrix2[(i1 - 1 + fsize) + (j1 - 1 + fsize) * fsize2] = S[i1][j1];
            }
          }
        }
        if (ik) {
          zhegv_(&itype, "V", "U", &fsize2, matrix1, &fsize2, matrix2, &fsize2,
                 EigenVal2[spin] + 1, work, &lwork, rwork, &info);
          for (i1 = 1; i1 <= fsize2; i1++) {
            for (j1 = 1; j1 <= fsize2; j1++) {
              Wk2[spin][j1][i1] = matrix1[(i1 - 1) + (j1 - 1) * fsize2];
            }
          }
        } else {
          zhegv_(&itype, "V", "U", &fsize2, matrix1, &fsize2, matrix2, &fsize2,
                 EigenVal1[spin] + 1, work, &lwork, rwork, &info);
          for (i1 = 1; i1 <= fsize2; i1++) {
            for (j1 = 1; j1 <= fsize2; j1++) {
              Wk1[spin][j1][i1] = matrix1[(i1 - 1) + (j1 - 1) * fsize2];
            }
          }
        }
      }
    } /* if (recalc[ik]) */
  }   /* ik */

  /* Added by N. Yamaguchi ***/
  if (calcBandMax > 0 && calcBandMin > 0 && calcBandMax >= calcBandMin) {
    /* ***/

    /****************************************************
      calculate an overlap matrix between one-particle
      wave functions  calculated at two k-points, k1 and k2
     ****************************************************/

    /* Modified by N. Yamaguchi ***/
    expOLP_Band(expOLP, H, MP, k2[1], k2[2], k2[3], dka, dkb, dkc, '-');
    for (spin = 0; spin < spinsize; spin++) {
      int numBand = calcBandMax - calcBandMin + 1;
      for (mu1 = 0; mu1 < numBand; mu1++) {
        for (j1 = 1; j1 <= fsize; j1++) {
          matrix1[(j1 - 1) + mu1 * fsize] = Wk1[spin][mu1 + calcBandMin][j1];
        }
      }
      for (i1 = 1; i1 <= fsize; i1++) {
        for (j1 = 1; j1 <= fsize; j1++) {
          matrix2[(i1 - 1) + (j1 - 1) * fsize] = H[i1][j1];
        }
      }
      beta = zero;
      zgemm_("C", "N", &numBand, &fsize, &fsize, &alpha, matrix1, &fsize,
             matrix2, &fsize, &beta, matrix3, &numBand);
      for (mu1 = 0; mu1 < numBand; mu1++) {
        for (j1 = 1; j1 <= fsize; j1++) {
          matrix1[(j1 - 1) + mu1 * fsize] = Wk2[spin][mu1 + calcBandMin][j1];
        }
      }
      zgemm_("N", "N", &numBand, &numBand, &fsize, &alpha, matrix3, &numBand,
             matrix1, &fsize, &beta, matrix2, &numBand);
      for (mu1 = 0; mu1 < numBand; mu1++) {
        for (mu2 = 0; mu2 < numBand; mu2++) {
          Sop[spin][mu1][mu2] = matrix2[mu1 + mu2 * numBand];
        }
      }
      if (SpinP_switch == 3) {
        for (mu1 = 0; mu1 < numBand; mu1++) {
          for (j1 = 1; j1 <= fsize; j1++) {
            matrix1[(j1 - 1) + mu1 * fsize] =
                Wk1[spin][mu1 + calcBandMin][j1 + fsize];
          }
        }
        for (i1 = 1; i1 <= fsize; i1++) {
          for (j1 = 1; j1 <= fsize; j1++) {
            matrix2[(i1 - 1) + (j1 - 1) * fsize] = H[i1][j1];
          }
        }
        zgemm_("C", "N", &numBand, &fsize, &fsize, &alpha, matrix1, &fsize,
               matrix2, &fsize, &beta, matrix3, &numBand);
        for (mu1 = 1; mu1 <= numBand; mu1++) {
          for (j1 = 1; j1 <= fsize; j1++) {
            matrix1[(j1 - 1) + (mu1 - 1) * fsize] = Wk2[spin][mu1][j1 + fsize];
          }
        }
        for (mu1 = 0; mu1 < numBand; mu1++) {
          for (mu2 = 0; mu2 < numBand; mu2++) {
            matrix2[mu1 + mu2 * numBand] = Sop[spin][mu1][mu2];
          }
        }
        beta.r = 1.0;
        zgemm_("N", "N", &numBand, &numBand, &fsize, &alpha, matrix3, &numBand,
               matrix1, &fsize, &beta, matrix2, &numBand);
        for (mu1 = 0; mu1 < numBand; mu1++) {
          for (mu2 = 0; mu2 < numBand; mu2++) {
            Sop[spin][mu1][mu2] = matrix2[mu1 + mu2 * numBand];
          }
        }
      }
    }
    /* ***/

    /* Added by N. Yamaguchi ***/
  }
  /* ***/

  /****************************************************
    deallocation of arrays:
    (Modified by N. Yamaguchi about this comment)
   ****************************************************/

  free(ko);
  free(M1);

  for (i = 0; i < fsize3; i++) {
    free(S[i]);
  }
  free(S);

  for (i = 0; i < fsize3; i++) {
    free(H[i]);
  }
  free(H);

  for (i = 0; i < fsize3; i++) {
    free(C[i]);
  }
  free(C);

  /* Added by N. Yamaguchi ***/
  free(matrix1);
  free(matrix2);
  free(matrix3);
  free(rwork);
  free(work);
  /* ***/
}
static void setexpOLP(double dka, double dkb, double dkc, int calcOrderMax,
                      dcomplex ****expOLP) {
  int version = 1;
  double dkx = dka * rtv[1][1] + dkb * rtv[2][1] + dkc * rtv[3][1];
  double dky = dka * rtv[1][2] + dkb * rtv[2][2] + dkc * rtv[3][2];
  double dkz = dka * rtv[1][3] + dkb * rtv[2][3] + dkc * rtv[3][3];
#ifdef HWF
  double ******OLPpo = OLPpo_HWF;
#endif
  int ct_AN;
  for (ct_AN = 1; ct_AN <= atomnum; ct_AN++) {
    double phase =
        -(dkx * Gxyz[ct_AN][1] + dky * Gxyz[ct_AN][2] + dkz * Gxyz[ct_AN][3]);
    double co = cos(phase);
    double si = sin(phase);
    int tnoA = Total_NumOrbs[ct_AN];
    int h_AN;
    for (h_AN = 0; h_AN <= FNAN[ct_AN]; h_AN++) {
      int i1;
      for (i1 = 0; i1 < tnoA; i1++) {
        int Gh_AN = natn[ct_AN][h_AN];
        int tnoB = Total_NumOrbs[Gh_AN];
        int j1;
        for (j1 = 0; j1 < tnoB; j1++) {
          double tmp3r, tmp3i;
          if (version > 0) {
            int order, fac = 1, sw = -1;
            double dkxn = 1, dkyn = 1, dkzn = 1;
            for (order = 0; order <= calcOrderMax; order++) {
              if (sw == -1) {
                tmp3r = OLP[ct_AN][h_AN][i1][j1];
                tmp3i = 0;
                fac *= -1;
                sw = 1;
                continue;
              }
              dkxn *= dkx;
              dkyn *= dky;
              dkzn *= dkz;
              fac *= order;
              if (sw == 1) {
                /* For openmx3.8 OLPpo[0,1,2][order-1] -> OLPpox, OLPpoy,
                 * OLPpoz*/
                tmp3i +=
                    +dkxn / fac * OLPpo[0][order - 1][ct_AN][h_AN][i1][j1] +
                    dkyn / fac * OLPpo[1][order - 1][ct_AN][h_AN][i1][j1] +
                    dkzn / fac * OLPpo[2][order - 1][ct_AN][h_AN][i1][j1];
                /* tmp3i += +dkxn / fac * OLPpox[ct_AN][h_AN][i1][j1] +
                         dkyn / fac * OLPpoy[ct_AN][h_AN][i1][j1] +
                         dkzn / fac * OLPpoz[ct_AN][h_AN][i1][j1];*/
                sw = 0;
              } else {
                tmp3r +=
                    +dkxn / fac * OLPpo[0][order - 1][ct_AN][h_AN][i1][j1] +
                    dkyn / fac * OLPpo[1][order - 1][ct_AN][h_AN][i1][j1] +
                    dkzn / fac * OLPpo[2][order - 1][ct_AN][h_AN][i1][j1];
                /* tmp3r+=
                            +dkxn/fac*OLPpox[ct_AN][h_AN][i1][j1]
                            +dkyn/fac*OLPpoy[ct_AN][h_AN][i1][j1]
                            +dkzn/fac*OLPpoz[ct_AN][h_AN][i1][j1]; */
                fac *= -1;
                sw = 1;
              }
            }
          }
          /* else {
            tmp3r=OLP[ct_AN][h_AN][i1][j1];
            tmp3i =-dkx*OLPpox[ct_AN][h_AN][i1][j1]
            -dky*OLPpoy[ct_AN][h_AN][i1][j1]
            -dkz*OLPpoz[ct_AN][h_AN][i1][j1];
          }*/
          expOLP[ct_AN][h_AN][i1][j1].r = co * tmp3r - si * tmp3i;
          expOLP[ct_AN][h_AN][i1][j1].i = co * tmp3i + si * tmp3r;
        }
      }
    }
  }
}
static dcomplex ****memoryAllocation_dcomplex(dcomplex ****in) {
  in = (dcomplex ****)malloc(sizeof(dcomplex ***) * (atomnum + 1));
  int ct_AN;
  in[0] = NULL;
  for (ct_AN = 1; ct_AN <= atomnum; ct_AN++) {
    int TNO1 = Total_NumOrbs[ct_AN];
    in[ct_AN] = (dcomplex ***)malloc(sizeof(dcomplex **) * (FNAN[ct_AN] + 1));
    int h_AN;
    for (h_AN = 0; h_AN <= FNAN[ct_AN]; h_AN++) {
      in[ct_AN][h_AN] = (dcomplex **)malloc(sizeof(dcomplex *) * TNO1);
      int TNO2;
      int Gh_AN = natn[ct_AN][h_AN];
      TNO2 = Total_NumOrbs[Gh_AN];
      int i;
      for (i = 0; i < TNO1; i++) {
        in[ct_AN][h_AN][i] = (dcomplex *)malloc(sizeof(dcomplex) * TNO2);
      }
    }
  }
  return in;
}
