/* Added by N. Yamaguchi ***/
#if defined SOField && !defined HWC || defined HWF
#define TEMP
#endif
#ifdef DEBUG_HWF_20190505
typedef float     Type_Orbs_Grid;       /* type of Orbs_Grid */
#define MPI_Type_Orbs_Grid  MPI_FLOAT   /* type of Orbs_Grid */
#else
typedef double    Type_Orbs_Grid;       /* type of Orbs_Grid */
#define MPI_Type_Orbs_Grid  MPI_DOUBLE  /* type of Orbs_Grid */
#endif
/* ***/

/*******************************************************
  int atomnun;
  the number of total atoms
 *******************************************************/
int atomnum;   

/*******************************************************
  int Catomnun;
  the number of atoms in the central region
 *******************************************************/
int Catomnum;

/*******************************************************
  int Latomnun;
  the number of atoms in the left lead
 *******************************************************/
int Latomnum;

/*******************************************************
  int Ratomnun;
  the number of atoms in the left lead
 *******************************************************/
int Ratomnum;   

/*******************************************************
  int SpinP_switch;
0: non-spin polarized 
1: spin polarized
 *******************************************************/
int SpinP_switch;

/*******************************************************
  int TCpyCell;
  the total number of periodic cells
 *******************************************************/
int TCpyCell;

/*******************************************************
  int Solver;
  method for solving eigenvalue problem
 *******************************************************/
int Solver;

/*******************************************************
  double ChemP;
  chemical potential
 *******************************************************/
double ChemP;

/*******************************************************
  int Valence_Electrons;
  total number of valence electrons
 *******************************************************/
int Valence_Electrons;

/*******************************************************
  double Total_SpinS;
  total value of Spin (2*Total_SpinS = muB)
 *******************************************************/
double Total_SpinS;

/* Added by N. Yamaguchi ***/
#ifndef TEMP
/* ***/

/*******************************************************
  double E_Temp;
  electronic temperature
 *******************************************************/
double E_Temp;

/* Added by N. Yamaguchi ***/
#else
/*******************************************************
  double E_Temp_HWF;
  electronic temperature
 *******************************************************/
double E_Temp_HWF;
#endif
/* ***/

/*******************************************************
  int *Total_NumOrbs; 
  the number of atomic orbitals in each atom
size: Total_NumOrbs[atomnum+1]
 *******************************************************/
int *Total_NumOrbs;

/*******************************************************
  int *FNAN; 
  the number of first neighboring atoms of each atom
size: FNAN[atomnum+1]
 *******************************************************/
int *FNAN;

/* Added by N. Yamaguchi ***/
#ifdef TEMP
/*******************************************************
  int *FNAN_HWF; 
  the number of first neighboring atoms of each atom
size: FNAN_HWF[atomnum+1]
 *******************************************************/
int *FNAN_HWF;
#endif
/* ***/

/*******************************************************
  int **natn; 
  grobal index of neighboring atoms of an atom ct_AN
size: natn[atomnum+1][FNAN[ct_AN]+1]
 *******************************************************/
int **natn;

/* Added by N. Yamaguchi ***/
#ifdef TEMP
/*******************************************************
  int **natn_HWF; 
  grobal index of neighboring atoms of an atom ct_AN
size: natn_HWF[atomnum+1][FNAN[ct_AN]+1]
 *******************************************************/
int **natn_HWF;
#endif
/* ***/

/* Added by N. Yamaguchi ***/
#ifndef TEMP
/* ***/

/*******************************************************
  int **ncn; 
  grobal index for cell of neighboring atoms of
  an atom ct_AN
size: ncn[atomnum+1][FNAN[ct_AN]+1]
 *******************************************************/
int **ncn;

/* Added by N. Yamaguchi ***/
#else
/*******************************************************
  int **ncn_HWF; 
  grobal index for cell of neighboring atoms of
  an atom ct_AN
size: ncn_HWF[atomnum+1][FNAN[ct_AN]+1]
 *******************************************************/
int **ncn_HWF;
#endif
/* ***/

/* Added by N. Yamaguchi ***/
#ifndef TEMP
/* ***/

/*******************************************************
  double **atv;
  x,y,and z-components of translation vector of  
  periodically copied cells
size: atv[TCpyCell+1][4];
 *******************************************************/
double **atv;

/* Added by N. Yamaguchi ***/
#else
/*******************************************************
  double **atv_HWF;
  x,y,and z-components of translation vector of  
  periodically copied cells
size: atv_HWF[TCpyCell+1][4];
 *******************************************************/
double **atv_HWF;
#endif
/* ***/

/* Added by N. Yamaguchi ***/
#ifndef TEMP
/* ***/

/*******************************************************
  int **atv_ijk;
  i,j,and k number of periodically copied cells
size: atv_ijk[TCpyCell+1][4];
 *******************************************************/
int **atv_ijk;

/* Added by N. Yamaguchi ***/
#else
/*******************************************************
  int **atv_ijk_HWF;
  i,j,and k number of periodically copied cells
size: atv_ijk_HWF[TCpyCell+1][4];
 *******************************************************/
int **atv_ijk_HWF;
#endif
/* ***/

/*******************************************************
  double tv[4][4];
  unit cell vectors in Bohr
 *******************************************************/
double tv[4][4];

/*******************************************************
  double rtv[4][4]:
  reciprocal unit cell vectors in Bohr^{-1}

note:
tv_i \dot rtv_j = 2PI * Kronecker's delta_{ij}
 *******************************************************/
double rtv[4][4];

/*******************************************************
  double Gxyz[atomnum+1][60];
  atomic coordinates in Bohr
 *******************************************************/
double **Gxyz;

/*******************************************************
  double *****Hks;
  Kohn-Sham matrix elements of basis orbitals
size: Hks[SpinP_switch+1]
[atomnum+1]
[FNAN[ct_AN]+1]
[Total_NumOrbs[ct_AN]]
[Total_NumOrbs[h_AN]] 
 *******************************************************/
double *****Hks;

/*******************************************************
  double *****iHks;
  imaginary Kohn-Sham matrix elements of basis orbitals
  for alpha-alpha, beta-beta, and alpha-beta spin matrices
  of which contributions come from spin-orbit coupling 
  and Hubbard U effective potential.
size: iHks[3]
[atomnum+1]
[FNAN[ct_AN]+1]
[Total_NumOrbs[ct_AN]]
[Total_NumOrbs[h_AN]] 
 *******************************************************/
double *****iHks;

/* Added by N. Yamaguchi ***/
#ifndef TEMP
/* ***/

/*******************************************************
  double ****OLP;
  overlap matrix
size: OLP[atomnum+1]
[FNAN[ct_AN]+1]
[Total_NumOrbs[ct_AN]]
[Total_NumOrbs[h_AN]]
 *******************************************************/
double ****OLP;

/* Added by N. Yamaguchi ***/
#else
/*******************************************************
  double ****OLP_HWF;
  overlap matrix
size: OLP_HWF[atomnum+1]
[FNAN[ct_AN]+1]
[Total_NumOrbs[ct_AN]]
[Total_NumOrbs[h_AN]]
 *******************************************************/
double ****OLP_HWF;
#endif
/* ***/

/*******************************************************
  double ****OLPpox;
  overlap matrix with position operator x
size: OLPpox[atomnum+1]
[FNAN[ct_AN]+1]
[Total_NumOrbs[ct_AN]]
[Total_NumOrbs[h_AN]] 
 *******************************************************/
double ****OLPpox;

/*******************************************************
  double ****OLPpoy;
  overlap matrix with position operator y
size: OLPpoy[atomnum+1]
[FNAN[ct_AN]+1]
[Total_NumOrbs[ct_AN]]
[Total_NumOrbs[h_AN]] 
 *******************************************************/
double ****OLPpoy;

/*******************************************************
  double ****OLPpoz;
  overlap matrix with position operator z
size: OLPpoz[atomnum+1]
[FNAN[ct_AN]+1]
[Total_NumOrbs[ct_AN]]
[Total_NumOrbs[h_AN]] 
 *******************************************************/
double ****OLPpoz;

/* Added by N. Yamaguchi ***/
#ifndef TEMP
/* ***/

/*******************************************************
  double *****DM;
  overlap matrix
size: DM[SpinP_switch+1]
[atomnum+1]
[FNAN[ct_AN]+1]
[Total_NumOrbs[ct_AN]]
[Total_NumOrbs[h_AN]]
 *******************************************************/
double *****DM;

/* Added by N. Yamaguchi ***/
#else
/*******************************************************
  double *****DM_HWF;
  overlap matrix
size: DM_HWF[SpinP_switch+1]
[atomnum+1]
[FNAN[ct_AN]+1]
[Total_NumOrbs[ct_AN]]
[Total_NumOrbs[h_AN]]
 *******************************************************/
double *****DM_HWF;
#endif

#ifndef TEMP
/*******************************************************
  int HWF_HWF; (added by N. Yamaguchi for HWC)
  int *atomicNumber;
  double Grid_Origin[4];
  double gtv[4][4];
  int Ngrid1;
  int Ngrid2;
  int Ngrid3;
  int TNumGrid;
  int *GridN_Atom;
  int **GridListAtom;
  int **CellListAtom;
  Type_Orbs_Grid ***Orbs_Grid;
 *******************************************************/
int HWF_HWF;
int *atomicNumber;
double Grid_Origin[4];
double gtv[4][4];
int Ngrid1;
int Ngrid2;
int Ngrid3;
int TNumGrid;
int *GridN_Atom;
int **GridListAtom;
int **CellListAtom;
Type_Orbs_Grid ***Orbs_Grid;
#else
/*******************************************************
  int HWF_HWF; (added by N. Yamaguchi for HWC)
  int *atomicNumber;
  double Grid_Origin_HWF[4];
  double gtv_HWF[4][4];
  int Ngrid1_HWF;
  int Ngrid2_HWF;
  int Ngrid3_HWF;
  int TNumGrid_HWF;
  int *GridN_Atom_HWF;
  int **GridListAtom_HWF;
  int **CellListAtom_HWF;
  Type_Orbs_Grid ***Orbs_Grid_HWF;
 *******************************************************/
int HWF_HWF;
int *atomicNumber;
double Grid_Origin_HWF[4];
double gtv_HWF[4][4];
int Ngrid1_HWF;
int Ngrid2_HWF;
int Ngrid3_HWF;
int TNumGrid_HWF;
int *GridN_Atom_HWF;
int **GridListAtom_HWF;
int **CellListAtom_HWF;
Type_Orbs_Grid ***Orbs_Grid_HWF;
#endif
/* ***/

/*******************************************************
  double dipole_moment_core[4];
 *******************************************************/
double dipole_moment_core[4];

/*******************************************************
  int version; (added by N. Yamaguchi)
 *******************************************************/
int version;

/*******************************************************
  int order_max; (added by N. Yamaguchi for HWC)
 *******************************************************/
int order_max;

/*******************************************************
  double *cc_vec; (added by N. Yamaguchi for HWC)
 *******************************************************/
double *cc_vec;

#if !defined HWF || defined DEBUG_HWF_A
/*******************************************************
  double ******OLPpo; (added by N. Yamaguchi for HWC)
 *******************************************************/
double ******OLPpo;
#else
/*******************************************************
  double ******OLPpo_HWF; (added by N. Yamaguchi for HWC)
 *******************************************************/
double ******OLPpo_HWF;
#endif

/*******************************************************
  double dipole_moment_background[4];
 *******************************************************/
double dipole_moment_background[4];

/* Added by N. Yamaguchi ***/
#ifndef SOField
/* ***/

void read_scfout(char *argv[]);

/* Added by N. Yamaguchi ***/
#else
void read_scfout(char *filename_wf, int Print_datFile);
void free_scfout();
#endif
/* ***/

