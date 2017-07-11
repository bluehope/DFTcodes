/**********************************************************************
 RestartSCFFileDFT.c
 RestartSCFFileDFT.c provies wayto restart OpenMX DFT routine from
 AAA.scfout or AAA.scfin (future)
**********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
/*  stat section */
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
/*  end stat section */
#include "openmx_common.h"
#include "mpi.h"

#include "RestartSCFFileDFT.h"
#define MAX_LINE_SIZE 256

#define USE_MALLOC  0
void Input(char *mode, FILE *fp, double *****H, double *****iHNL,
           double ****CntOLP,
           double *****CDM);
int RestartSCFFileDFT(char *mode,int SpinP_switch, int MD_iter, double *****H,double *****iHNL,
                      double ****CntOLP,
                      double *****CDM, double *etime);

int RestartSCFFileDFT(char *mode, int SpinP_switch, int MD_iter, double *****H,double *****iHNL,
                      double ****CntOLP,
                      double *****CDM, double *etime)
{
    int numprocs, myid;
    char fname[YOUSO10];
    static FILE *fp;

    MPI_Comm_size(mpi_comm_level1, &numprocs);
    MPI_Comm_rank(mpi_comm_level1, &myid);

    //if (ML_flag) return 0;
    if (0 != myid ) return 0;

    sprintf(fname, "%s%s.scfout", filepath, filename);
    //if (myid == Host_ID) {
    MPI_Barrier(mpi_comm_level1);
    if ((fp = fopen(fname, "r")) != NULL) {
        printf(" Read the scfout file (%s)\n", fname);
        fflush(stdout);

        Input(mode,fp,H,iHNL,CntOLP,CDM);


        fclose(fp);
        printf(" Done reading the scfout file (%s)\n", fname);
    }
    else {
        printf("Failure of reading the scfout file (%s).\n", fname);
        fflush(stdout);
    }
    MPI_Barrier(mpi_comm_level1);
    //}
    return 0;
}




void Input(char *mode, FILE *fp, double *****Hks, double *****iHks,
           double ****OLP,
           double *****DM)
{
    static int Gc_AN, ct_AN, h_AN, i, j, can, Gh_AN;
    static int wan1, wan2, TNO1, TNO2, spin, Rn, num_lines;
    static int k, q_AN, Gq_AN;
    static int i_vec[20], *p_vec;
    static double d_vec[20];
    static char makeinp[100];
    static char strg[MAX_LINE_SIZE];
    FILE *fp_makeinp;
    char buf[fp_bsize];          /* setvbuf */

    static int readH = 1;
    if(0 == strcasecmp(mode,"readonlydm")) {
        readH = 0;
    }
    /****************************************************
      atomnum
      spinP_switch
    ****************************************************/

    fread(i_vec, sizeof(int), 6, fp);
    atomnum = i_vec[0];
    SpinP_switch = i_vec[1];
    Catomnum = i_vec[2];
    Latomnum = i_vec[3];
    Ratomnum = i_vec[4];
    TCpyCell = i_vec[5];

    /****************************************************
    allocation of arrays:

    double atv[TCpyCell+1][4];
    ****************************************************/
#if USE_MALLOC
    atv = (double**)malloc(sizeof(double*)*(TCpyCell + 1));
    for (Rn = 0; Rn <= TCpyCell; Rn++) {
        atv[Rn] = (double*)malloc(sizeof(double) * 4);
    }
#endif

    /****************************************************
    read atv[TCpyCell+1][4];
    ****************************************************/

    for (Rn = 0; Rn <= TCpyCell; Rn++) {
        fread(atv[Rn], sizeof(double), 4, fp);
    }

    /****************************************************
    allocation of arrays:

    int atv_ijk[TCpyCell+1][4];
    ****************************************************/
#if USE_MALLOC
    atv_ijk = (int**)malloc(sizeof(int*)*(TCpyCell + 1));
    for (Rn = 0; Rn <= TCpyCell; Rn++) {
        atv_ijk[Rn] = (int*)malloc(sizeof(int) * 4);
    }
#endif

    /****************************************************
    read atv_ijk[TCpyCell+1][4];
    ****************************************************/

    for (Rn = 0; Rn <= TCpyCell; Rn++) {
        fread(atv_ijk[Rn], sizeof(int), 4, fp);
    }

    /****************************************************
    allocation of arrays:

    int Total_NumOrbs[atomnum+1];
    int FNAN[atomnum+1];
    ****************************************************/

    Total_NumOrbs = (int*)malloc(sizeof(int)*(atomnum + 1));
#if USE_MALLOC
    FNAN = (int*)malloc(sizeof(int)*(atomnum + 1));
#endif
    /****************************************************
    the number of orbitals in each atom
    ****************************************************/

    p_vec = (int*)malloc(sizeof(int)*atomnum);
    fread(p_vec, sizeof(int), atomnum, fp);
    Total_NumOrbs[0] = 1;
    for (ct_AN = 1; ct_AN <= atomnum; ct_AN++) {
        Total_NumOrbs[ct_AN] = p_vec[ct_AN - 1];
    }
    free(p_vec);

    /****************************************************
    FNAN[]:
    the number of first nearest neighbouring atoms
    ****************************************************/

    p_vec = (int*)malloc(sizeof(int)*atomnum);
    fread(p_vec, sizeof(int), atomnum, fp);
    FNAN[0] = 0;
    for (ct_AN = 1; ct_AN <= atomnum; ct_AN++) {
        FNAN[ct_AN] = p_vec[ct_AN - 1];
    }
    free(p_vec);

    /****************************************************
    allocation of arrays:

    int natn[atomnum+1][FNAN[ct_AN]+1];
    int ncn[atomnum+1][FNAN[ct_AN]+1];
    ****************************************************/
#if USE_MALLOC
    natn = (int**)malloc(sizeof(int*)*(atomnum + 1));
    for (ct_AN = 0; ct_AN <= atomnum; ct_AN++) {
        natn[ct_AN] = (int*)malloc(sizeof(int)*(FNAN[ct_AN] + 1));
    }

    ncn = (int**)malloc(sizeof(int*)*(atomnum + 1));
    for (ct_AN = 0; ct_AN <= atomnum; ct_AN++) {
        ncn[ct_AN] = (int*)malloc(sizeof(int)*(FNAN[ct_AN] + 1));
    }
#endif
    /****************************************************
    natn[][]:
    grobal index of neighboring atoms of an atom ct_AN
    ****************************************************/

    for (ct_AN = 1; ct_AN <= atomnum; ct_AN++) {
        fread(natn[ct_AN], sizeof(int), FNAN[ct_AN] + 1, fp);
    }

    /****************************************************
    ncn[][]:
    grobal index for cell of neighboring atoms
    of an atom ct_AN
    ****************************************************/

    for (ct_AN = 1; ct_AN <= atomnum; ct_AN++) {
        fread(ncn[ct_AN], sizeof(int), FNAN[ct_AN] + 1, fp);
    }

    /****************************************************
    tv[4][4]:
    unit cell vectors in Bohr
    ****************************************************/

    fread(tv[1], sizeof(double), 4, fp);
    fread(tv[2], sizeof(double), 4, fp);
    fread(tv[3], sizeof(double), 4, fp);

    /****************************************************
    rtv[4][4]:
    unit cell vectors in Bohr
    ****************************************************/

    fread(rtv[1], sizeof(double), 4, fp);
    fread(rtv[2], sizeof(double), 4, fp);
    fread(rtv[3], sizeof(double), 4, fp);

    /****************************************************
    Gxyz[][1-3]:
    atomic coordinates in Bohr
    ****************************************************/
#if USE_MALLOC
    Gxyz = (double**)malloc(sizeof(double*)*(atomnum + 1));
    for (i = 0; i<(atomnum + 1); i++) {
        Gxyz[i] = (double*)malloc(sizeof(double) * 60);
    }
#endif
    for (ct_AN = 1; ct_AN <= atomnum; ct_AN++) {
        fread(Gxyz[ct_AN], sizeof(double), 4, fp);
    }

    /****************************************************
    allocation of arrays:

    Kohn-Sham Hamiltonian

    dooble Hks[SpinP_switch+1]
    [atomnum+1]
    [FNAN[ct_AN]+1]
    [Total_NumOrbs[ct_AN]]
    [Total_NumOrbs[h_AN]];

    Overlap matrix

    dooble OLP[atomnum+1]
    [FNAN[ct_AN]+1]
    [Total_NumOrbs[ct_AN]]
    [Total_NumOrbs[h_AN]];

    Overlap matrix with position operator x, y, z

    dooble OLPpox,y,z
    [atomnum+1]
    [FNAN[ct_AN]+1]
    [Total_NumOrbs[ct_AN]]
    [Total_NumOrbs[h_AN]];

    Density matrix

    dooble DM[SpinP_switch+1]
    [atomnum+1]
    [FNAN[ct_AN]+1]
    [Total_NumOrbs[ct_AN]]
    [Total_NumOrbs[h_AN]];
    ****************************************************/
#if USE_MALLOC
    Hks = (double*****)malloc(sizeof(double****)*(SpinP_switch + 1));
    for (spin = 0; spin <= SpinP_switch; spin++) {

        Hks[spin] = (double****)malloc(sizeof(double***)*(atomnum + 1));
        for (ct_AN = 0; ct_AN <= atomnum; ct_AN++) {
            TNO1 = Total_NumOrbs[ct_AN];
            Hks[spin][ct_AN] = (double***)malloc(sizeof(double**)*(FNAN[ct_AN] + 1));
            for (h_AN = 0; h_AN <= FNAN[ct_AN]; h_AN++) {
                Hks[spin][ct_AN][h_AN] = (double**)malloc(sizeof(double*)*TNO1);

                if (ct_AN == 0) {
                    TNO2 = 1;
                }
                else {
                    Gh_AN = natn[ct_AN][h_AN];
                    TNO2 = Total_NumOrbs[Gh_AN];
                }
                for (i = 0; i<TNO1; i++) {
                    Hks[spin][ct_AN][h_AN][i] = (double*)malloc(sizeof(double)*TNO2);
                }
            }
        }
    }


    iHks = (double*****)malloc(sizeof(double****) * 3);
    for (spin = 0; spin<3; spin++) {

        iHks[spin] = (double****)malloc(sizeof(double***)*(atomnum + 1));
        for (ct_AN = 0; ct_AN <= atomnum; ct_AN++) {
            TNO1 = Total_NumOrbs[ct_AN];
            iHks[spin][ct_AN] = (double***)malloc(sizeof(double**)*(FNAN[ct_AN] + 1));
            for (h_AN = 0; h_AN <= FNAN[ct_AN]; h_AN++) {
                iHks[spin][ct_AN][h_AN] = (double**)malloc(sizeof(double*)*TNO1);

                if (ct_AN == 0) {
                    TNO2 = 1;
                }
                else {
                    Gh_AN = natn[ct_AN][h_AN];
                    TNO2 = Total_NumOrbs[Gh_AN];
                }
                for (i = 0; i<TNO1; i++) {
                    iHks[spin][ct_AN][h_AN][i] = (double*)malloc(sizeof(double)*TNO2);
                    for (j = 0; j<TNO2; j++) iHks[spin][ct_AN][h_AN][i][j] = 0.0;
                }
            }
        }
    }

    scfOLP = (double****)malloc(sizeof(double***)*(atomnum + 1));
    for (ct_AN = 0; ct_AN <= atomnum; ct_AN++) {
        TNO1 = Total_NumOrbs[ct_AN];
        scfOLP[ct_AN] = (double***)malloc(sizeof(double**)*(FNAN[ct_AN] + 1));
        for (h_AN = 0; h_AN <= FNAN[ct_AN]; h_AN++) {
            scfOLP[ct_AN][h_AN] = (double**)malloc(sizeof(double*)*TNO1);

            if (ct_AN == 0) {
                TNO2 = 1;
            }
            else {
                Gh_AN = natn[ct_AN][h_AN];
                TNO2 = Total_NumOrbs[Gh_AN];
            }
            for (i = 0; i<TNO1; i++) {
                scfOLP[ct_AN][h_AN][i] = (double*)malloc(sizeof(double)*TNO2);
            }
        }
    }
#endif


    OLPpox = (double****)malloc(sizeof(double***)*(atomnum + 1));
    for (ct_AN = 0; ct_AN <= atomnum; ct_AN++) {
        TNO1 = Total_NumOrbs[ct_AN];
        OLPpox[ct_AN] = (double***)malloc(sizeof(double**)*(FNAN[ct_AN] + 1));
        for (h_AN = 0; h_AN <= FNAN[ct_AN]; h_AN++) {
            OLPpox[ct_AN][h_AN] = (double**)malloc(sizeof(double*)*TNO1);

            if (ct_AN == 0) {
                TNO2 = 1;
            }
            else {
                Gh_AN = natn[ct_AN][h_AN];
                TNO2 = Total_NumOrbs[Gh_AN];
            }
            for (i = 0; i<TNO1; i++) {
                OLPpox[ct_AN][h_AN][i] = (double*)malloc(sizeof(double)*TNO2);
            }
        }
    }

    OLPpoy = (double****)malloc(sizeof(double***)*(atomnum + 1));
    for (ct_AN = 0; ct_AN <= atomnum; ct_AN++) {
        TNO1 = Total_NumOrbs[ct_AN];
        OLPpoy[ct_AN] = (double***)malloc(sizeof(double**)*(FNAN[ct_AN] + 1));
        for (h_AN = 0; h_AN <= FNAN[ct_AN]; h_AN++) {
            OLPpoy[ct_AN][h_AN] = (double**)malloc(sizeof(double*)*TNO1);

            if (ct_AN == 0) {
                TNO2 = 1;
            }
            else {
                Gh_AN = natn[ct_AN][h_AN];
                TNO2 = Total_NumOrbs[Gh_AN];
            }
            for (i = 0; i<TNO1; i++) {
                OLPpoy[ct_AN][h_AN][i] = (double*)malloc(sizeof(double)*TNO2);
            }
        }
    }

    OLPpoz = (double****)malloc(sizeof(double***)*(atomnum + 1));
    for (ct_AN = 0; ct_AN <= atomnum; ct_AN++) {
        TNO1 = Total_NumOrbs[ct_AN];
        OLPpoz[ct_AN] = (double***)malloc(sizeof(double**)*(FNAN[ct_AN] + 1));
        for (h_AN = 0; h_AN <= FNAN[ct_AN]; h_AN++) {
            OLPpoz[ct_AN][h_AN] = (double**)malloc(sizeof(double*)*TNO1);

            if (ct_AN == 0) {
                TNO2 = 1;
            }
            else {
                Gh_AN = natn[ct_AN][h_AN];
                TNO2 = Total_NumOrbs[Gh_AN];
            }
            for (i = 0; i<TNO1; i++) {
                OLPpoz[ct_AN][h_AN][i] = (double*)malloc(sizeof(double)*TNO2);
            }
        }
    }
#if USE_MALLOC
    CDM = (double*****)malloc(sizeof(double****)*(SpinP_switch + 1));
    for (spin = 0; spin <= SpinP_switch; spin++) {

        CDM[spin] = (double****)malloc(sizeof(double***)*(atomnum + 1));
        for (ct_AN = 0; ct_AN <= atomnum; ct_AN++) {
            TNO1 = Total_NumOrbs[ct_AN];
            CDM[spin][ct_AN] = (double***)malloc(sizeof(double**)*(FNAN[ct_AN] + 1));
            for (h_AN = 0; h_AN <= FNAN[ct_AN]; h_AN++) {
                CDM[spin][ct_AN][h_AN] = (double**)malloc(sizeof(double*)*TNO1);

                if (ct_AN == 0) {
                    TNO2 = 1;
                }
                else {
                    Gh_AN = natn[ct_AN][h_AN];
                    TNO2 = Total_NumOrbs[Gh_AN];
                }
                for (i = 0; i<TNO1; i++) {
                    CDM[spin][ct_AN][h_AN][i] = (double*)malloc(sizeof(double)*TNO2);
                }
            }
        }
    }
#endif
    /****************************************************
    Hamiltonian matrix
    ****************************************************/

    for (spin = 0; spin <= SpinP_switch; spin++) {
        for (ct_AN = 1; ct_AN <= atomnum; ct_AN++) {
            TNO1 = Total_NumOrbs[ct_AN];
            for (h_AN = 0; h_AN <= FNAN[ct_AN]; h_AN++) {
                Gh_AN = natn[ct_AN][h_AN];
                TNO2 = Total_NumOrbs[Gh_AN];
                // readH
                if(1 == readH) {
                    for (i = 0; i < TNO1; i++)
                    {
                        fread(Hks[spin][ct_AN][h_AN][i], sizeof(double), TNO2, fp);
#if 0
                        for (int j = 0; j < TNO2; j++)
                            Hks[spin][ct_AN][h_AN][i][j] = 0.0; //test
#endif
                    }
                } else {
                    /* do not update Hks */
                    double* tmpread = (double*)malloc(sizeof(double)*TNO2);
                    for (i = 0; i < TNO1; i++)
                        fread(tmpread, sizeof(double), TNO2, fp);
                    free(tmpread);
                }
            }
        }
    }

    /****************************************************
    iHks:
    imaginary Kohn-Sham matrix elements of basis orbitals
    for alpha-alpha, beta-beta, and alpha-beta spin matrices
    of which contributions come from spin-orbit coupling
    and Hubbard U effective potential.
    ****************************************************/

    if (SpinP_switch == 3) {
        for (spin = 0; spin<3; spin++) {
            for (ct_AN = 1; ct_AN <= atomnum; ct_AN++) {
                TNO1 = Total_NumOrbs[ct_AN];
                for (h_AN = 0; h_AN <= FNAN[ct_AN]; h_AN++) {
                    Gh_AN = natn[ct_AN][h_AN];
                    TNO2 = Total_NumOrbs[Gh_AN];
                    for (i = 0; i<TNO1; i++) {
                        fread(iHks[spin][ct_AN][h_AN][i], sizeof(double), TNO2, fp);
                    }
                }
            }
        }
    }

    /****************************************************
    Overlap matrix
    ****************************************************/

    for (ct_AN = 1; ct_AN <= atomnum; ct_AN++) {
        TNO1 = Total_NumOrbs[ct_AN];
        for (h_AN = 0; h_AN <= FNAN[ct_AN]; h_AN++) {
            Gh_AN = natn[ct_AN][h_AN];
            TNO2 = Total_NumOrbs[Gh_AN];
            for (i = 0; i<TNO1; i++) {
                fread(OLP[ct_AN][h_AN][i], sizeof(double), TNO2, fp);
            }
        }
    }

    /****************************************************
    Overlap matrix with position operator x
    ****************************************************/

    for (ct_AN = 1; ct_AN <= atomnum; ct_AN++) {
        TNO1 = Total_NumOrbs[ct_AN];
        for (h_AN = 0; h_AN <= FNAN[ct_AN]; h_AN++) {
            Gh_AN = natn[ct_AN][h_AN];
            TNO2 = Total_NumOrbs[Gh_AN];
            for (i = 0; i<TNO1; i++) {
                fread(OLPpox[ct_AN][h_AN][i], sizeof(double), TNO2, fp);
            }
        }
    }

    /****************************************************
    Overlap matrix with position operator y
    ****************************************************/

    for (ct_AN = 1; ct_AN <= atomnum; ct_AN++) {
        TNO1 = Total_NumOrbs[ct_AN];
        for (h_AN = 0; h_AN <= FNAN[ct_AN]; h_AN++) {
            Gh_AN = natn[ct_AN][h_AN];
            TNO2 = Total_NumOrbs[Gh_AN];
            for (i = 0; i<TNO1; i++) {
                fread(OLPpoy[ct_AN][h_AN][i], sizeof(double), TNO2, fp);
            }
        }
    }

    /****************************************************
    Overlap matrix with position operator z
    ****************************************************/

    for (ct_AN = 1; ct_AN <= atomnum; ct_AN++) {
        TNO1 = Total_NumOrbs[ct_AN];
        for (h_AN = 0; h_AN <= FNAN[ct_AN]; h_AN++) {
            Gh_AN = natn[ct_AN][h_AN];
            TNO2 = Total_NumOrbs[Gh_AN];
            for (i = 0; i<TNO1; i++) {
                fread(OLPpoz[ct_AN][h_AN][i], sizeof(double), TNO2, fp);
            }
        }
    }

    /****************************************************
    Density matrix
    ****************************************************/

    for (spin = 0; spin <= SpinP_switch; spin++) {
        for (ct_AN = 1; ct_AN <= atomnum; ct_AN++) {
            TNO1 = Total_NumOrbs[ct_AN];
            for (h_AN = 0; h_AN <= FNAN[ct_AN]; h_AN++) {
                Gh_AN = natn[ct_AN][h_AN];
                TNO2 = Total_NumOrbs[Gh_AN];
                for (i = 0; i<TNO1; i++) {
                    fread(DM[spin][ct_AN][h_AN][i], sizeof(double), TNO2, fp);
                }
            }
        }
    }

    /****************************************************
    Solver
    ****************************************************/

    fread(i_vec, sizeof(int), 1, fp);
    Solver = i_vec[0];

    /****************************************************
    ChemP
    Temp
    ****************************************************/

    fread(d_vec, sizeof(double), 10, fp);
    ChemP = d_vec[0];
    E_Temp = d_vec[1];
    dipole_moment_core[1] = d_vec[2];
    dipole_moment_core[2] = d_vec[3];
    dipole_moment_core[3] = d_vec[4];
    dipole_moment_background[1] = d_vec[5];
    dipole_moment_background[2] = d_vec[6];
    dipole_moment_background[3] = d_vec[7];
    Valence_Electrons = d_vec[8];
    Total_SpinS = d_vec[9];

    /****************************************************
    input file
    ****************************************************/

    fread(i_vec, sizeof(int), 1, fp);
    num_lines = i_vec[0];

    sprintf(makeinp, "temporal_123456.input");

    if ((fp_makeinp = fopen(makeinp, "w")) != NULL) {

#ifdef xt3
        setvbuf(fp_makeinp, buf, _IOFBF, fp_bsize);  /* setvbuf */
#endif

        for (i = 1; i <= num_lines; i++) {
            fread(strg, sizeof(char), MAX_LINE_SIZE, fp);
            fprintf(fp_makeinp, "%s", strg);
        }

        fclose(fp_makeinp);
    }
    else {
        printf("error in making temporal_12345.input\n");
    }

    /* Free up */
    for (ct_AN = 0; ct_AN <= atomnum; ct_AN++) {
        TNO1 = Total_NumOrbs[ct_AN];
        for (h_AN = 0; h_AN <= FNAN[ct_AN]; h_AN++) {

            for (i = 0; i<TNO1; i++) {
                free(OLPpox[ct_AN][h_AN][i]);
                free(OLPpoy[ct_AN][h_AN][i]);
                free(OLPpoz[ct_AN][h_AN][i]);
            }
            free(OLPpox[ct_AN][h_AN]);
            free(OLPpoy[ct_AN][h_AN]);
            free(OLPpoz[ct_AN][h_AN]);
        }
        free(OLPpox[ct_AN]);
        free(OLPpoy[ct_AN]);
        free(OLPpoz[ct_AN]);
    }
    free(OLPpox);
    free(OLPpoy);
    free(OLPpoz);

}



