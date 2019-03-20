/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015, by the GROMACS development team, led by
 * Mark Abraham, David van der Spoel, Berk Hess, and Erik Lindahl,
 * and including many others, as listed in the AUTHORS file in the
 * top-level source directory and at http://www.gromacs.org.
 *
 * GROMACS is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 *
 * GROMACS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with GROMACS; if not, see
 * http://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at http://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out http://www.gromacs.org.
 */

/* 
 * This functions deal with the LIO-GROMACS interface for QM-MM calculations.
 * Main things involved are unit conversions from GROMACS MD-units to LIO 
 * Atomic Units (A.U.) and vice-versa. In addition, 2D arrays need to be
 * reordered since LIO is coded in FORTRAN, thus allocating arrays differently.
 */

#include "gmxpre.h"
#include "config.h"

#if GMX_QMMM_LIO

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "gromacs/fileio/confio.h"
#include "gromacs/gmxlib/network.h"
#include "gromacs/gmxlib/nrnb.h"
#include "gromacs/math/units.h"
#include "gromacs/math/vec.h"
#include "gromacs/mdlib/force.h"
#include "gromacs/mdlib/ns.h"
#include "gromacs/mdlib/qmmm.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/smalloc.h"

extern "C" 
{
    void 
    F77_FUNC(init_lio_gromacs, INIT_LIO_GROMACS) (int *nqmat, int  Zqmat[], 
                                                  int *nmmat, int *charge);
    void
    F77_FUNC(scf_gro, SCF_GRO) (double *energy  , double qmcoord[] , 
                                double mmcoord[], double mmcharge[], 
                                int    *nmmat  );
    void   
    F77_FUNC(dft_get_qm_forces, DFT_GET_QM_FORCES) (double qmgrad[]);
 
    void
    F77_FUNC(dft_get_mm_forces, DFT_GET_MM_FORCES) (double mmgrad[] ,
                                                    double qmgrad[]);
}


void init_lio(t_QMrec *qm, t_MMrec *mm)
{
/* 
* Calls Lio initialization subroutines. This does not need any kind of unit   
*   conversion or the like, since arguments passed are integers and Lio options
*   are read from an additional input file in the work directory.             
*/
    F77_FUNC(init_lio_gromacs, INIT_LIO_GROMACS) (&qm->nrQMatoms, 
             qm->atomicnumberQM, &mm->nrMMatoms, &qm->QMcharge );
}; /* init_lio */


real call_lio(const t_forcerec *fr, t_QMrec *qm, t_MMrec *mm, rvec f[], rvec fshift[])
{
/*
* Lio works in atomic units (Bohr, Hartree, Electron mass) and as such Gromacs
* coordinates, forces and other variables need to be converted. In addition,
* since Lio is written in FORTRAN, arrays must be reordered.
*/
    int
        i, j;
    double
        qmener = 0, energ = 0, *qmcoord, *mmcoord, *mmgrad, *qmgrad, *mmcharge;
/*
* Reorders coordinates array to the format used in FORTRAN, treating it as a
* single-dimension array.
*/
    snew(qmcoord, 3*(qm->nrQMatoms));
    snew(mmcoord, 3*(mm->nrMMatoms));
    for (j = 0; j < 3; j++)
    {
        for (i = 0; i < qm->nrQMatoms; i++)
        {
            qmcoord[j+3*i] = qm->xQM[i][j]/(BOHR2NM);
        };
        for (i = 0; i < mm->nrMMatoms; i++)
        {
            mmcoord[j+3*i] = mm->xMM[i][j]/(BOHR2NM);
        };
    };

/* 
* The following lines are to be sure that MM partial charges are in a 
* double-precision variable. Then the actual SCF calculation is performed
* via the SCF_GRO function.
*/
    snew(mmcharge, (mm->nrMMatoms));
    for (i = 0; i < mm->nrMMatoms; i++)
    {
         mmcharge[i] = mm->MMcharges[i];
    };

    F77_FUNC(scf_gro, SCF_GRO) (&qmener, qmcoord, mmcoord, mmcharge,
                                &mm->nrMMatoms);
/*
* LIO separates forces from the QM and the MM parts, which Gromacs takes
* all together. As such, two routines need to be called and forces
* obtained are joined in f[] and fshift[] arrays.
*/
    snew(qmgrad, 3*(qm->nrQMatoms));
    snew(mmgrad, 3*(mm->nrMMatoms));

    F77_FUNC(dft_get_qm_forces, DFT_GET_QM_FORCES) (qmgrad);
    F77_FUNC(dft_get_mm_forces, DFT_GET_MM_FORCES) (mmgrad, qmgrad);

    for (j = 0; j < 3; j++)
    {
        for (i = 0; i < qm->nrQMatoms; i++)
        {
            f[i][j]      = qmgrad[3*i+j]*HARTREE_BOHR2MD;
            fshift[i][j] = qmgrad[3*i+j]*HARTREE_BOHR2MD;
        };
        for (i = 0; i < mm->nrMMatoms; i++)
        {
            f[i+qm->nrQMatoms][j]      = mmgrad[3*i+j]*HARTREE_BOHR2MD;
            fshift[i+qm->nrQMatoms][j] = mmgrad[3*i+j]*HARTREE_BOHR2MD;
        };
    };

    energ = qmener*HARTREE2KJ*AVOGADRO;

    free(mmcharge);
    free(qmcoord);
    free(mmcoord);
    free(qmgrad);
    free(mmgrad);
    
    return(energ);

}; /* call_lio */
#endif
