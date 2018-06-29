/**
 * @file
 * @author Keenan Crane <keenan@cs.caltech.edu>
 * @version 1.0
 *
 * @section LICENSE
 *
 * Copyright 2012 Keenan Crane. All rights reserved.
 * 
 * Redistribution and use in source and binary forms, with or without modification,
 * are permitted provided that the following conditions are met:
 * 
 * 1. Redistributions of source code must retain the above copyright notice, this
 *    list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 * 
 * THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS OR IMPLIED
 * WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
 * MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT
 * SHALL THE FREEBSD PROJECT OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
 * INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE
 * OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
 * ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 * 
 * The views and conclusions contained in the software and documentation are those
 * of the author and should not be interpreted as representing official policies,
 * either expressed or implied, of any other person or institution.
 *
 * @section DESCRIPTION
 *
 * This file contains global constants.
 *
 */

#define hmVersionString "version 1.0 (July 2012)"

#ifndef LIBGEODESIC_HMCONSTANTS_H
#define LIBGEODESIC_HMCONSTANTS_H

/* select a linear solver */
#define HM_USE_CHOLMOD
#undef  HM_USE_HSLMA87

/* linear solver parameters */
#ifdef HM_USE_HSLMA87

/* reordering schemes */
#define hmReorderAMD   1 /* approximate minimum degree */
#define hmReorderMD    2 /* minimum degree */
#define hmReorderMeTiS 3 /* MeTiS with default settings */
#define hmReorderIndef 4 /* MA47 ordering for indefinite matrices */

/* reordering scheme used for Cholesky factorization */
static const int hmReorderingScheme = hmReorderMeTiS;

#endif


/** \brief Default storage size for variable-length vectors. */
static const size_t hmVectorDefaultStorage = 16;

/** \brief Magnitude of diagonal regularization added to Laplacian to get numerical positive-definiteness (needed for CHOLMOD). */
static const double hmRegularization = 1e-10;

#endif /* LIBGEODESIC_HMCONSTANTS_H */

