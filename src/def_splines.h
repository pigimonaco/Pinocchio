/*****************************************************************
 *                        PINOCCHIO  V5.1                        *
 *  (PINpointing Orbit-Crossing Collapsed HIerarchical Objects)  *
 *****************************************************************

 This code was written by
 Pierluigi Monaco, Tom Theuns, Giuliano Taffoni, Marius Lepinzan,
 Chiara Moretti, Luca Tornatore, David Goz, Tiago Castro
 Copyright (C) 2025

 github: https://github.com/pigimonaco/Pinocchio
 web page: http://adlibitum.oats.inaf.it/monaco/pinocchio.html

 This program is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation; either version 2 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program; if not, write to the Free Software
 Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

/* #if defined(SCALE_DEPENDENT) && defined(ELL_CLASSIC) */
/* #error Trying to compile with ELL_CLASSIC and SCALE_DEPENDENT together */
/* #endif */

#if defined(MOD_GRAV_FR) && !defined(SCALE_DEPENDENT)
#define SCALE_DEPENDENT
#endif

#ifndef SCALE_DEPENDENT
#define NkBINS 1
#else
#define NkBINS 10
#define LOGKMIN ((double)-3.0)
#define DELTALOGK ((double)0.5)
#endif

#define SP_TIME 0
#define SP_INVTIME 1
#define SP_COMVDIST 2
#define SP_DIAMDIST 3
#define SP_INVGROW 4
#define SP_MASSVAR 5
#define SP_DISPVAR 6
#define SP_RADIUS 7
#define SP_DVARDR 8

#define SP_GROW1 (9)
#define SP_GROW2 (9 + NkBINS)
#define SP_GROW31 (9 + 2 * NkBINS)
#define SP_GROW32 (9 + 3 * NkBINS)
#define SP_FOMEGA1 (9 + 4 * NkBINS)
#define SP_FOMEGA2 (9 + 5 * NkBINS)
#define SP_FOMEGA31 (9 + 6 * NkBINS)
#define SP_FOMEGA32 (9 + 7 * NkBINS)

#define SP_PK (9 + 8 * NkBINS)
#define SP_EOS (9 + 8 * NkBINS + 1)
#define SP_INTEOS (9 + 8 * NkBINS + 2)

#define NSPLINES_BASE (9 + 8 * NkBINS + 3)

#ifdef READ_HUBBLE_TABLE
#define SP_EXT_HUBBLE (NSPLINES_BASE)
#define NSPLINES (NSPLINES_BASE + 1)
#else
#define NSPLINES (NSPLINES_BASE)
#endif
