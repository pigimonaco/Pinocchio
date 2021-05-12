/*****************************************************************
 *                        PINOCCHIO  V4.1                        *
 *  (PINpointing Orbit-Crossing Collapsed HIerarchical Objects)  *
 *****************************************************************
 
 This code was written by
 Pierluigi Monaco
 Copyright (C) 2016
 
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

#define NV 6
#define FILAMENT 1
#define SHIFT 0.5

#define ORDER_FOR_GROUPS 2
#define ORDER_FOR_CATALOG 3

extern int ngroups;

typedef struct
{
  int M,i;
  PRODFLOAT R,q[3],v[3],D,Dv;
  double z;
#ifdef TWO_LPT
  PRODFLOAT D2,D2v,v2[3];
#ifdef THREE_LPT
  PRODFLOAT D31,D31v,v31[3],D32,D32v,v32[3];
#endif
#endif
#ifdef RECOMPUTE_DISPLACEMENTS
  PRODFLOAT w;
  PRODFLOAT v_aft[3];
#ifdef TWO_LPT
  PRODFLOAT v2_aft[3];
#ifdef THREE_LPT
  PRODFLOAT v31_aft[3],v32_aft[3];
#endif
#endif
#endif
} pos_data;
pos_data obj, obj1, obj2;

typedef struct
{
  unsigned long long int name;
#ifdef LIGHT_OUTPUT
  PRODFLOAT M,x[3],v[3];
#else
  double M,q[3],x[3],v[3];
  int n,pad;
#endif
}  catalog_data;

typedef struct
{
  int n, mine;
  double z[MAXOUTPUTS];
} Segment_data;
Segment_data Segment;

void condition_for_accretion(int, int, int, int, int, PRODFLOAT, int, double *, double *); // LEVARE primo argomento
void condition_for_merging(PRODFLOAT, int, int, int *);
void set_obj(int, PRODFLOAT, pos_data *);
void set_obj_vel(int, PRODFLOAT, pos_data *);
void set_point(int, int, int, int, PRODFLOAT, pos_data *);
void set_group(int, pos_data *);
PRODFLOAT q2x(int, pos_data *, int);
PRODFLOAT vel(int, pos_data *);
PRODFLOAT distance(int, pos_data *, pos_data *);
void clean_list(int *);
PRODFLOAT virial(int, PRODFLOAT, int);
void merge_groups(int, int, PRODFLOAT);
void update_history(int, int, PRODFLOAT);
void accretion(int, int, int, int, int, PRODFLOAT);
void update(pos_data *, pos_data *);
int write_catalog(int);
int write_histories(void);
int compute_mf(int);
int find_location(int, int, int);
#ifdef PLC
int write_PLC();
void coord_transformation_cartesian_polar(PRODFLOAT *, double *, double *, double *);
#endif
