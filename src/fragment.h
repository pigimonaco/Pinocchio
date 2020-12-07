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
  double R,q[3],v[3],D,Dv,z;
#ifdef TWO_LPT
  double D2,D2v,v2[3];
#ifdef THREE_LPT
  double D31,D31v,v31[3],D32,D32v,v32[3];
#endif
#endif
#ifdef RECOMPUTE_DISPLACEMENTS
  double w;
  double v_aft[3];
#ifdef TWO_LPT
  double v2_aft[3];
#ifdef THREE_LPT
  double v31_aft[3],v32_aft[3];
#endif
#endif
#endif
} pos_data;
pos_data obj, obj1, obj2;

typedef struct
{
  unsigned long long int name;
  double M,q[3],x[3],v[3];
  int n, pad;
}  catalog_data;

#ifdef RECOMPUTE_DISPLACEMENTS
typedef struct
{
  int n, mine;
  double z[MAXOUTPUTS];
} Segment_data;
Segment_data Segment;
#endif

void condition_for_accretion(int, int, int, int, double, int, double *, double *);
void condition_for_merging(double, int, int, int *);
void set_obj(int, double, pos_data *);
void set_obj_vel(int, double, pos_data *);
void set_point(int, int, int, int, double, pos_data *);
void set_group(int, pos_data *);
double q2x(int, pos_data *, int);
double vel(int, pos_data *);
double distance(int, pos_data *, pos_data *);
void clean_list(int *);
double virial(int, double, int);
void merge_groups(int, int, double);
void update_history(int, int, double);
void accretion(int, int, int, int, int, double);
void update(pos_data *, pos_data *);
int write_catalog(int);
int write_histories(void);
int compute_mf(int);
#ifdef PLC
int write_PLC();
#endif
