/*
** coarseGridMaker.h
**
** Automatically generated header file.
** Source file:   coarseGridMaker.c
** Last modified: Wednesday, June 06, 2007 at 16:58:48
*/


int	makeOptimalCoarseLayeredECL	(struct MASSAGER*);
void getDepthArray(struct MASSAGER *massager);
double hypercube_loop(int *ks_indx,int Nz, int *kFine, struct MASSAGER *massager);
double vertical_mismatch();
gint depth_info_compareXYLZ(gconstpointer a, gconstpointer b);
gint depth_info_compareLXYZ(gconstpointer a, gconstpointer b);
int readEclipseGridToSubgrid(struct MASSAGER *massager);
int exportCoarseGrid(struct MASSAGER *massager);
gboolean GHRCleanSegmentsFunc(gpointer key,gpointer value, gpointer user_data);
void flush_plateau_list();
void append_to_plateau_list(int N);
void string_hash_free(gpointer gs);



