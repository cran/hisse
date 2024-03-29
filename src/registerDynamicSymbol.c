#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>


/* .C calls */
extern void initmod_bisse(void (* odeparms)(int *, double *));
extern void initmod_hisse(void (* odeparms)(int *, double *));
extern void initmod_hisse_null(void (* odeparms)(int *, double *));
extern void initmod_geosse(void (* odeparms)(int *, double *));
extern void initmod_geohisse(void (* odeparms)(int *, double *));
extern void initmod_noclass(void (* odeparms)(int *, double *));
extern void initmod_hinoclass(void (* odeparms)(int *, double *));

/* Introduced with new version with new models */
extern void initmod_musse(void (* odeparms)(int *, double *));
extern void initmod_muhisse(void (* odeparms)(int *, double *));
extern void initmod_fgeohisse(void (* odeparms)(int *, double *));
extern void initmod_fhinoclass(void (* odeparms)(int *, double *));
extern void initmod_fgeosse(void (* odeparms)(int *, double *));
extern void initmod_fnoclass(void (* odeparms)(int *, double *));
extern void initmod_misse(void (* odeparms)(int *, double *));
extern void initmod_fbisse(void (* odeparms)(int *, double *));
extern void initmod_fhisse(void (* odeparms)(int *, double *));

extern void maddison_DE_bisse(void *, void *, void *, void *, void *, void *);
extern void maddison_DE_hisse(void *, void *, void *, void *, void *, void *);
extern void maddison_DE_hisse_null(void *, void *, void *, void *, void *, void *);
extern void classe_geosse_equivalent_derivs(void *, void *, void *, void *, void *, void *);
extern void geohisse_derivs(void *, void *, void *, void *, void *, void *);
extern void notclasse_derivs(void *, void *, void *, void *, void *, void *);
extern void notclasse_more_derivs(void *, void *, void *, void *, void *, void *);
extern void musse_derivs(void *, void *, void *, void *, void *, void *);
extern void muhisse_derivs(void *, void *, void *, void *, void *, void *);
extern void fgeohisse_derivs(void *, void *, void *, void *, void *, void *);
extern void fnotclasse_more_derivs(void *, void *, void *, void *, void *, void *);
extern void fclasse_geosse_equivalent_derivs(void *, void *, void *, void *, void *, void *);
extern void fnotclasse_derivs(void *, void *, void *, void *, void *, void *);
extern void misse_derivs(void *, void *, void *, void *, void *, void *);
extern void misse_strat_derivs(void *, void *, void *, void *, void *, void *);
extern void maddison_DE_fbisse(void *, void *, void *, void *, void *, void *);
extern void maddison_DE_strat_fbisse(void *, void *, void *, void *, void *, void *);
extern void maddison_DE_fhisse(void *, void *, void *, void *, void *, void *);
extern void maddison_DE_strat_fhisse(void *, void *, void *, void *, void *, void *);

extern void set_birth_bisse_void(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void set_birth_hisse_null_void(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void set_birth_hisse_void(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);

static const R_CMethodDef CEntries[] = {
    {"initmod_bisse",                   (DL_FUNC) &initmod_bisse,                   1},
    {"initmod_hisse",                   (DL_FUNC) &initmod_hisse,                   1},
    {"initmod_hisse_null",              (DL_FUNC) &initmod_hisse_null,              1},
    {"initmod_geosse",		            (DL_FUNC) &initmod_geosse,                  1},
    {"initmod_geohisse",		        (DL_FUNC) &initmod_geohisse,                1},
    {"initmod_noclass",		            (DL_FUNC) &initmod_noclass,                 1},
    {"initmod_hinoclass",		        (DL_FUNC) &initmod_hinoclass,               1},
    {"initmod_musse",		            (DL_FUNC) &initmod_musse,                   1},
    {"initmod_muhisse",		            (DL_FUNC) &initmod_muhisse,                 1},
    {"initmod_fgeohisse",		        (DL_FUNC) &initmod_fgeohisse,               1},
    {"initmod_fhinoclass",		        (DL_FUNC) &initmod_fhinoclass,              1},
    {"initmod_fgeosse",		            (DL_FUNC) &initmod_fgeosse,                 1},
    {"initmod_fnoclass",		        (DL_FUNC) &initmod_fnoclass,                1},
    {"initmod_misse",		            (DL_FUNC) &initmod_misse,                   1},
    {"initmod_fbisse",                  (DL_FUNC) &initmod_fbisse,                  1},
    {"initmod_fhisse",                  (DL_FUNC) &initmod_fhisse,                  1},
    {"maddison_DE_bisse",               (DL_FUNC) &maddison_DE_bisse,               6},
    {"maddison_DE_hisse",               (DL_FUNC) &maddison_DE_hisse,               6},
    {"maddison_DE_hisse_null",          (DL_FUNC) &maddison_DE_hisse_null,          6},
    {"classe_geosse_equivalent_derivs", (DL_FUNC) &classe_geosse_equivalent_derivs, 6},
    {"geohisse_derivs",                 (DL_FUNC) &geohisse_derivs,                 6},
    {"notclasse_derivs",                (DL_FUNC) &notclasse_derivs,                6},
    {"notclasse_more_derivs",           (DL_FUNC) &notclasse_more_derivs,           6},
    {"musse_derivs",                    (DL_FUNC) &musse_derivs,                    6},
    {"muhisse_derivs",                  (DL_FUNC) &muhisse_derivs,                  6},
    {"fgeohisse_derivs",                (DL_FUNC) &fgeohisse_derivs,                6},
    {"fnotclasse_more_derivs",          (DL_FUNC) &fnotclasse_more_derivs,          6},
    {"fclasse_geosse_equivalent_derivs",(DL_FUNC) &fclasse_geosse_equivalent_derivs,6},
    {"fnotclasse_derivs",               (DL_FUNC) &fnotclasse_derivs,               6},
    {"misse_derivs",                    (DL_FUNC) &misse_derivs,                    6},
    {"misse_strat_derivs",              (DL_FUNC) &misse_strat_derivs,              6},
    {"maddison_DE_fbisse",              (DL_FUNC) &maddison_DE_fbisse,              6},
    {"maddison_DE_strat_fbisse",        (DL_FUNC) &maddison_DE_strat_fbisse,        6},
    {"maddison_DE_fhisse",              (DL_FUNC) &maddison_DE_fhisse,              6},
    {"maddison_DE_strat_fhisse",        (DL_FUNC) &maddison_DE_strat_fhisse,        6},
    {"set_birth_bisse_void",            (DL_FUNC) &set_birth_bisse_void,           29},
    {"set_birth_hisse_null_void",       (DL_FUNC) &set_birth_hisse_null_void,      54},
    {"set_birth_hisse_void",            (DL_FUNC) &set_birth_hisse_void,           59},
    {NULL, NULL, 0}
};

void R_init_hisse(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
