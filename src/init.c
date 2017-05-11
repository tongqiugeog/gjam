#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* .Call calls */
extern SEXP gjam_byRcpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP gjam_condMVNRcpp(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP gjam_trMVNmatrixRcpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP gjam_rmvnormRcpp(SEXP, SEXP, SEXP);
extern SEXP gjam_solveRcpp(SEXP);
extern SEXP gjam_getPmatKRcpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP gjam_fnZRcpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP gjam_invWbyRcpp(SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
  {"gjam_byRcpp",           (DL_FUNC) &gjam_byRcpp,          6},
  {"gjam_condMVNRcpp",      (DL_FUNC) &gjam_condMVNRcpp,     5},
  {"gjam_trMVNmatrixRcpp",  (DL_FUNC) &gjam_trMVNmatrixRcpp, 7},
  {"gjam_rmvnormRcpp",      (DL_FUNC) &gjam_rmvnormRcpp,     3},
  {"gjam_solveRcpp",        (DL_FUNC) &gjam_solveRcpp,       1},
  {"gjam_getPmatKRcpp",     (DL_FUNC) &gjam_getPmatKRcpp,    7},
  {"gjam_fnZRcpp",          (DL_FUNC) &gjam_fnZRcpp,         8},
  {"gjam_invWbyRcpp",       (DL_FUNC) &gjam_invWbyRcpp,      2},
	{NULL, NULL, 0}
};

void R_init_gjam(DllInfo *dll) {
	R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
	R_useDynamicSymbols(dll, FALSE);
}
