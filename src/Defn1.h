/*
 *  R : A Computer Language for Statistical Data Analysis
 *  Copyright (C) 1995, 1996  Robert Gentleman and Ross Ihaka
 *  Copyright (C) 1998--2008  The R Development Core Team.
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, a copy is available at
 *  http://www.r-project.org/Licenses/
 */
 /*
 *  Modified by Vasyl Zhabotynsky on 04/05/2011
 *  using part of Defn.h with slight alterations
 */
#ifndef DEFN_H_
#define DEFN_H_
#endif

/* Localization */

#ifdef ENABLE_NLS
#include <libintl.h>
#ifdef Win32
#define _(String) libintl_gettext (String)
#undef gettext /* needed for graphapp */
#else
#define _(String) gettext (String)
#endif
#define gettext_noop(String) String
#define N_(String) gettext_noop (String)
#else /* not NLS */
#define _(String) (String)
#define N_(String) String
#define ngettext(String, StringP, N) (N > 1 ? StringP: String)
#endif

// #include "Rinternals.h"
#include <Rinternals.h>
typedef double optimfn1(int, double *, void *, SEXP);
typedef void optimgr1(int, double *, double *, void *, SEXP);
