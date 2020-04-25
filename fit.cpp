/***************************************************************************
 *   Copyright (C) 2020 by Nuno Cardoso                                    *
 *   nmrcardoso@gmail.com                                                  *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/



#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <limits>
#include <cstdlib>
#include <algorithm>
#include <numeric>
#include <cstring>
#include <cassert>
#include <map>
#include <vector>
#include <cmath>
#include <float.h>

#include <stdlib.h>
#include <stdio.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multifit_nlinear.h>


#include "fit.h"
#include "log.h"
#include "jackerr.h"


#ifndef uint
#define uint unsigned int
#endif


using namespace std;




struct data {
  size_t n;
  double * t;
  double * y;
};

int expb_f (const gsl_vector * x, void *data, gsl_vector * f);

int expb_df (const gsl_vector * x, void *data, gsl_matrix * J);

void callback(const size_t iter, void *params, const gsl_multifit_nlinear_workspace *w);


GSLfitRes GSLfit1(DataLine data);






int expb_f (const gsl_vector * x, void *data, gsl_vector * f){
  size_t n = ((struct data *)data)->n;
  double *t = ((struct data *)data)->t;
  double *y = ((struct data *)data)->y;

  double A = gsl_vector_get (x, 0);
  double B = gsl_vector_get (x, 1);
  double C = gsl_vector_get (x, 2);

  size_t i;

  for (i = 0; i < n; i++){
      /* Model Yi = A * exp(-lambda * t_i) + b */
      double Yi = A + B/t[i] +C*t[i];
      gsl_vector_set (f, i, Yi - y[i]);
    }

  return GSL_SUCCESS;
}

int expb_df (const gsl_vector * x, void *data, gsl_matrix * J){
  size_t n = ((struct data *)data)->n;
  double *t = ((struct data *)data)->t;

  double A = gsl_vector_get (x, 0);
  double B = gsl_vector_get (x, 1);
  double C = gsl_vector_get (x, 2);

  size_t i;

  for (i = 0; i < n; i++) {
      /* Jacobian matrix J(i,j) = dfi / dxj, */
      /* where fi = (Yi - yi)/sigma[i],      */
      /*       Yi = A * exp(-lambda * t_i) + b  */
      /* and the xj are the parameters (A,lambda,b) */
      gsl_matrix_set (J, i, 0, 1.0);
      gsl_matrix_set (J, i, 1, 1.0/t[i]);
      gsl_matrix_set (J, i, 2, t[i]);
    }

  return GSL_SUCCESS;
}

void callback(const size_t iter, void *params, const gsl_multifit_nlinear_workspace *w){
	gsl_vector *f = gsl_multifit_nlinear_residual(w);
	gsl_vector *x = gsl_multifit_nlinear_position(w);
	double rcond;

	// compute reciprocal condition number of J(x)
	gsl_multifit_nlinear_rcond(&rcond, w);

	qlog << "iter " << iter;
	qlog << ": A = " << gsl_vector_get(x, 0);
	qlog << ", B = " << gsl_vector_get(x, 1);
	qlog << ", C = " << gsl_vector_get(x, 2);
	qlog << ", cond(J) = " << 1.0 / rcond;
	qlog << ", |f(x)| = " << gsl_blas_dnrm2(f) << endl;
}

	

GSLfitRes GSLfit1(DataLine data, bool print_info){
	const gsl_multifit_nlinear_type *T = gsl_multifit_nlinear_trust;
	gsl_multifit_nlinear_workspace *w;
	gsl_multifit_nlinear_fdf fdf;
	gsl_multifit_nlinear_parameters fdf_params =
	gsl_multifit_nlinear_default_parameters();
	const size_t n = data.size();
	const size_t p = 3;

	gsl_vector *f;
	gsl_matrix *J;
	gsl_matrix *covar = gsl_matrix_alloc (p, p);
	double * t = new double[n];
	double * y = new double[n];
	double * weights = new double[n];

	struct data d = { n, t, y };
	double x_init[3] = { 1.0, -0.6, .1 }; /* starting values */
	gsl_vector_view x = gsl_vector_view_array (x_init, p);
	gsl_vector_view wts = gsl_vector_view_array(weights, n);
	gsl_rng * r;
	double chisq, chisq0;
	int status, info;
	size_t i;

	const double xtol = 1e-10;
	const double gtol = 1e-10;
	const double ftol = 0.0;

	gsl_rng_env_setup();
	r = gsl_rng_alloc(gsl_rng_default);

	/* define the function to be minimized */
	fdf.f = expb_f;
	fdf.df = expb_df;   /* set to NULL for finite-difference Jacobian */
	fdf.fvv = NULL;     /* not using geodesic acceleration */
	fdf.n = n;
	fdf.p = p;
	fdf.params = &d;

	/* this is the data to be fitted */
	for (i = 0; i < n; i++){
		double ti = data[i].x;
		double yi = data[i].y;
		double si = data[i].error;
		double dy = gsl_ran_gaussian(r, si);

		t[i] = ti;
		y[i] = yi;
		weights[i] = 1.0 / (si * si);
		//printf ("data: %g %g %g\n", ti, y[i], si);
	}

	/* allocate workspace with default parameters */
	w = gsl_multifit_nlinear_alloc (T, &fdf_params, n, p);

	/* initialize solver with starting point and weights */
	gsl_multifit_nlinear_winit (&x.vector, &wts.vector, &fdf, w);

	/* compute initial cost function */
	f = gsl_multifit_nlinear_residual(w);
	gsl_blas_ddot(f, f, &chisq0);

	/* solve the system with a maximum of 100 iterations */
	if(print_info) status = gsl_multifit_nlinear_driver(100, xtol, gtol, ftol, callback, NULL, &info, w);
	else status = gsl_multifit_nlinear_driver(100, xtol, gtol, ftol, NULL, NULL, &info, w);

	/* compute covariance of best fit parameters */
	J = gsl_multifit_nlinear_jac(w);
	gsl_multifit_nlinear_covar (J, 0.0, covar);

	/* compute final cost */
	gsl_blas_ddot(f, f, &chisq);

	#define FIT(i) gsl_vector_get(w->x, i)
	#define ERR(i) sqrt(gsl_matrix_get(covar,i,i))

	if(print_info){
		qlog << "summary from method '" << gsl_multifit_nlinear_name(w) << "/" << gsl_multifit_nlinear_trs_name(w) << "'" << endl;
		qlog << "number of iterations: " << gsl_multifit_nlinear_niter(w) << endl;
		qlog << "function evaluations: " << fdf.nevalf << endl;
		qlog << "Jacobian evaluations: " << fdf.nevaldf << endl;
		qlog << "reason for stopping: " << ((info == 1) ? "small step size" : "small gradient") << endl;
		qlog << "initial |f(x)| = " << sqrt(chisq0) << endl;
		qlog << "final   |f(x)| = " << sqrt(chisq) << endl;
	}

	GSLfitRes res;
	{
		double dof = n - p;
		double c = GSL_MAX_DBL(1, sqrt(chisq / dof));

		if(print_info) {
			qlog << "GSL info (not jackknife errors)" << endl;
			qlog << "chisq/dof = " << chisq / dof << endl;
			qlog << "A      = " << FIT(0) << " ± " << c*ERR(0) << endl;
			qlog << "B      = " << FIT(1) << " ± " << c*ERR(1) << endl;
			qlog << "C      = " << FIT(2) << " ± " << c*ERR(2) << endl;
		}
		for(int i=0;i<3;++i){
			res.val[i] = FIT(i);
			res.error[i] = c*ERR(i);
		}
		res.chi2_dof = chisq / dof;
	}
	gsl_strerror (status);
	if(print_info) qlog << "status = " << gsl_strerror (status) << endl;

	gsl_multifit_nlinear_free (w);
	gsl_matrix_free (covar);
	gsl_rng_free (r);
	delete[] t,y,weights;
	return res;
}





GSLfitRes GSLfit(DataLine data, bool err_mean_main_fit){
	qlog << "------------------- GSL Fit -------------------" << endl;
	//Fit data
	GSLfitRes res = GSLfit1(data, true);

	//Jackknife fit error
	vector<GSLfitRes> err;
	vector<double> r0_jack;
	for( uint ii = 0; ii < data.size(); ++ii ) {
		DataLine data1;
		for( uint i = 0; i < data.size(); ++i ) {			
			if ( ii != i )
				data1.push_back(data[i]);
		}
		GSLfitRes res1 = GSLfit1(data1, false);
		err.push_back(res1);
		double r0 = std::sqrt( (1.65 + res1.val[1]) / res1.val[2] );
		r0_jack.push_back(r0);
		//qlog << ii << "\t" << res1.val[0] << "\t" << res1.val[1] << "\t" << res1.val[2] << endl;
	}

	double err_mean[3];
	double errJ_mean[3];
	for( uint i = 0; i < 3; ++i ){
		vector<double> vals;
		for(int j = 0; j < err.size(); ++j)
			vals.push_back(err[j].val[i]);
		if( err_mean_main_fit ) res.errorJ[i] = jackerr(vals, res.val[i]);
		res.errorJ[i] = jackerr(vals);
		//qlog << jackerr(vals) << "::" << jackerr(vals, res.val[i]) << endl;
	}
	
	//chi2/dof
	double qui = 0.0;
	for( uint i = 0; i < data.size(); ++i ) {	
		double yy = res.val[0] + res.val[1] / data[i].x + res.val[2] * data[i].x;
		double dif = data[i].y - yy;
		double inverr2 = 1.0 / (data[i].error * data[i].error);
		qui += dif * dif * inverr2;
	}
	//qlog << "qui: " << qui << endl;
	res.chi2_dofJ = qui /(data.size() - 3);
	
	
	double r0 = std::sqrt( (1.65 + res.val[1]) / res.val[2] );
	double r0_error = jackerr(r0_jack);
	qlog << "Using jackknife:" << endl;
	qlog << "r0 = " << r0 << " ± " << r0_error << endl;
	//qlog << "r0 = " << r0 << " ± " << jackerr(r0_jack, r0) << endl;
	
	double eb = res.errorJ[1] / (2.0 * res.val[2] * r0);
	double ec = res.errorJ[2] * r0 / (2.0 * res.val[2]);
	double roerr = std::sqrt( eb * eb + ec * ec);
	qlog << "Using propagation error:" << endl;
	qlog << "r0 = " << r0 << " ± " << roerr << endl;
	

	qlog << "-----------------------------------------------" << endl;
	return res;
}


