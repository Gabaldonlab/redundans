// Copyright 2015 Yutaro Konta

#include "LambdaCalculator.hh"
#include "cbrc_linalg.hh"
#include <cassert>
#include <cfloat>
#include <cmath>
#include <numeric>
//#include <iostream>
using namespace std;

static bool calculate_inv_sum(double **matrix, int alpha_size, double tau, double* inv_sum, double **tmpMat, double *tmpVec)
{
  for (int i=0; i<alpha_size; i++)
    for (int j=0; j<alpha_size; j++)
      tmpMat[i][j] = exp(tau * matrix[i][j]);

  std::fill_n(tmpVec, alpha_size, 1.0);

  try {
    cbrc::linalgSolve(tmpMat, tmpVec, alpha_size);
  } catch(...) {
    return false;
  }

  *inv_sum = std::accumulate(tmpVec, tmpVec + alpha_size, 0.0);
  return true;
}

namespace cbrc{

void LambdaCalculator::setBad(){
  lambda_ = -1;
  letterProbs1_.clear();
  letterProbs2_.clear();
}

bool LambdaCalculator::find_ub(double **matrix, int alpha_size, double *ub)
{
  double r_max_min = DBL_MAX;
  double c_max_min = DBL_MAX;
  double r_max;
  double c_max;

  double r_min;
  double c_min;

  int l_r = 0;
  int l_c = 0;

  for (int i=0; i<alpha_size; i++)
    {
      r_max = -DBL_MAX;
      r_min = DBL_MAX;
      for (int j=0; j<alpha_size; j++)
        {
          if (matrix[i][j] > r_max)
            r_max = matrix[i][j];

          if (matrix[i][j] < r_min)
            r_min = matrix[i][j];
        }
      if (r_max == 0 && r_min == 0)
        l_r++;
      else if (r_max <= 0 || r_min >= 0)
        return false;
      else if (r_max < r_max_min)
        r_max_min = r_max;
    }

  for (int j=0; j<alpha_size; j++)
    {
      c_max = -DBL_MAX;
      c_min = DBL_MAX;
      for (int i=0; i<alpha_size; i++)
        {
          if (matrix[i][j] > c_max)
            c_max = matrix[i][j];

          if (matrix[i][j] < c_min)
            c_min = matrix[i][j];
        }
      if (c_max == 0 && c_min == 0)
        l_c++;
      else if (c_max <= 0 || c_min >= 0)
        return false;
      else if (c_max < c_max_min)
        c_max_min = c_max;
    }

  if (l_r == alpha_size) return false;

  // the multiplication by 1.1 is sometimes necessary, presumably to
  // prevent the upper bound from being too tight:
  if (r_max_min > c_max_min)
    *ub = 1.1 * log(1.0 * (alpha_size - l_r))/r_max_min;
  else
    *ub = 1.1 * log(1.0 * (alpha_size - l_c))/c_max_min;

  return true;
}

bool LambdaCalculator::binary_search(double** matrix, int alpha_size, double lb, double ub, vector<double>& letprob1, vector<double>& letprob2, double* lambda, int maxiter)
{
  double l=0;
  double r=0;
  double l_sum=0;
  double r_sum=0;
  int iter=0;

  std::vector<double> tmpVals(alpha_size * alpha_size);
  std::vector<double *> tmpPtrs(alpha_size);
  for (int i = 0; i < alpha_size; ++i)
    tmpPtrs[i] = &tmpVals[i * alpha_size];
  double **tmpMat = &tmpPtrs[0];

  double *tmpVec = &letprob2[0];

  while (iter < maxiter && (l>=r || (l_sum < 1.0 && r_sum < 1.0) || (l_sum > 1.0 && r_sum > 1.0)))
    {
      l = lb + (ub - lb)*rand()/RAND_MAX;
      r = lb + (ub - lb)*rand()/RAND_MAX;

      if (!calculate_inv_sum(matrix, alpha_size, l, &l_sum, tmpMat, tmpVec) ||
	  !calculate_inv_sum(matrix, alpha_size, r, &r_sum, tmpMat, tmpVec))
        {
          l = 0;
          r = 0;
        }
      iter++;
    }

  if (iter >= maxiter)
    return false;

  while (l_sum != 1.0 && r_sum != 1.0 && (l+r)/2.0 != l && (l+r)/2.0 != r)
    {
      double mid = (l + r)/2.0;
      double mid_sum;
      if (!calculate_inv_sum(matrix, alpha_size, mid, &mid_sum, tmpMat, tmpVec))
        return false;

      if (fabs(mid_sum) >= DBL_MAX)
        return false;

      if ((l_sum < 1.0 && mid_sum >= 1.0) || (l_sum > 1.0 && mid_sum <= 1.0))
        {
          r = mid;
          r_sum = mid_sum;
        }

      else if ((r_sum < 1.0 && mid_sum >= 1.0) || (r_sum > 1.0 && mid_sum <= 1.0))
        {
          l = mid;
          l_sum = mid_sum;
        }

      else
        return false;
    }

  if (fabs(l_sum - 1.0) < fabs(r_sum - 1.0))
    {
      if (check_lambda(matrix, l, alpha_size, letprob1, letprob2, tmpMat))
        {
          *lambda = l;
          return true;
        }
      return false;
    }

  if (check_lambda(matrix, r, alpha_size, letprob1, letprob2, tmpMat))
    {
      *lambda = r;
      return true;
    }
  return false;
}

double LambdaCalculator::calculate_lambda(double** matrix, int alpha_size, vector<double>& letprob1, vector<double>& letprob2, int maxiter, int max_boundary_search_iter, double lb_ratio)
{
  double ub;

  if (!find_ub(matrix, alpha_size, &ub))
    return -1;

  double lb = ub*lb_ratio;
  double lambda = -1;
  int iter = 0;
  bool flag = false;

  while (!flag && iter < maxiter)
    {
      flag = binary_search(matrix, alpha_size, lb, ub, letprob1, letprob2, &lambda, max_boundary_search_iter);
      iter++;
    }

  return lambda;
}

bool LambdaCalculator::check_lambda(double** matrix, double lambda, int alpha_size, vector<double>& letprob1, vector<double>& letprob2, double** tmpMat)
{
  for (int i=0; i<alpha_size; i++)
    for (int j=0; j<alpha_size; j++)
      tmpMat[i][j] = exp(lambda * matrix[i][j]);

  std::fill_n(&letprob2[0], alpha_size, 1.0);
  cbrc::linalgSolve(tmpMat, &letprob2[0], alpha_size);

  for (int i=0; i<alpha_size; i++)
    {
      double p = letprob2[i];
      if (p < 0 || p > 1)
	return false;
      letprob2[i] = roundToFewDigits(p);
    }

  for (int i=0; i<alpha_size; i++)
    for (int j=0; j<alpha_size; j++)
      tmpMat[i][j] = exp(lambda * matrix[j][i]);

  std::fill_n(&letprob1[0], alpha_size, 1.0);
  cbrc::linalgSolve(tmpMat, &letprob1[0], alpha_size);

  for (int j=0; j<alpha_size; j++)
    {
      double q = letprob1[j];
      if (q < 0 || q > 1)
	return false;
      letprob1[j] = roundToFewDigits(q);
    }

  return true;
}

void LambdaCalculator::calculate(const const_int_ptr *matrix, int alphSize) {
  assert(alphSize >= 0);
  lambda_ = -1;
  letterProbs1_.resize(alphSize);
  letterProbs2_.resize(alphSize);

  int maxiter = 1000;
  int max_boundary_search_iter = 100;
  double lb_ratio = 1e-6;

  std::vector<double> matVals(alphSize * alphSize);
  std::vector<double *> matPtrs(alphSize);
  for (int i = 0; i < alphSize; ++i)
    matPtrs[i] = &matVals[i * alphSize];
  double** mat = &matPtrs[0];

  for (int i=0; i<alphSize; i++)
    for (int j=0; j<alphSize; j++)
      mat[i][j] = matrix[i][j];

  // xxx srand ?
  lambda_ = calculate_lambda(mat, alphSize, letterProbs1_, letterProbs2_, maxiter, max_boundary_search_iter, lb_ratio);
}
}
