////////////////////////////////////////////////////////////////////// 
// libsrc/MathStats.cpp 
// (c) 2000-2011 Goncalo Abecasis
// 
// This file is distributed as part of the Goncalo source code package   
// and may not be redistributed in any form, without prior written    
// permission from the author. Permission is granted for you to       
// modify this file for your own personal use, but modified versions  
// must retain this copyright notice and must not be distributed.     
// 
// Permission is granted for you to use this file to compile Goncalo.    
// 
// All computer programs have bugs. Use this file at your own risk.   
// 
// Friday March 25, 2011
// 
 
#include "MathConstant.h"
#include "MathStats.h"
#include "Error.h"

#include <math.h>
#include <stdlib.h>
#include <stdio.h>

// Approximates the sqrt of an integer
// (should be within +/- 1 for inputs < 1,000,000)

int introot(int n)
   {
   int root = 1, scale = 1;

   while (scale < n)
       root *= 3, scale *= 10;

   for (int i = 1; i < 4; i++)
       root = (root + n / root) / 2;

   return root;
   }

// The normal distribution function
//

#define LOWER_TAIL_ONE     7.5
#define UPPER_TAIL_ZERO    37.5

double ndist(double z, bool upper)
   {
   // C version of ID Hill, "The Normal Integral"
   // Applied Statistics, Vol 22, pp. 424-427

   // If 7 digit accuracy is enough, alternative is
   // return erfcc(x / M_SQRT2) * 0.5;

   if (z < 0)
      {
      upper = !upper;
      z = -z;
      }

   if (z > LOWER_TAIL_ONE && !upper || z > UPPER_TAIL_ZERO)
      return (upper) ? 0.0 : 1.0;

   double p, y = 0.5 * z * z;

   if (z < 1.28)
      {
      p = 0.5 - z * (0.398942280444 - 0.399903438504 * y /
                    (y + 5.75885480458 - 29.8213557808 /
                    (y + 2.62433121679 + 48.6959930692 /
                    (y + 5.92885724438))));
      }
   else
      {
      p = 0.398942270385 * exp (-y) /
          (z - 2.8052e-8 + 1.00000615302 /
          (z + 3.98064794e-4 + 1.98615381364 /
          (z - 0.151679116635 + 5.29330324926 /
          (z + 4.8385912808 - 15.1508972451 /
          (z + 0.742380924027 + 30.789933034 /
          (z + 3.99019417011))))));
      }

   return (upper) ? p : 1 - p;
   }

double logndist(double z, bool upper)
   {
   // Based on original FORTRAN code for ANORM by W. Cody

   if (z < 0)
      {
      upper = !upper;
      z = -z;
      }

   if (!upper || z < sqrt(32))
      return log(ndist(z, upper));

   double p[] = {0.21589853405795699,  0.1274011611602473639,
                 0.022235277870649807,  0.001421619193227893466,
                 2.9112874951168792E-5,  0.02307344176494017303 };

   double q[] = {1.28426009614491121,    0.468238212480865118,
                 0.0659881378689285515,  0.00378239633202758244,
                 7.29751555083966205E-5};

   double onexsq = 1.0 / (z * z);
   double xnum = p[5] * onexsq;
   double xden = onexsq;

   for (int i = 0; i < 4; i++)
      {
      xnum = (xnum + p[i]) * onexsq;
      xden = (xden + q[i]) * onexsq;
      }

   double result = onexsq * (xnum + p[4]) / (xden + q[4]);
   result = (0.39894228040143267794 - result) / z;

   onexsq = floor (z * 16) / 16;
   double delta = (z - onexsq) * (z + onexsq);

   return (-onexsq * onexsq * 0.5) + (-delta * 0.5) + log(result);
   }


// The standard normal distribution
//

double ninv( long double p )
/****************************************************
   C Equivalent of Wichura's PPND16, Algorithm AS241
   Applied Statistics Vol 37 1988 pp 477 - 484
*****************************************************/
{
   const double   SPLIT1 = 0.425,
                  SPLIT2 = 5.0,
                  CONST1 = 0.180625,
                  CONST2 = 1.6;

   static const double a[8] = {
                           3.3871328727963666080E0,
                           1.3314166789178437745E2,
                           1.9715909503065514427E3,
                           1.3731693765509461125E4,
                           4.5921953931549871457E4,
                           6.7265770927008700853E4,
                           3.3430575583588128105E4,
                           2.5090809287301226727E3
                        } ;

   static const double b[7] = {
                           4.2313330701600911252E1,
                           6.8718700749205790830E2,
                           5.3941960214247511077E3,
                           2.1213794301586595867E4,
                           3.9307895800092710610E4,
                           2.8729085735721942674E4,
                           5.2264952788528545610E3
                        } ;

   static const double c[8] = {
                           1.42343711074968357734E0,
                           4.63033784615654529590E0,
                           5.76949722146069140550E0,
                           3.64784832476320460504E0,
                           1.27045825245236838258E0,
                           2.41780725177450611770E-1,
                           2.27238449892691845833E-2,
                           7.74545014278341407640E-4
                        } ;

   static const double d[7] = {
                          2.05319162663775882187E0,
                           1.67638483018380384940E0,
                           6.89767334985100004550E-1,
                           1.48103976427480074590E-1,
                           1.51986665636164571966E-2,
                           5.47593808499534494600E-4,
                           1.05075007164441684324E-9
                        } ;

   static const long double e[8] = {
                           6.65790464350110377720E0,
                           5.46378491116411436990E0,
                           1.78482653991729133580E0,
                           2.96560571828504891230E-1,
                           2.65321895265761230930E-2,
                           1.24266094738807843860E-3,
                           2.71155556874348757815E-5,
                           2.01033439929228813265E-7
                        } ;

   static const long double f[7] = {
                           5.99832206555887937690E-1,
                           1.36929880922735805310E-1,
                           1.48753612908506148525E-2,
                           7.86869131145613259100E-4,
                           1.84631831751005468180E-5,
                           1.42151175831644588870E-7,
                           2.04426310338993978564E-15
                        } ;

   long double q = p - 0.5;
   long double r, x ;

   if ( fabs( q ) < SPLIT1 ) {
      r = CONST1 - q * q ;
      return q * ((((((( a[7] * r + a[6] ) * r + a[5] ) * r + a[4] ) * r
                        + a[3] ) * r + a[2] ) * r + a[1] ) * r + a[0] ) /
                 ((((((( b[6] * r + b[5] ) * r + b[4] ) * r + b[3] ) * r
                        + b[2] ) * r + b[1] ) * r + b[0] ) * r + 1.0 ) ;
   } else {
      if ( q < 0 )
         r = p ;
      else
         r = 1.0 - p ;

#ifdef __BORLANDC__
   #define NINV_LONG_DOUBLE
#endif

#ifdef __GNUC__
   #define NINV_LONG_DOUBLE
#endif

#ifdef NINV_LONG_DOUBLE
      if ( r < 1e-1000L )
#else
      if ( r < 1e-300 )
#endif
         error("p-value [%.2g] outside range in ninv()", p );

      if ( r > 0.0 )
         {
#ifdef NINV_LONG_DOUBLE
         r = sqrtl( -logl( r ) ) ;
#else
         r = sqrt( -log( r ) );
#endif
         if ( r <= SPLIT2 )
            {
            r -= CONST2 ;
            x = ((((((( c[7] * r + c[6] ) * r + c[5] ) * r + c[4] ) * r
                        + c[3] ) * r + c[2] ) * r + c[1] ) * r + c[0] ) /
                  ((((((( d[6] * r + d[5] ) * r + d[4] ) * r + d[3] ) * r
                        + d[2] ) * r + d[1] ) * r + d[0] ) * r + 1.0 ) ;
            }
         else
            {
            r -= SPLIT2 ;
            x =  ((((((( e[7] * r + e[6] ) * r + e[5] ) * r + e[4] ) * r
                        + e[3] ) * r + e[2] ) * r + e[1] ) * r + e[0] ) /
                  ((((((( f[6] * r + f[5] ) * r + f[4] ) * r + f[3] ) * r
                        + f[2] ) * r + f[1] ) * r + f[0] ) * r + 1.0L ) ;
            }
         }
      else
         x = HUGE_VAL;

      if ( q < 0 )
         x = -x ;
      return x ;
   }
}

// The chi-squared distribution
//
double chidist(double x, double v)
   { return gammq (0.5 * v, 0.5 * x); }

// The non-central chi-squared distribution
//
double chidist(double x, double f, double theta)
   {
   if (x < 0.0 || f < 0.0 || theta < 0.0)
      error("Invalid arguments in chidist function");

   if (x == 0.0)
      return 1.0;

   // Evaluate the first term in series
   int n = 1;

   double lambda = theta * 0.5;
   double u = exp(-lambda);
   double v = u;
   double x2 = x * 0.5, f2 = f * 0.5;
   double t = pow(x2, f2) * exp(-x2) / exp(gammln(f2 + 1.0));

   double result = v * t;

   // Initial approximation
   while ( f + 2.0 * n < x )
      {
      u *= lambda / n;
      v += u;
      t *= x / (f + 2.0 * n);
      result += v * t;
      n++;
      }

   // Loop until we have accurate result or exceed ITMAX
   while (t * x / (f + 2.0 * n - x) > 1e-10)
      {
      if (n > ITMAX)
         error("chidist function did not converge within %d iterations\n");

      u *= lambda / n;
      v += u;
      t *= x / (f + 2.0 * n);
      result += v * t;
      n++;
      }

   return 1.0 - result;
   }

// The error function
//
double erff(double x)
   {
   return x < 0.0 ? -gammp(0.5, x*x) : gammp(0.5, x*x);
   }

double erffc(double x)
   {
   return x < 0.0 ? 1.0 + gammp(0.5, x*x) : gammq(0.5, x*x);
   }

double erfcc(double x)
// returns the complementary of the error function erfc(x),
// with a fractional error everywhere of less than
// 1.2 * 10-7
   {
   double t, z, ans;

   z = fabs(x);
   t = 1.0 / (1.0 + 0.5 *    z);
   ans = t*exp(-z*z -1.26551223 +t*(1.00002368 +t*(0.37409196 +t*(0.09678418
         +t*(-0.18628806 +t*(0.27886807 +t*(-1.13520398 +t*(1.48851587
         +t*(-0.82215223 +t*0.17087277)))))))));
   return (x >= 0.0 ? ans : 2.0 - ans);
   }

// The f-distribution
//

double fdist(double x, double v1, double v2)
   {
   return betai(v2/2, v1/2, v2/(v2+v1*x));
   }

// The student's T-distribution
double tdist(double x, double df)
   {
   return betai(df * 0.5, 0.5, df/(df + x*x));
   }

// Gamma distribution functions
//

double gammln ( double xx )
   {
   double x, y, tmp, ser;
   static double cof[6] = { 76.18009172947146, -86.50532032941677, 24.01409824083091,
                      -1.231739572450155, 0.1208650973866179e-2, -0.5395239384953e-5 };
   int j;

   y = x = xx;

   tmp = x + 5.5;

   tmp -= ( x + 0.5 ) * log ( tmp );

   ser = 1.000000000190015;

   for ( j=0; j<=5; j++) ser += cof[j]/++y;

   return - tmp + log ( 2.5066282746310005 * ser/x);
   }

double gammp ( double a, double x )
   {
   double gamser, gamcf, gln;

   if ( x < 0.0 || a <= 0.0)
      error("Invalid arguments in routine gammp");
   if ( x < (a + 1.0))
      {
      gser(&gamser, a, x, &gln);    // use series representation
      return gamser;
      }
   else
      {
      gcf(&gamcf, a, x, &gln);     // use the continued fraction representation
      return 1.0 - gamcf;          // and take its complement
      }
   }

double gammq ( double a, double x )
   {
   double gammser, gammcf, gln;

   if (x < 0.0 || a <= 0.0) error("Invalid arguments in routine gammq");
   if (x < (a + 1.0))
      {
      gser (&gammser, a, x, &gln);     // use the series representation
      return (1.0 - gammser);          // and take its complement
      }
   else
      {
      gcf ( &gammcf, a, x, &gln);      // use the continued fraction representation
      return gammcf;
      }
   }

void gser ( double * gamser, double a, double x, double * gln)
   {
   int   n;
   double sum, del, ap;

   *gln=gammln(a);
   if (x <= 0.0)
      {
      if (x < 0.0) error("x less than 0 in gamma series routine (gser)");
      *gamser = 0.0;
      return;
      }
   else
      {
      ap = a;
      del = sum = 1.0 / a;
      for (n = 1; n <= ITMAX; n++)
         {
         ++ ap;
         del *= x / ap;
         sum += del;
         if ( fabs(del) < fabs(sum) * EPS )
            {
            *gamser = sum * exp ( -x + a * log (x) - (*gln));
            return;
            }
         }
      error("a too large, ITMAX too small in gamma series routine (gser)");
      return;
      }
   }

void gcf ( double * gammcf, double a, double x, double * gln)
   {
   int i;
   double an, b, c, d, del, h;

   *gln = gammln(a);

   b = x + 1.0 - a;           // Setup for evaluating continued fraction by
   c = 1.0 / FPMIN;           // Lentz method (cf NRC 5.2)
   d = 1.0 / b;
   h = d;

   for ( i = 1; i <= ITMAX; i++ )   // Iterate to convergence
      {
      an = -i * ( i - a );
      b += 2.0;
      d = an * d + b;
      if ( fabs(d) < FPMIN ) d = FPMIN;
      c = b + an / c;
      if ( fabs(c) < FPMIN ) c = FPMIN;
      d = 1.0 / d;
      del = d * c;
      h *= del;
      if ( fabs(del-1.0) < EPS ) break;
      }
   if ( i > ITMAX ) error ("a too large, ITMAX too small in gamma countinued fraction (gcf)");
   *gammcf = exp (-x + a*log(x) - (*gln)) * h;
   }

// Beta functions
//

double betai(double a, double b, double x)
   {
   double bt;

   if ( x < 0.0 || x > 1.0) error("betai: Bad x");
   if ( x == 0.0 || x == 1.0)
      bt = 0.0;
   else
      bt = exp(gammln(a+b)-gammln(a)-gammln(b)+a*log(x)+b*log(1.0-x));
   if (x < (a + 1.0)/(a + b + 2.0))
      // use continued fraction directly
      return bt*betacf(a,b,x)/a;
   else
      // use continued fraction after making the symmetry transformation
      return 1.0-bt*betacf(b,a,1.0-x)/b;
   }

double betacf(double a, double b, double x)
   {
   int m, m2;
   double aa, c, d, del, h, qab, qam, qap;

   // these q's will be used in factors that appear in coefficients
   qab = a + b;
   qap = a + 1.0;
   qam = a - 1.0;

   // First step of Lentz's method
   c = 1.0;
   d = 1.0 - qab*x/qap;
   if (fabs(d) < FPMIN) d=FPMIN;
   d = 1.0 / d;
   h = d;
   for (m=1; m<=ITMAX; m++)
      {
      m2 = 2*m;

      // The even step of the recurrence
      aa = m * (b-m)*x/((qam+m2)*(a+m2));
      d = 1.0 + aa*d;
      if (fabs(d) < FPMIN) d=FPMIN;
      c = 1.0 + aa/c;
      if (fabs(c) < FPMIN) c=FPMIN;
      d = 1.0/d;
      h *= d*c;

      // The odd step of the recurrence
      aa = -(a+m)*(qab+m)*x/((a+m2)*(qap+m2));
      d = 1.0+aa*d;
      if (fabs(d) < FPMIN) d = FPMIN;
      c = 1.0+aa/c;
      if (fabs(c) < FPMIN) c = FPMIN;
      d = 1.0/d;
      del = d*c;
      h *= del;

      // Are we done?
      if (fabs(del - 1.0) < EPS) break;
      }
   if (m > ITMAX)
      error("betacf: a or b too big or ITMAX too small");
   return h;
   }


// Bivariate normal distribution (assuming mean = 0 and sd = 1) functions
//

//  Translated from Alanz Genz's FORTRAN code in mvtnorm R package.
//  Normal distribution probabilities accurate to 1d-15.
//  Reference: J.L. Schonfelder, Math Comp 32(1978), pp 1232-1240.
//
double normp(double z) {
    int i, im;

    double bm, b, bp, p, rtwo, t, xa;

    rtwo = 1.414213562373095048801688724209;
    im = 24;

    double a[44] = { 6.10143081923200417926465815756e-1, -4.34841272712577471828182820888e-1,
                     1.76351193643605501125840298123e-1, -6.0710795609249414860051215825e-2,
                     1.7712068995694114486147141191e-2, -4.321119385567293818599864968e-3,
                     8.54216676887098678819832055e-4, -1.27155090609162742628893940e-4,
                     1.1248167243671189468847072e-5, 3.13063885421820972630152e-7,
                     -2.70988068537762022009086e-7, 3.0737622701407688440959e-8,
                     2.515620384817622937314e-9, -1.028929921320319127590e-9,
                     2.9944052119949939363e-11, 2.6051789687266936290e-11,
                     -2.634839924171969386e-12, -6.43404509890636443e-13,
                     1.12457401801663447e-13, 1.7281533389986098e-14,
                     -4.264101694942375e-15, -5.45371977880191e-16,
                     1.58697607761671e-16, 2.0899837844334e-17,
                     -5.900526869409e-18, -9.41893387554e-19, 2.14977356470e-19,
                     4.6660985008e-20, -7.243011862e-21, -2.387966824e-21,
                     1.91177535e-22, 1.20482568e-22, -6.72377e-25, -5.747997e-24,
                     -4.28493e-25, 2.44856e-25, 4.3793e-26, -8.151e-27, -3.089e-27,
                     9.3e-29, 1.74e-28, 1.6e-29, -8.0e-30, -2.0e-30
    };

    xa = fabs(z) / rtwo;
    if (xa > 100) {
        p = 0;
    } else {
        t = ( 8 * xa - 30 ) / ( 4 * xa + 15 );
        bm = 0;
        b  = 0;
        for (int i = im; i >= 0; --i) {
            bp = b;
            b = bm;
            bm = t * b - bp  + a[i];
        }
        p = exp( -xa * xa ) * ( bm - bp ) / 4;
    }
    if (z > 0) {
        p = 1 - p;
    }
    return p;
}

// bivariate normal probability density function
double binormp(double x, double y, double rho) {
    if (fabs(rho) == 1.0) { // degenerate case
        if ((x == 0.0) && (y == 0.0)) {
            return std::numeric_limits<double>::infinity();
        }
        return 0.0;
    }
    double c = 1.0 - rho * rho ;
    return (1.0 / (2.0 * M_PI * sqrt(c))) * exp(-(x * x - 2.0 * rho * x * y + y * y) / (2 * c));
}

// Translated from Alanz Genz's FORTRAN code in mvtnorm R package.
// cumulative distribution function
// calculates probability that X > lx and Y > ly
double binormq(double lx, double ly, double rho) {
    double bvn = 0.0;
    double twopi = 6.283185307179586;

    int lg, ng;

    double as, a, b, c, d, rs, xs;
    double sn, asr, h, k, bs, hs, hk;

    // Gauss Legendre Points and Weights
    double x[3][10] = {
            {-0.9324695142031522, -0.6612093864662647, -0.2386191860831970, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
            {-0.9815606342467191, -0.9041172563704750, -0.7699026741943050, -0.5873179542866171, -0.3678314989981802, -0.1252334085114692, 0.0, 0.0, 0.0, 0.0},
            {-0.9931285991850949, -0.9639719272779138, -0.9122344282513259, -0.8391169718222188, -0.7463319064601508, -0.6360536807265150, -0.5108670019508271, -0.3737060887154196, -0.2277858511416451, -0.7652652113349733e-1}
    };

    // Gauss Legendre Points and Weights
    double w[3][10] = {
            {0.1713244923791705, 0.3607615730481384, 0.4679139345726904, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
            {0.4717533638651177e-1, 0.1069393259953183, 0.1600783285433464, 0.2031674267230659, 0.2334925365383547, 0.2491470458134029, 0.0, 0.0, 0.0, 0.0},
            {0.1761400713915212e-1, 0.4060142980038694e-1, 0.6267204833410906e-1, 0.8327674157670475e-1, 0.1019301198172404, 0.1181945319615184, 0.1316886384491766, 0.1420961093183821, 0.1491729864726037, 0.1527533871307259}
    };

    if (fabs(rho) < 0.3) {
        ng = 0;
        lg = 3;
    } else if (fabs(rho) < 0.75) {
        ng = 1;
        lg = 6;
    } else {
        ng = 2;
        lg = 10;
    }

    h = lx;
    k = ly;
    hk = h * k;

    if (fabs(rho) < 0.925) {
        hs = ( h * h + k * k ) / 2;
        asr = asin(rho);
        for (int i = 0; i < lg; ++i) {
            sn = sin(asr * ( x[ng][i] + 1.0 ) / 2.0 );
            bvn = bvn + w[ng][i] * exp( ( sn * hk - hs ) / ( 1.0 - sn * sn ) );
            sn = sin (asr * (-x[ng][i] + 1 ) / 2.0);
            bvn = bvn + w[ng][i] * exp( ( sn * hk - hs ) / ( 1.0 - sn * sn ) );
        }
        return bvn * asr / (2.0 * twopi) + normp(-h) * normp(-k);
    } else {
        if (rho < 0) {
            k = -k;
            hk = -hk;
        }
        if (fabs(rho) < 1) {
            as = ( 1.0 - rho )*( 1.0 + rho );
            a = sqrt(as);
            bs = pow( h - k, 2.0);
            c = ( 4 - hk ) / 8;
            d = ( 12 - hk ) / 16;
            bvn = a * exp( -(bs / as + hk) / 2 ) * ( 1 - c * (bs - as) * (1 - d * bs / 5) / 3 + c * d * as * as / 5 );
            if (hk > -160) {
                b = sqrt(bs);
                bvn = bvn - exp(-hk / 2) * sqrt(twopi) * normp(-b / a) * b * ( 1 - c * bs *( 1 - d * bs / 5 ) / 3 );
            }
            a = a / 2;
            for (int i = 0; i < lg; ++i) {
                xs = pow( a * (x[ng][i] + 1), 2.0);
                rs = sqrt( 1 - xs );
                bvn = bvn + a * w[ng][i] * ( exp( -bs/ ( 2 * xs) - hk / (1 + rs) )/ rs - exp( -(bs/xs + hk) / 2 )*( 1 + c * xs *( 1 + d * xs ) ) );
                xs = as * pow(-x[ng][i] +1, 2) / 4;
                rs = sqrt( 1 - xs );
                bvn = bvn + a * w[ng][i] * exp( -(bs / xs + hk) / 2 ) * ( exp( -hk * xs / (2 * pow(1 + rs, 2)) ) / rs - ( 1 + c * xs * ( 1 + d * xs ) ) );
            }
            bvn = -bvn / twopi;
        }
        if (rho > 0) {
            bvn =  bvn + normp( -max( h, k ) );
        } else {
            bvn = -bvn;
            if (k > h) {
                if (h < 0) {
                    bvn = bvn + normp(k)  - normp(h);
                } else {
                    bvn = bvn + normp(-h) - normp(-k);
                }
            }
        }
    }

    return bvn;
}


 
