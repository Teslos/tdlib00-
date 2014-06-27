#include <fstream.h>
#include <iomanip.h>
#include <cuox_pl.h>

const char *HELP =
"\nEvgenii Rudnyi (C) All rights reserved, 1998 \n\n\
usage:\n\n\
\t\tY123 model.mod p -T z Tmin Tmax Tstep\n\
or\n\
\t\tY123 model.mod p -z T zmin zmax zstep\n\n";

CuOx_plane Y123;

void print_all(const CuOx_plane &Y123, const StateTp &Tp, const double &z)
{
  StateTp Tpr;
  Tpr.T() = 298.15;
  Tpr.p() = Tp.p();
  StateX xx(2);
  xx[0] = 1. - z;
  xx[1] = z;
  double H298 = Y123.Z(::H, ::full, Tpr, xx);
  const vec_double &x_ = Y123.IntVarsEq(Tp, xx);
  cout << setw(8) << Tp.T() << setw(8) << z << setw(8) << x_[0];
  cout << setw(10) << Y123.Z(::Cp, ::full, Tp, xx);
  cout << setw(10) << Y123.Z(::S, ::full, Tp, xx);
  cout << setw(10) << (Y123.Z(::H, ::full, Tp, xx)  - H298)/1000.;
  cout << setw(10) << Y123.Z(::H, ::mix, Tp, xx)/1000.;
  const vec_double &x = Y123.dZdx(::G, ::mix, Tp, xx);
  double tmp = 2.*(x[1] - x[0])/global::R/Tp.T();
  if (tmp > 1e6)
    cout << setw(10)<< "inf";
  else if (tmp < -1e6)
    cout << setw(10) << "-inf";
  else
    cout << setw(10) << tmp;
//  cout << setw(10) << Y123.Z(::V, ::full, Tp, xx)
//       << setw(10) << Y123.Z(::dVdT, ::full, Tp, xx)
//       << setw(10) << Y123.Z(::dVdp, ::full, Tp, xx);
   cout << endl;
}

int main(int argc, char *argv[])
{
  try
  {
    global::neg_log = false;
    cout.precision(2);
    cout.setf(ios::showpoint);
    cout.setf(ios::fixed, ios::floatfield);
    if (argc < 5 || argc > 8)
    {
      cout << HELP;
      return 0;
    }
    ifstream in(argv[1]);
    in >> Y123;
    StateTp Tp;
    Tp.p() = atof(argv[2]);
    cout << "p= " << Tp.p() << endl;
    cout <<
  "    T        z       x        Cpz        S      H - H298  Del_ox_H ln{pO2/po}"
      << endl;
    if (strcmp(argv[3], "-T") == 0)
    {
      double z = 0., Tmin = 0., Tmax = 0., Tstep = 0.;
      z = atof(argv[4]);
      if (argc > 5)
        Tmin = atof(argv[5]);
      if (argc > 6)
        Tmax = atof(argv[6]);
      if (argc > 7)
        Tstep = atof(argv[7]);
      if (Tmin < 100.) Tmin = 100.;
      if (Tmax == 0.)
        Tmax = Tmin;
      else if (Tmax < 100.)
        Tmax = 100.;
      if (Tstep == 0.) Tstep = (Tmax - Tmin)/21.;
      if (Tstep == 0.) Tstep = 1.;
      if ((Tmax - Tmin)*Tstep < 0.) Tstep = -Tstep;
      for (Tp.T() = Tmin;
           Tstep >=0 ? Tp.T() <= Tmax + Tstep/2. : Tp.T() >= Tmax - Tstep/2.;
           Tp.T() += Tstep)
        print_all(Y123, Tp, z);
      return 0;
    }
    else if (strcmp(argv[3], "-z") == 0)
    {
      double T = 0., zmin = 0., zmax = -1., zstep = 0.;
      T = atof(argv[4]);
      if (T < 100)
      {
        cout << HELP;
        return 1;
      }
      if (argc > 5)
        zmin = atof(argv[5]);
      if (argc > 6)
        zmax = atof(argv[6]);
      if (argc > 7)
        zstep = atof(argv[7]);
      if (zmin < 0.) zmin = 0.;
      if (zmax == -1.)
        zmax = 1.;
      else if (zmax < 0.)
        zmax = 0.;
      if (zmax > 1.) zmax = 1.;
      if (zstep == 0.) zstep = (zmax - zmin)/21.;
      if (zstep == 0.) zstep = 1.;
      if ((zmax - zmin)*zstep < 0.) zstep = -zstep;
      Tp.T() = T;
      for (double z = zmin;
           zstep >=0 ? z <= zmax + zstep/2. : z >= zmax - zstep/2.;
           z += zstep)
      {
        if (z > 1.) z = 1.;
        print_all(Y123, Tp, z);
      }
      return 0;
    }
    else
    {
      cout << HELP;
      return 0;
    }
  }
  catch (const gError &mes)
  {
    cout << mes.message << endl;
    return 1;
  }
}
