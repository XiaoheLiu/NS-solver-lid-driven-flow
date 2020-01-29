#include "VectorField.h"

VectorField::VectorField(double lx, double ly, int ni, int nj) : Lx(lx), Ly(ly), Ni(ni), Nj(nj)
{
  array = new double[(ni + 2) * (nj + 2) * Nk];
  for (int k = 0; k < Nk; k++)
  {
    for (int i = 0; i < ni + 2; i++)
    {
      for (int j = 0; j < nj + 2; j++)
      {
        array[index(k, i, j)] = 0.;
      }
    }
  }
}

VectorField::VectorField(double lx, double ly, int ni, int nj, double (*f[3])(double, double)) : Lx(lx), Ly(ly), Ni(ni), Nj(nj)
{
  array = new double[(ni + 2) * (nj + 2) * Nk];
  set_fields(f);
}

VectorField::~VectorField() { delete[] array; }

VectorField &VectorField::operator+=(const VectorField &dU)
{
  for (int k = 0; k < Nk; k++)
  {
    for (int i = 0; i < Ni + 2; i++)
    {
      for (int j = 0; j < Nj + 2; j++)
      {
        this->array[index(k, i, j)] += dU.at(k, i, j);
      }
    }
  }
  return *this;
}

VectorField &VectorField::operator-=(const VectorField &dU)
{
  for (int k = 0; k < Nk; k++)
  {
    for (int i = 0; i < Ni + 2; i++)
    {
      for (int j = 0; j < Nj + 2; j++)
      {
        this->array[index(k, i, j)] -= dU.at(k, i, j);
      }
    }
  }
  return *this;
}

VectorField &VectorField::operator*=(const double scalar)
{
  for (int k = 0; k < Nk; k++)
  {
    for (int i = 0; i < Ni + 2; i++)
    {
      for (int j = 0; j < Nj + 2; j++)
      {
        this->array[index(k, i, j)] *= scalar;
      }
    }
  }
  return *this;
}

double VectorField::at(int k, int i, int j) const
{
  return array[index(k, i, j)];
}

double VectorField::P(int i, int j) const
{
  return at(0, i, j);
}

double VectorField::u(int i, int j) const
{
  return at(1, i, j);
}

double VectorField::v(int i, int j) const
{
  return at(2, i, j);
}

void VectorField::set(int k, int i, int j, double value) { array[index(k, i, j)] = value; }

void VectorField::set(int i, int j, double values[3])
{
  for (int k = 0; k < Nk; k++)
    array[index(k, i, j)] = values[k];
}

void VectorField::set_fields(double (*f[3])(double, double))
{
  for (int k = 0; k < Nk; k++)
  {
    for (int i = 0; i < Ni + 2; i++)
    {
      for (int j = 0; j < Nj + 2; j++)
      {
        array[index(k, i, j)] = f[k](x(i), y(j));
      }
    }
  }
}

void VectorField::calculate_L2_error(VectorField &Exact)
{
  double rms[3]{0., 0., 0.};
  cout << "- L2 norm:{\n";
  cout.precision(7);
  for (int k = 0; k < Nk; k++)
  {
    for (int i = 1; i < Ni + 1; i++)
    {
      for (int j = 1; j < Nj + 1; j++)
      {
        double diff = (at(k, i, j) - Exact.at(k, i, j));
        rms[k] += diff * diff;
      }
    }
    rms[k] = sqrt(rms[k] / Ni / Nj);
    cout << scientific << rms[k] << endl;
  }
  cout << "}\n\n";
}

double VectorField::calculate_L2_norm(int k)
{
  double rms = 0.0;
  for (int i = 1; i < Ni + 1; i++)
  {
    for (int j = 1; j < Nj + 1; j++)
    {
      double diff = at(k, i, j);
      rms += diff * diff;
    }
  }
  rms = sqrt(rms / Ni / Nj);
  return rms;
}

void VectorField::print(int k)
{
  for (int i = 0; i < Ni + 2; i++)
  {
    for (int j = 0; j < Nj + 2; j++)
    {
      cout << at(k, i, j) << ",";
    }
    cout << "\n";
  }
  cout << "\n";
}

void VectorField::print_partial(int k)
{
  for (int i = 9; i < 12; i++)
  {
    for (int j = 9; j < 12; j++)
    {
      cout << at(k, i, j) << ",";
    }
    cout << "\n";
  }
  cout << "\n";
}

void VectorField::write_dat(string fileName)
{
  ofstream file;
  file.open("./outputs/" + fileName + ".dat");
  if (file.is_open())
  {
    file << setprecision(5);
    for (int i = 1; i < Ni + 1; i++)
    {
      for (int j = 1; j < Nj + 1; j++)
      {
        file << "I:\t" << setw(2) << i << " J:\t" << setw(2) << j;
        file << " P:\t" << setw(9) << at(0, i, j);
        file << " U:\t" << setw(9) << at(1, i, j);
        file << " V:\t" << setw(9) << at(2, i, j) << "\n";
      }
      file << "\n";
    }
  }
  else
    cout << "Unable to open file";
}

void VectorField::write_csv(string fileName)
{
  string cpns[3]{"P", "u", "v"};
  for (int k = 0; k < Nk; k++)
  {
    ofstream file;
    file.open("./outputs/" + fileName + "-" + cpns[k] + ".csv");
    if (file.is_open())
    {
      for (int i = 0; i < Ni + 2; i++)
      {
        for (int j = 0; j < Nj + 1; j++)
        {
          file << at(k, i, j) << ",";
        }
        file << at(k, i, Nj + 1) << "\n";
      }
    }
    else
      cout << "Unable to open file";
  }
}

void VectorField::output_midline_u()
{
  int k = 1;
  ofstream file;
  file.open("./outputs/u_midline-" + to_string(Ni) + "by" + to_string(Nj) + ".csv");
  if (file.is_open())
  {
    for (int j = 0; j < Nj + 2; j++)
    {
      double mid = (at(k, Ni / 2, j) + at(k, Ni / 2 + 1, j)) / 2.;
      file << mid << "\n";
    }
  }
  else
    cout << "Unable to open file";
}

void VectorField::output_vorticity(string fileName)
{
  ofstream file;
  file.open("./outputs/vorticity-" + fileName + ".csv");
  if (file.is_open())
  {
    for (int i = 1; i < Ni + 1; i++)
    {
      for (int j = 1; j < Nj; j++)
      {
        double vorticity = (v(i + 1, j) - v(i - 1, j)) / (2. * Dx) - (u(i, j + 1) - u(i, j - 1)) / (2. * Dy);
        file << vorticity << ",";
      }
      int j = Nj;
      double vorticity = (v(i + 1, j) - v(i - 1, j)) / (2. * Dx) - (u(i, j + 1) - u(i, j - 1)) / (2. * Dy);
      file << vorticity << "\n";
    }
  }
  else
    cout << "Unable to open file";
}

double VectorField::x(int i)
{
  return (double(i) - 0.5) * Dx;
}

double VectorField::y(int j)
{
  return (double(j) - 0.5) * Dy;
}

int VectorField::index(int k, int i, int j) const
{
  return k * (Nj + 2) * (Ni + 2) + ((Nj + 2) * i + j);
}
