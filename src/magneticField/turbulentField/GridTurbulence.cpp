#include "crpropa/magneticField/turbulentField/GridTurbulence.h"
#include "crpropa/GridTools.h"
#include "crpropa/Random.h"

#ifdef CRPROPA_HAVE_FFTW3F
#include "fftw3.h"

namespace crpropa {

/* Helper functions for synthetic turbulent field models */

std::vector<std::pair<int, float>> gridPowerSpectrum(ref_ptr<Grid3f> grid) {

  double rms = rmsFieldStrength(grid);
  size_t n = grid->getNx(); // size of array

  // arrays to hold the complex vector components of the B(k)-field
  fftwf_complex *Bkx, *Bky, *Bkz;
  Bkx = (fftwf_complex *)fftwf_malloc(sizeof(fftwf_complex) * n * n * n);
  Bky = (fftwf_complex *)fftwf_malloc(sizeof(fftwf_complex) * n * n * n);
  Bkz = (fftwf_complex *)fftwf_malloc(sizeof(fftwf_complex) * n * n * n);

  fftwf_complex *Bx = (fftwf_complex *)Bkx;
  fftwf_complex *By = (fftwf_complex *)Bky;
  fftwf_complex *Bz = (fftwf_complex *)Bkz;

  // save to temp
  int i;
  for (size_t ix = 0; ix < n; ix++) {
    for (size_t iy = 0; iy < n; iy++) {
      for (size_t iz = 0; iz < n; iz++) {
        i = ix * n * n + iy * n + iz;
        Vector3<float> &b = grid->get(ix, iy, iz);
        Bx[i][0] = b.x / rms;
        By[i][0] = b.y / rms;
        Bz[i][0] = b.z / rms;
      }
    }
  }

  // in-place, real to complex, inverse Fourier transformation on each component
  // note that the last elements of B(x) are unused now
  fftwf_plan plan_x =
      fftwf_plan_dft_3d(n, n, n, Bx, Bkx, FFTW_FORWARD, FFTW_ESTIMATE);
  fftwf_execute(plan_x);
  fftwf_destroy_plan(plan_x);

  fftwf_plan plan_y =
      fftwf_plan_dft_3d(n, n, n, By, Bky, FFTW_FORWARD, FFTW_ESTIMATE);
  fftwf_execute(plan_y);
  fftwf_destroy_plan(plan_y);

  fftwf_plan plan_z =
      fftwf_plan_dft_3d(n, n, n, Bz, Bkz, FFTW_FORWARD, FFTW_ESTIMATE);
  fftwf_execute(plan_z);
  fftwf_destroy_plan(plan_z);

  float power;
  std::map<size_t, std::pair<float, int>> spectrum;
  int k;

  for (size_t ix = 0; ix < n; ix++) {
    for (size_t iy = 0; iy < n; iy++) {
      for (size_t iz = 0; iz < n; iz++) {
        i = ix * n * n + iy * n + iz;
        k = static_cast<int>(
            std::floor(std::sqrt(ix * ix + iy * iy + iz * iz)));
        if (k > n / 2. || k == 0)
          continue;
        power = ((Bkx[i][0] * Bkx[i][0] + Bkx[i][1] * Bkx[i][1]) +
                 (Bky[i][0] * Bky[i][0] + Bky[i][1] * Bky[i][1]) +
                 (Bkz[i][0] * Bkz[i][0] + Bkz[i][1] * Bkz[i][1]));
        if (spectrum.find(k) == spectrum.end()) {
          spectrum[k].first = power;
          spectrum[k].second = 1;
        } else {
          spectrum[k].first += power;
          spectrum[k].second += 1;
        }
      }
    }
  }

  fftwf_free(Bkx);
  fftwf_free(Bky);
  fftwf_free(Bkz);

  std::vector<std::pair<int, float>> points;
  for (std::map<size_t, std::pair<float, int>>::iterator it = spectrum.begin();
       it != spectrum.end(); ++it) {
    points.push_back(
        std::make_pair(it->first, (it->second).first / (it->second).second));
  }

  return points;
}

// Check the grid properties before the FFT procedure
void checkGridRequirementsTEMP(ref_ptr<Grid3f> grid, double lMin, double lMax) {
  size_t Nx = grid->getNx();
  size_t Ny = grid->getNy();
  size_t Nz = grid->getNz();
  Vector3d spacing = grid->getSpacing();

  if ((Nx != Ny) or (Ny != Nz))
    throw std::runtime_error("turbulentField: only cubic grid supported");
  if ((spacing.x != spacing.y) or (spacing.y != spacing.z))
    throw std::runtime_error("turbulentField: only equal spacing suported");
  if (lMin < 2 * spacing.x)
    throw std::runtime_error("turbulentField: lMin < 2 * spacing");
  if (lMin >= lMax)
    throw std::runtime_error("turbulentField: lMin >= lMax");
  if (lMax > Nx * spacing.x) // before was (lMax > Nx * spacing.x / 2)
    throw std::runtime_error("turbulentField: lMax > size");
}

// Execute inverse discrete FFT in-place for a 3D grid, from complex to real
// space
void executeInverseFFTInplaceTEMP(ref_ptr<Grid3f> grid, fftwf_complex *Bkx,
                                  fftwf_complex *Bky, fftwf_complex *Bkz) {

  size_t n = grid->getNx(); // size of array
  size_t n2 = (size_t)floor(n / 2) +
              1; // size array in z-direction in configuration space

  // in-place, complex to real, inverse Fourier transformation on each component
  // note that the last elements of B(x) are unused now
  float *Bx = (float *)Bkx;
  fftwf_plan plan_x = fftwf_plan_dft_c2r_3d(n, n, n, Bkx, Bx, FFTW_ESTIMATE);
  fftwf_execute(plan_x);
  fftwf_destroy_plan(plan_x);

  float *By = (float *)Bky;
  fftwf_plan plan_y = fftwf_plan_dft_c2r_3d(n, n, n, Bky, By, FFTW_ESTIMATE);
  fftwf_execute(plan_y);
  fftwf_destroy_plan(plan_y);

  float *Bz = (float *)Bkz;
  fftwf_plan plan_z = fftwf_plan_dft_c2r_3d(n, n, n, Bkz, Bz, FFTW_ESTIMATE);
  fftwf_execute(plan_z);
  fftwf_destroy_plan(plan_z);

  // save to grid
  int i;
  for (size_t ix = 0; ix < n; ix++) {
    for (size_t iy = 0; iy < n; iy++) {
      for (size_t iz = 0; iz < n; iz++) {
        i = ix * n * 2 * n2 + iy * 2 * n2 + iz;
        Vector3f &b = grid->get(ix, iy, iz);
        b.x = Bx[i];
        b.y = By[i];
        b.z = Bz[i];
      }
    }
  }
}


GridTurbulence::GridTurbulence(double Brms, double sindex, double qindex,
                               double lBendover, double lMin, double lMax,
                               int gridSize, double boxSize, unsigned int seed)
    : TurbulentField(Brms, sindex, qindex, lBendover), lMin(lMin), lMax(lMax),
      gridSize(gridSize), boxSize(boxSize), seed(seed),
      spacing(boxSize / gridSize) {
  initGrid();
  initTurbulence(gridPtr, Brms, lMin, lMax, -sindex - 2, seed);
}

Vector3d GridTurbulence::getField(const Vector3d &pos) const {
  return gridPtr->interpolate(pos);
}

void GridTurbulence::initGrid() {
  gridPtr = new Grid3f(Vector3d(-boxSize / 2), gridSize, spacing);
}

void GridTurbulence::initTurbulence(ref_ptr<Grid3f> grid, double Brms,
                                    double lMin, double lMax, double alpha,
                                    int seed) {

  checkGridRequirementsTEMP(grid, lMin, lMax);

  Vector3d spacing = grid->getSpacing();
  size_t n = grid->getNx(); // size of array
  size_t n2 = (size_t)floor(n / 2) +
              1; // size array in z-direction in configuration space

  // arrays to hold the complex vector components of the B(k)-field
  fftwf_complex *Bkx, *Bky, *Bkz;
  Bkx = (fftwf_complex *)fftwf_malloc(sizeof(fftwf_complex) * n * n * n2);
  Bky = (fftwf_complex *)fftwf_malloc(sizeof(fftwf_complex) * n * n * n2);
  Bkz = (fftwf_complex *)fftwf_malloc(sizeof(fftwf_complex) * n * n * n2);

  Random random;
  if (seed != 0)
    random.seed(seed); // use given seed

  // calculate the n possible discrete wave numbers
  double K[n];
  for (int i = 0; i < n; i++)
    K[i] = ((double)i / n - i / (n / 2));

  // construct the field in configuration space
  int i;
  double k;

  // parameters goes for non-helical calculations
  double theta, phase, cosPhase, sinPhase, amplitude;

  // double kMin = 2*M_PI / lMax; // * 2 * spacing.x; // spacing.x / lMax;
  // double kMax = 2*M_PI / lMin; // * 2 * spacing.x; // spacing.x / lMin;
  double kMin = spacing.x / lMax;
  double kMax = spacing.x / lMin;
  auto lambda = lBendover / spacing.x * 2 * M_PI;

  Vector3f b;           // real b-field vector
  Vector3f ek, e1, e2;  // orthogonal base
  Vector3f n0(1, 1, 1); // arbitrary vector to construct orthogonal base

  for (size_t ix = 0; ix < n; ix++) {
    for (size_t iy = 0; iy < n; iy++) {
      for (size_t iz = 0; iz < n2; iz++) {

        i = ix * n * n2 + iy * n2 + iz;
        ek.setXYZ(K[ix], K[iy], K[iz]);
        k = ek.getR();

        // wave outside of turbulent range -> B(k) = 0
        if ((k < kMin) || (k > kMax)) {
          Bkx[i][0] = 0;
          Bkx[i][1] = 0;
          Bky[i][0] = 0;
          Bky[i][1] = 0;
          Bkz[i][0] = 0;
          Bkz[i][1] = 0;
          continue;
        }

        // construct an orthogonal base ek, e1, e2
        if (ek.isParallelTo(n0, float(1e-3))) {
          // ek parallel to (1,1,1)
          e1.setXYZ(-1., 1., 0);
          e2.setXYZ(1., 1., -2.);
        } else {
          // ek not parallel to (1,1,1)
          e1 = n0.cross(ek);
          e2 = ek.cross(e1);
        }
        e1 /= e1.getR();
        e2 /= e2.getR();

        // random orientation perpendicular to k
        theta = 2 * M_PI * random.rand();
        b = e1 * cos(theta) + e2 * sin(theta);

        // normal distributed amplitude with mean = 0 and sigma = k^alpha/2
        amplitude = k * lambda *
                    pow(1.0 + k * k * lambda * lambda, (alpha / 2 - 1) / 2.);
        b *= amplitude;

        // uniform random phase
        phase = 2 * M_PI * random.rand();
        cosPhase = cos(phase); // real part
        sinPhase = sin(phase); // imaginary part

        Bkx[i][0] = b.x * cosPhase;
        Bkx[i][1] = b.x * sinPhase;
        Bky[i][0] = b.y * cosPhase;
        Bky[i][1] = b.y * sinPhase;
        Bkz[i][0] = b.z * cosPhase;
        Bkz[i][1] = b.z * sinPhase;
      } // for iz
    }   // for iy
  }     // for ix

  executeInverseFFTInplaceTEMP(grid, Bkx, Bky, Bkz);

  fftwf_free(Bkx);
  fftwf_free(Bky);
  fftwf_free(Bkz);

  scaleGrid(grid, Brms / rmsFieldStrength(grid)); // normalize to Brms
}

} // namespace crpropa

#endif // CRPROPA_HAVE_FFTW3F
