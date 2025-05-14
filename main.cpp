#include "Connectivity.hpp"
#include "Gmsh.hpp"
#include "Nodes.hpp"
#include "Numerics.hpp"
#include "Tools.hpp"

int main(int argc, char **argv) {

  // ======================================================================
  // 1A) Get simulation parameters
  // ======================================================================

  string meshFile;               // Mesh file
  int N;                         // Degree of polynomial bases
  int iOutputGmsh;               // Number of time steps between gmsh output
  double FinalTime;              // Final time of the simulation
  vector<int> _CoefDataBaseNum;  // Database for parameters
  vector<double> _CoefDataBase1; // Database for parameters
  vector<double> _CoefDataBase2; // Database for parameters

  // Read setup file
  ifstream setupfile("setup");
  if (!setupfile.is_open()) {
    cerr << "ERROR - Setup file not opened.\n";
    exit(1);
  }
  string line;
  while (setupfile) {
    setupfile >> line;
    if (line == "$Mesh")
      setupfile >> meshFile;
    if (line == "$PolynomialDegree")
      setupfile >> N;
    if (line == "$OutputStep")
      setupfile >> iOutputGmsh;
    if (line == "$Duration")
      setupfile >> FinalTime;
    if (line == "$Param") {
      int PARAM;
      setupfile >> PARAM;
      _CoefDataBaseNum.resize(PARAM);
      _CoefDataBase1.resize(PARAM);
      _CoefDataBase2.resize(PARAM);
      for (int i = 0; i < PARAM; i++) {
        setupfile >> _CoefDataBaseNum[i];
        setupfile >> _CoefDataBase1[i];
        setupfile >> _CoefDataBase2[i];
      }
    }
  }

  // ======================================================================
  // 1B) Load mesh
  // ======================================================================

  int K;              // Number of elements in the mesh
  vector<double> _VX; // Coordinate 'x' of the vertices of the mesh [#vert]
  vector<double> _VY; // Coordinate 'y' of the vertices of the mesh [#vert]
  vector<double> _VZ; // Coordinate 'z' of the vertices of the mesh [#vert]
  vector<int> _EMsh;  // List of orginal gmsh numbers of elements [K]
  vector<int> _ETag;  // List of element tag (gmsh: physical tag) [K]
  vector<int> _EToV;  // Element-to-vertex connectivity array [K,NVertTet]

  // Load mesh
  loadMeshGmsh(meshFile, K, _VX, _VY, _VZ, _EMsh, _ETag, _EToV);

  // ======================================================================
  // 1C) Get finite element nodes
  // ======================================================================

  int Np = (N + 1) * (N + 2) * (N + 3) / 6; // Number of nodes per element
  int Nfp = (N + 1) * (N + 2) / 2;          // Number of nodes per face
  vector<double> _r, _s, _t;                // Coordinates of reference nodes [Np]
  vector<double> _x, _y, _z;                // Coordinates of physical nodes [K,Np]

  // Get local coordinates of nodes on the reference element
  getReferenceNodes(N, _r, _s, _t);

  // Get global coordinates of nodes on the physical elements
  getPhysicalNodes(_EToV, _VX, _VY, _VZ, _r, _s, _t, _x, _y, _z);

  // ======================================================================
  // 1D) Build connectivity matrices, elemental matrices & geometric factors
  // ======================================================================

  vector<int> _EToE;  // Element-to-Element [K,NFacesTet]
  vector<int> _EToF;  // Element-to-Face [K,NFacesTet]
  vector<int> _NfToN; // LocalFaceNode-to-LocalNode [NFacesTet,Nfp]
  vector<int> _mapP;  // GlobalFaceNode-to-NeighborGlobalNode [K,Nfp,NFacesTet]

  // Build connectivity matrices
  buildConnectivity(K, Nfp, Np, _r, _s, _t, _x, _y, _z, _EToV, _EToE, _EToF, _NfToN, _mapP);

  vector<double> _Dr;   // r-derivative matrix [Np,Np]
  vector<double> _Ds;   // s-derivative matrix [Np,Np]
  vector<double> _Dt;   // t-derivative matrix [Np,Np]
  vector<double> _LIFT; // Lift matrix [NFacesTet,Nfp,Np]

  // Build elemental matrices
  buildElemMatrices(N, _r, _s, _t, _NfToN, _Dr, _Ds, _Dt, _LIFT);

  vector<double> _rstxyz; // ... for volume terms  [K,9]
  vector<double> _Fscale; // ... for surface terms [K,NFacesTet]
  vector<double> _nx;     // x-component of normal to the face [K,NFacesTet]
  vector<double> _ny;     // y-component of normal to the face [K,NFacesTet]
  vector<double> _nz;     // z-component of normal to the face [K,NFacesTet]

  // Build geometric factors
  buildGeomFactors(_EToV, _VX, _VY, _VZ, _rstxyz, _Fscale, _nx, _ny, _nz);

  // ======================================================================
  // 1E) Build physical coefficient maps
  // ======================================================================

  vector<double> _c(K);
  vector<double> _rho(K);

  for (int k = 0; k < K; k++) {
    int i = 0;
    while (i < _CoefDataBaseNum.size()) {
      if (_CoefDataBaseNum[i] == _ETag[k]) {
        _c[k] = _CoefDataBase1[i];
        _rho[k] = _CoefDataBase2[i];
        break;
      }
      i++;
    }
    if (i == _CoefDataBaseNum.size()) {
      cerr << "ERROR - ETag " << _ETag[k] << " not found in parameter database\n";
      exit(1);
    }
  }
  exportParamGmsh("info_c", _EMsh, _c);
  exportParamGmsh("info_rho", _EMsh, _rho);

  // ======================================================================
  // 1F) Build time stepping scheme
  // ======================================================================

  double dt = 1e9;
  for (int k = 0; k < K; k++) {
    for (int f = 0; f < NFacesTet; ++f) {
      double val = 1. / (_c[k] * (N + 1) * (N + 1) * _Fscale[k * NFacesTet + f]);
      if (val < dt)
        dt = val;
    }
  }
  double CFL = 0.75;                 // CFL
  //double CFL = 1;
  dt *= CFL;                         // Time step
  int Nsteps = ceil(FinalTime / dt); // Number of global time steps

  vector<double> rk4a, rk4b, rk4c;
  getRungeKuttaCoefficients(rk4a, rk4b, rk4c);

  cout << "INFO - DT = " << dt << " and Nsteps = " << Nsteps << "\n";

  // ======================================================================
  // 1G) Memory storage for fields, RHS and residual at nodes
  // ======================================================================

  int Nfields = 4; // Number of unknown fields

  // Memory storage at each node of each element
  vector<double> _valQ(K * Np * Nfields, 0.); // Values of fields
  vector<double> _rhsQ(K * Np * Nfields, 0.); // RHS
  vector<double> _resQ(K * Np * Nfields, 0.); // residual

  // Initialization
  for (int k = 0; k < K; ++k) {
    for (int n = 0; n < Np; ++n) {
      double x = _x[k * Np + n];
      double y = _y[k * Np + n];
      double z = _z[k * Np + n];
      _valQ[k * Np * Nfields + n * Nfields + 0] = exp(-(x * x + y * y + z * z) / 0.1);
    }
  }

  cout << "INFO - Total number of DOF = " << K * Np * Nfields << "\n";

  // Export initial solution
  if (iOutputGmsh > 0)
    exportSolGmsh(N, _r, _s, _t, _EMsh, 0, 0., _valQ);

  // ======================================================================
  // 2) RUN
  // ======================================================================


  double time_total_begin = dsecnd();
  double time = 0;
  double average = 0;
  // Global time iteration
  for (int nGlo = 0; nGlo < Nsteps; ++nGlo) {
    double time_begin = dsecnd();
    double runTime = nGlo * dt; // Time at the beginning of the step

    // Local time iteration
    for (int nLoc = 0; nLoc < 5; ++nLoc) {
      double a = rk4a[nLoc];
      double b = rk4b[nLoc];

      // ======================== (1) UPDATE RHS

      for (int k = 0; k < K; ++k) {

        double c = _c[k];
        double rho = _rho[k];

        // ======================== (1.1) UPDATE PENALTY VECTOR

        vector<double> s_p_flux(Nfp * NFacesTet);
        vector<double> s_u_flux(Nfp * NFacesTet);
        vector<double> s_v_flux(Nfp * NFacesTet);
        vector<double> s_w_flux(Nfp * NFacesTet);

        for (int f = 0; f < NFacesTet; f++) {

          // Fetch normal
          double nx = _nx[k * NFacesTet + f];
          double ny = _ny[k * NFacesTet + f];
          double nz = _nz[k * NFacesTet + f];

          // Fetch medium parameters
          int k2 = _EToE[k * NFacesTet + f];
          double cP = _c[k2];
          double rhoP = _rho[k2];
          
          double rho_c = rho * c;
          double un_sur_rhoPcP = 1. / (rhoP * cP);
          double rhoPcP_rhoc = rhoP * cP + rho_c;

          // Compute penalty terms
          for (int nf = f * Nfp; nf < (f + 1) * Nfp; nf++) {

            // Index of node in current element
            int n1 = _NfToN[nf];

            // Index of node in neighbor element
            int n2 = _mapP[nf + k * NFacesTet * Nfp];
            int coeff1 = k * Np * Nfields + n1 * Nfields;

            // Load values 'minus' corresponding to current element
            double pM = _valQ[coeff1 + 0];
            double uM = _valQ[coeff1 + 1];
            double vM = _valQ[coeff1 + 2];
            double wM = _valQ[coeff1 + 3];
            double nMdotuM = (nx * uM + ny * vM + nz * wM);

            if (n2 >= 0) { // ... if there is a neighbor element ...

              int coeff2 = k2 * Np * Nfields + n2 * Nfields;
              // Load values 'plus' corresponding to neighbor element
              double pP = _valQ[coeff2 + 0];
              double uP = _valQ[coeff2 + 1];
              double vP = _valQ[coeff2 + 2];
              double wP = _valQ[coeff2 + 3];
              double nMdotuP = (nx * uP + ny * vP + nz * wP);

              // Penalty terms for interface between two elements
              s_p_flux[nf] = c / (un_sur_rhoPcP + 1. / (rho_c)) * ((nMdotuP - nMdotuM) - un_sur_rhoPcP * (pP - pM));
              s_u_flux[nf] = nx * c / (rhoPcP_rhoc) * ((pP - pM) - rhoP * cP * (nMdotuP - nMdotuM));
              s_v_flux[nf] = ny * c / (rhoPcP_rhoc) * ((pP - pM) - rhoP * cP * (nMdotuP - nMdotuM));
              s_w_flux[nf] = nz * c / (rhoPcP_rhoc) * ((pP - pM) - rhoP * cP * (nMdotuP - nMdotuM));

            } else {

              // Homogeneous Dirichlet on 'p'
              double tmp = -2. / (rhoP * cP) * pM;
              // Homogeneous Dirichlet on 'u'
              // double tmp = 2*nMdotuM;
              // ABC
              // double tmp = nMdotuM - 1./(rhoP*cP) * pM;

              // Penalty terms for boundary of the domain
              s_p_flux[nf] = -c / (1. / (rhoP * cP) + 1. / (rho_c)) * tmp;
              s_u_flux[nf] = nx * c / (rhoPcP_rhoc) * tmp;
              s_v_flux[nf] = ny * c / (rhoPcP_rhoc) * tmp;
              s_w_flux[nf] = nz * c / (rhoPcP_rhoc) * tmp;
            }
          }
        }

        // ======================== (1.2) COMPUTING VOLUME TERMS

        // Load geometric factors
        double rx = _rstxyz[k * 9 + 0];
        double ry = _rstxyz[k * 9 + 1];
        double rz = _rstxyz[k * 9 + 2];
        double sx = _rstxyz[k * 9 + 3];
        double sy = _rstxyz[k * 9 + 4];
        double sz = _rstxyz[k * 9 + 5];
        double tx = _rstxyz[k * 9 + 6];
        double ty = _rstxyz[k * 9 + 7];
        double tz = _rstxyz[k * 9 + 8];

        // Load fields
        vector<double> s_p(Np);
        vector<double> s_u(Np);
        vector<double> s_v(Np);
        vector<double> s_w(Np);
        for (int n = 0; n < Np; ++n) {
          int coeff3 = k * Np * Nfields + n * Nfields;
          s_p[n] = _valQ[coeff3 + 0];
          s_u[n] = _valQ[coeff3 + 1];
          s_v[n] = _valQ[coeff3 + 2];
          s_w[n] = _valQ[coeff3 + 3];
        }
        // Compute mat-vec product for surface term
        for (int n = 0; n < Np; ++n) {
          double dpdr = 0, dpds = 0, dpdt = 0;
          double dudr = 0, duds = 0, dudt = 0;
          double dvdr = 0, dvds = 0, dvdt = 0;
          double dwdr = 0, dwds = 0, dwdt = 0;
          for (int m = 0; m < Np; ++m) {
            double Dr = _Dr[n + m * Np];
            dpdr += Dr * s_p[m];
            dudr += Dr * s_u[m];
            dvdr += Dr * s_v[m];
            dwdr += Dr * s_w[m];
            double Ds = _Ds[n + m * Np];
            dpds += Ds * s_p[m];
            duds += Ds * s_u[m];
            dvds += Ds * s_v[m];
            dwds += Ds * s_w[m];
            double Dt = _Dt[n + m * Np];
            dpdt += Dt * s_p[m];
            dudt += Dt * s_u[m];
            dvdt += Dt * s_v[m];
            dwdt += Dt * s_w[m];
          }

          double dpdx = rx * dpdr + sx * dpds + tx * dpdt;
          double dpdy = ry * dpdr + sy * dpds + ty * dpdt;
          double dpdz = rz * dpdr + sz * dpds + tz * dpdt;
          double dudx = rx * dudr + sx * duds + tx * dudt;
          double dvdy = ry * dvdr + sy * dvds + ty * dvdt;
          double dwdz = rz * dwdr + sz * dwds + tz * dwdt;
          double divU = dudx + dvdy + dwdz;
          // Compute RHS (only part corresponding to volume terms)
          int coeff4 = k * Np * Nfields + n * Nfields;
          _rhsQ[coeff4 + 0] = -c * c * rho * divU;
          _rhsQ[coeff4 + 1] = -1. / rho * dpdx;
          _rhsQ[coeff4 + 2] = -1. / rho * dpdy;
          _rhsQ[coeff4 + 3] = -1. / rho * dpdz;
        }

        // ======================== (1.3) COMPUTING SURFACE TERMS

        for (int n = 0; n < Np; ++n) {
          for (int f = 0; f < NFacesTet; f++) {

            // Compute mat-vec product for surface term
            double p_lift = 0.;
            double u_lift = 0.;
            double v_lift = 0.;
            double w_lift = 0.;
            for (int m = f * Nfp; m < (f + 1) * Nfp; m++) {
              double tmp = _LIFT[n * NFacesTet * Nfp + m];
              p_lift += tmp * s_p_flux[m];
              u_lift += tmp * s_u_flux[m];
              v_lift += tmp * s_v_flux[m];
              w_lift += tmp * s_w_flux[m];
            }

            // Load geometric factor
            double Fscale = _Fscale[k * NFacesTet + f];

            // Update RHS (with part corresponding to surface terms)
            int coeff5 = k * Np * Nfields + n * Nfields;
            _rhsQ[coeff5 + 0] -= p_lift * Fscale;
            _rhsQ[coeff5 + 1] -= u_lift * Fscale;
            _rhsQ[coeff5 + 2] -= v_lift * Fscale;
            _rhsQ[coeff5 + 3] -= w_lift * Fscale;
          }
        }
      }

      // ======================== (2) UPDATE RESIDUAL + FIELDS

      for (int k = 0; k < K; ++k) {
        for (int n = 0; n < Np; ++n) {
          for (int iField = 0; iField < Nfields; ++iField) {
            int id = k * Np * Nfields + n * Nfields + iField;
            _resQ[id] = a * _resQ[id] + dt * _rhsQ[id];
            _valQ[id] = _valQ[id] + b * _resQ[id];
          }
        }
      }
    }


    double time_end = dsecnd();
    average += time_end - time_begin;
    //cout << "Temps d'exécution maximal: " << time << "secondes. Itérations:" << nGlo << endl;
    // Export solution
    if ((iOutputGmsh > 0) && ((nGlo + 1) % iOutputGmsh == 0))
      exportSolGmsh(N, _r, _s, _t, _EMsh, nGlo + 1, runTime + dt, _valQ);
      
  }


  // Temps d'exécution moyen par itérations, temps total d'exécution
  double time_total_end = dsecnd();
  double time_total = time_total_end - time_total_begin;
  average /= Nsteps;
  //double alpha = Nsteps*5*K* ((NFacesTet+1)*Nfp*80 + Np*Np*60 + Np*NFacesTet*8*(1 + (NFacesTet+1)*Nfp) + Np*Nfields*5);
  double flops_par_iter = 1600.0 + 68.0*Np;
  double total_flops = 5.0 * Nsteps * K * flops_par_iter;
  double gflops = total_flops/(time_total*1e9);
  cout << "Temps d'exécution moyen d'une itération: " << average << " secondes." << endl;
  cout << "Temps d'exécution total: " << time_total << " secondes." << endl;
  cout << "Performance arithmétique: " << gflops << " GFLOP/s." << endl;;

  // ======================================================================
  // 3) Post-processing
  // ======================================================================*


  
  cout << "Degré des Polynômes: " << N << endl;
  cout << "OutputStep: " << iOutputGmsh << endl;
  // Export final solution
  if (iOutputGmsh > 0)
    exportSolGmsh(N, _r, _s, _t, _EMsh, Nsteps, Nsteps * dt, _valQ);

  return 0;
}
