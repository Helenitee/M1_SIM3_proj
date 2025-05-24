#include "Connectivity.hpp"
#include "Gmsh.hpp"
#include "Nodes.hpp"
#include "Numerics.hpp"
#include "Tools.hpp"
#include <iostream>
#include <iomanip>
#include <omp.h>
#include <mkl.h>

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
  int totalDOF =  K * Np * Nfields;

  cout << "INFO - Total number of DOF = " << totalDOF << "\n";

  // Export initial solution
  if (iOutputGmsh > 0)
    exportSolGmsh(N, _r, _s, _t, _EMsh, 0, 0., _valQ);

  // ======================================================================
  // 2) RUN
  // ======================================================================
  omp_set_num_threads(4);


  // Define quantities for runtime and average runtime
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
      //#pragma omp parallel for
      #pragma omp parallel for schedule(static)
      for (int k = 0; k < K; ++k) {

        int size_flux_vect = Nfp * NFacesTet; // we compute it once!
        vector<double> s_p_flux(size_flux_vect);
        vector<double> s_u_flux(size_flux_vect);
        vector<double> s_v_flux(size_flux_vect);
        vector<double> s_w_flux(size_flux_vect);


        double c = _c[k];
        double rho = _rho[k];

        // ======================== (1.1) UPDATE PENALTY VECTOR

        for (int f = 0; f < NFacesTet; f++) {

          // Fetch normal
          int comp_normal_indx = k * NFacesTet + f;
          double nx = _nx[comp_normal_indx];
          double ny = _ny[comp_normal_indx];
          double nz = _nz[comp_normal_indx];

          // Fetch medium parameters
          int k2 = _EToE[comp_normal_indx];
          double cP = _c[k2];
          double rhoP = _rho[k2];
          

          // Compute penalty terms
          for (int nf = f * Nfp; nf < (f + 1) * Nfp; nf++) {

            // Index of node in current element
            int n1 = _NfToN[nf];

            // Index of node in neighbor element
            int n2 = _mapP[nf + k * NFacesTet * Nfp];

            // Load values 'minus' corresponding to current element
            int elem_k = k * Np * Nfields + n1 * Nfields;
            double pM = _valQ[elem_k + 0];
            double uM = _valQ[elem_k + 1];
            double vM = _valQ[elem_k + 2];
            double wM = _valQ[elem_k + 3];
            double nMdotuM = (nx * uM + ny * vM + nz * wM);
            // NbOp += 5;

            double coeff_rho_c = rhoP * cP + rho * c;
            if (n2 >= 0) { // ... if there is a neighbor element ...

              // Load values 'plus' corresponding to neighbor element
              int elem_k2 = k2 * Np * Nfields + n2 * Nfields;
              double pP = _valQ[elem_k2 + 0];
              double uP = _valQ[elem_k2 + 1];
              double vP = _valQ[elem_k2 + 2];
              double wP = _valQ[elem_k2 + 3];
              double nMdotuP = (nx * uP + ny * vP + nz * wP);
              // NbOp += 5;

              // Penalty terms for interface between two elements
              double penalty_term_int = c / coeff_rho_c * ((pP - pM) - rhoP * cP * (nMdotuP - nMdotuM));
              
              s_p_flux[nf] = c / (1. / (rhoP * cP) + 1. / (rho * c)) * ((nMdotuP - nMdotuM) - 1. / (rhoP * cP) * (pP - pM));
              s_u_flux[nf] = nx * penalty_term_int;
              s_v_flux[nf] = ny * penalty_term_int;
              s_w_flux[nf] = nz * penalty_term_int;
              // NbOp += (13+11*3);
            } else { // conditions limites 

              // Homogeneous Dirichlet on 'p'
              double tmp = -2. / (rhoP * cP) * pM;
              // NbOp += 3;
              // Homogeneous Dirichlet on 'u'
              // double tmp = 2*nMdotuM;
              // ABC
              // double tmp = nMdotuM - 1./(rhoP*cP) * pM;

              // Penalty terms for boundary of the domain
              double penalty_term_ext = c / coeff_rho_c * tmp;
              s_p_flux[nf] = -c / (1. / (rhoP * cP) + 1. / (rho * c)) * tmp;
              s_u_flux[nf] = nx * penalty_term_ext;
              s_v_flux[nf] = ny * penalty_term_ext;
              s_w_flux[nf] = nz * penalty_term_ext;
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

        double ck = _c[k];
        double rhok = _rho[k];

        const double* valQk = &_valQ[k * Np * Nfields];

        // Vecteurs corrigés pour être à la suite (car valQk ne stocke pas les champs de la bonne manière pour mkl)
        vector<double> q_p(Np), q_u(Np), q_v(Np), q_w(Np);
        for (int n = 0; n < Np; ++n) {
          q_p[n] = valQk[n * Nfields + 0];
          q_u[n] = valQk[n * Nfields + 1];
          q_v[n] = valQk[n * Nfields + 2];
          q_w[n] = valQk[n * Nfields + 3];
        }

        
        // Vecteurs temporaires pour les dérivées
        vector<double> dpdr(Np), dpds(Np), dpdt(Np);
        vector<double> dudr(Np), duds(Np), dudt(Np);
        vector<double> dvdr(Np), dvds(Np), dvdt(Np);
        vector<double> dwdr(Np), dwds(Np), dwdt(Np);


        // Produits Matrice-Vecteur pour les termes de surface
        // 
        // Dr * Q = dérivées en r
        cblas_dgemv(CblasColMajor, CblasNoTrans, Np, Np, 1.0, _Dr.data(), Np, q_p.data(), 1, 0.0, dpdr.data(), 1);
        cblas_dgemv(CblasColMajor, CblasNoTrans, Np, Np, 1.0, _Dr.data(), Np, q_u.data(), 1, 0.0, dudr.data(), 1);
        cblas_dgemv(CblasColMajor, CblasNoTrans, Np, Np, 1.0, _Dr.data(), Np, q_v.data(), 1, 0.0, dvdr.data(), 1);
        cblas_dgemv(CblasColMajor, CblasNoTrans, Np, Np, 1.0, _Dr.data(), Np, q_w.data(), 1, 0.0, dwdr.data(), 1);

        // Ds * Q = dérivées en s
        cblas_dgemv(CblasColMajor, CblasNoTrans, Np, Np, 1.0, _Ds.data(), Np, q_p.data(), 1, 0.0, dpds.data(), 1);
        cblas_dgemv(CblasColMajor, CblasNoTrans, Np, Np, 1.0, _Ds.data(), Np, q_u.data(), 1, 0.0, duds.data(), 1);
        cblas_dgemv(CblasColMajor, CblasNoTrans, Np, Np, 1.0, _Ds.data(), Np, q_v.data(), 1, 0.0, dvds.data(), 1);
        cblas_dgemv(CblasColMajor, CblasNoTrans, Np, Np, 1.0, _Ds.data(), Np, q_w.data(), 1, 0.0, dwds.data(), 1);

        // Dt * Q = dérivées en t
        cblas_dgemv(CblasColMajor, CblasNoTrans, Np, Np, 1.0, _Dt.data(), Np, q_p.data(), 1, 0.0, dpdt.data(), 1);
        cblas_dgemv(CblasColMajor, CblasNoTrans, Np, Np, 1.0, _Dt.data(), Np, q_u.data(), 1, 0.0, dudt.data(), 1);
        cblas_dgemv(CblasColMajor, CblasNoTrans, Np, Np, 1.0, _Dt.data(), Np, q_v.data(), 1, 0.0, dvdt.data(), 1);
        cblas_dgemv(CblasColMajor, CblasNoTrans, Np, Np, 1.0, _Dt.data(), Np, q_w.data(), 1, 0.0, dwdt.data(), 1);

        // Construction des RHS à partir des dérivées
        for (int n = 0; n < Np; ++n) {
            double dpdx = rx * dpdr[n] + sx * dpds[n] + tx * dpdt[n];
            double dpdy = ry * dpdr[n] + sy * dpds[n] + ty * dpdt[n];
            double dpdz = rz * dpdr[n] + sz * dpds[n] + tz * dpdt[n];

            double dudx = rx * dudr[n] + sx * duds[n] + tx * dudt[n];
            double dvdy = ry * dvdr[n] + sy * dvds[n] + ty * dvdt[n];
            double dwdz = rz * dwdr[n] + sz * dwds[n] + tz * dwdt[n];

            double divU = dudx + dvdy + dwdz;

            int id = k * Np * Nfields + n * Nfields;

            _rhsQ[id + 0] = -ck * ck * rhok * divU;
            _rhsQ[id + 1] = -1. / rhok * dpdx;
            _rhsQ[id + 2] = -1. / rhok * dpdy;
            _rhsQ[id + 3] = -1. / rhok * dpdz;
        }

        // ======================== (1.3) COMPUTING SURFACE TERM       
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
            int elem_k = k * Np * Nfields + n * Nfields;
            _rhsQ[elem_k + 0] -= p_lift * Fscale;
            _rhsQ[elem_k + 1] -= u_lift * Fscale;
            _rhsQ[elem_k + 2] -= v_lift * Fscale;
            _rhsQ[elem_k + 3] -= w_lift * Fscale;
          }
        }
      }
      

      // ======================== (2) UPDATE RESIDUAL + FIELDS
      // _resQ[id] = a * _resQ[id]
      cblas_dscal(totalDOF, a, _resQ.data(), 1);
      // resQ += dt * rhsQ
      cblas_daxpy(totalDOF, dt, _rhsQ.data(), 1, _resQ.data(), 1);
      // valQ += b * resQ
      cblas_daxpy(totalDOF, b, _resQ.data(), 1, _valQ.data(), 1);
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

  // Performance arithmétique
  double total_flops = 5.0 * K * Nsteps * (20 * Nfp + Np * (24 * Np + 53 + 8*Nfp + 5*Nfields));
  double gflops = total_flops/(time_total*1e9);

  // Interface
  cout << "Temps d'exécution moyen d'une itération: " << average << " secondes." << endl;
  cout << "Temps d'exécution total: " << time_total << " secondes." << endl;
  cout << "Performance arithmétique: " << gflops << " GFLOP/s." << endl;
  

  cout << "mesh " <<  " Nsteps " <<  " DOF " <<  " Ttot " <<  " Tmoy " <<  " GFLP " << endl;
  cout << N << " " << Nsteps << " " << K * Np * Nfields << " " << time_total << " " << average << " " << gflops << endl;

  // ======================================================================
  // 3) Post-processing
  // ======================================================================*

  // Export final solution
  if (iOutputGmsh > 0)
    exportSolGmsh(N, _r, _s, _t, _EMsh, Nsteps, Nsteps * dt, _valQ);

  return 0;
}