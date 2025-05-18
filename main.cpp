// Compilation : make
#include "Connectivity.hpp"
#include "Gmsh.hpp"
#include "Nodes.hpp"
#include "Numerics.hpp"
#include "Tools.hpp"

#include <iostream>
#include <iomanip>  // needed for std::setprecision
#include <omp.h> // OpenMP

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
  // double CFL = 1;                 // CFL
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
  omp_set_num_threads(4);

  double duration_tot = 0.0;
  double duration_moy = 0.0;
  cout << " Nsteps DOF Ttot Tmoy GFLOPs" << endl;
  // Global time iteration
  for (int nGlo = 0; nGlo < Nsteps; ++nGlo) {
    int NbOp = 0; // compteur de FLOP par itération
    double timeBegin = dsecnd(); // time 
    double runTime = nGlo * dt; // Time at the beginning of the step

    // Local time iteration
    for (int nLoc = 0; nLoc < 5; ++nLoc) {
      double a = rk4a[nLoc];
      double b = rk4b[nLoc];

      #pragma omp parallel for schedule(static)

      // ======================== (1) UPDATE RHS
      for (int k = 0; k < K; ++k) { // boucle sur les élements du maillage

        double c = _c[k];
        double rho = _rho[k];

        // ======================== (1.1) UPDATE PENALTY VECTOR 


        int size_flux_vect = Nfp * NFacesTet; // we compute it once!
        vector<double> s_p_flux(size_flux_vect);
        vector<double> s_u_flux(size_flux_vect);
        vector<double> s_v_flux(size_flux_vect);
        vector<double> s_w_flux(size_flux_vect);

        for (int f = 0; f < NFacesTet; f++) { // boucle sur faces 

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
          for (int nf = f * Nfp; nf < (f + 1) * Nfp; nf++) { // boucle sur les points de f

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
              // NbOp += (8+6*3); //* (f + 1) * NFacesTet * K * nLoc * nGlo;
            } // end if 
          } // end for nf (nf < (f + 1))
        } // end for f (f < NFacesTet)

        // ======================== (1.2) COMPUTING VOLUME TERMS
        // calcul des dérivées spatiales 

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
          
          int vol_ter = Nfields * (k * Np + n); // we just compute it once!
          s_p[n] = _valQ[vol_ter + 0];
          s_u[n] = _valQ[vol_ter + 1];
          s_v[n] = _valQ[vol_ter + 2];
          s_w[n] = _valQ[vol_ter + 3];
        } // end for n (n < Np)

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
            
            // NbOp += 8*3;
          } // end for m (m < Np)

          double dpdx = rx * dpdr + sx * dpds + tx * dpdt;
          double dpdy = ry * dpdr + sy * dpds + ty * dpdt;
          double dpdz = rz * dpdr + sz * dpds + tz * dpdt;
          double dudx = rx * dudr + sx * duds + tx * dudt;
          double dvdy = ry * dvdr + sy * dvds + ty * dvdt;
          double dwdz = rz * dwdr + sz * dwds + tz * dwdt;
          double divU = dudx + dvdy + dwdz;

          // NbOp += 5*5 + 2;

          // Compute RHS (only part corresponding to volume terms)
          int vol_ter = Nfields * (k * Np + n); // we just compute it once!
          _rhsQ[vol_ter + 0] = -c * c * rho * divU;
          _rhsQ[vol_ter + 1] = -1. / rho * dpdx;
          _rhsQ[vol_ter + 2] = -1. / rho * dpdy;
          _rhsQ[vol_ter + 3] = -1. / rho * dpdz;
          // NbOp += 3 + 2*3;
        } // end for n (n < Np)

        // ======================== (1.3) COMPUTING SURFACE TERMS
        // opérateur lift

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
              // NbOp +=8;
            } // end for m  (m < (f + 1))

            // Load geometric factor
            double Fscale = _Fscale[k * NFacesTet + f];

            // Update RHS (with part corresponding to surface terms)
            int elem_k = k * Np * Nfields + n * Nfields;
            _rhsQ[elem_k + 0] -= p_lift * Fscale;
            _rhsQ[elem_k + 1] -= u_lift * Fscale;
            _rhsQ[elem_k + 2] -= v_lift * Fscale;
            _rhsQ[elem_k + 3] -= w_lift * Fscale;
            // NbOp += 8;
          } // end for f (f < NFacesTet)
        } // end for n (n < Np)
      } // end for k (k < K)

      // fin RHS update
      // ======================== (2) UPDATE RESIDUAL + FIELDS

      #pragma omp parallel for schedule(static)

      for (int k = 0; k < K; ++k) {
        for (int n = 0; n < Np; ++n) {
          for (int iField = 0; iField < Nfields; ++iField) {
            int id = k * Np * Nfields + n * Nfields + iField; // 5 opérations il y en a qu'on peut éviter!

            _resQ[id] = a * _resQ[id] + dt * _rhsQ[id];
            _valQ[id] = _valQ[id] + b * _resQ[id];
            // NbOp += 3 + 2;
          } // end for iField (iField < Nfields)
        } // end for n (n < Np)
      } // end for k (k < K)
    } // end for nLoc (nLoc < 5)

    double timeEnd = dsecnd(); // time 
    double duration = timeEnd - timeBegin;
    duration_tot += duration;
    // cout << "dsecnd" << timeEnd - timeBegin << endl;

    // Export solution
    if ((iOutputGmsh > 0) && ((nGlo + 1) % iOutputGmsh == 0))
      // exportSolGmsh(N, _r, _s, _t, _EMsh, nGlo + 1, runTime + dt, _valQ);
      exportSolGmsh(N, _r, _s, _t, _EMsh, nGlo + 1, duration, _valQ);
    
    // if (nGlo % 100 == 0) {
    //   cout << "iter : " << nGlo << endl;
    //   // cout << "FLOP : " << NbOp << endl;
    // }
    } // end for nGlo (nGlo < Nsteps)

  
  // Temps d'exécution moyen par itérations, temps total d'exécution
  
  double average = duration_tot/Nsteps;
  
  double flops_par_iter = 1600.0 + 68.0*Np;
  double total_flops = 5.0 * Nsteps * K * flops_par_iter;
  double gflops = total_flops/(duration_tot*1e9);
  //                         Nsteps           DOF                        Ttot                   Tmoy           GFLOP
  cout << fixed << setprecision(6) << Nsteps << " " << K * Np * Nfields << " " << duration_tot << " " << average << " " << gflops << endl;
  // cout << Nsteps * K * (N^6 + N^5 * NFacesTet + N^3 * Nfields)/(duration_tot*1e9) << endl;

  // ======================================================================
  // 3) Post-processing
  // ======================================================================
  // cout << "Duration totale : " << duration_tot << endl;

  // Export final solution
  if (iOutputGmsh > 0){
    exportSolGmsh(N, _r, _s, _t, _EMsh, Nsteps, Nsteps * dt, _valQ);
  }

  return 0;
}
