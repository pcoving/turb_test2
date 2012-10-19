#include "CliffNoEdSolver.hpp"
#include "PenaltyTurbulence.hpp"

// =========================================================
// Cliff
// 
// Solves the incompressible constant density Navier-Stokes
// equations using a unstructured node-based fractional 
// step method.
// =========================================================


// the following class is an example of how you can customize cliff.

class MyCliff : public CliffNoEdSolver {
  
private:
  
  AverageOperator * average;

  double (*u_avg_zavg)[3];
  double (*u_rms_zavg)[3];
  double (*u_rey_zavg)[3];

  int nno_turbplane;
  double (*u_turbplane)[3];
  double (*x_no_plane)[3];
  int * ino_turbplane;
  
  MySyntheticTurbulence * myInflowTurb;

  double um_it[3];
  double Rd_it[3];
  double Rod_it[3];
  
 public:
  
  MyCliff() {
    
    if (mpi_rank == 0)
      cout << "MyCliff()" << endl;
        
    u_avg_zavg = NULL; registerData(u_avg_zavg,"U_AVG_ZAVG",NO_DATA);
    u_rms_zavg = NULL; registerData(u_rms_zavg,"U_RMS_ZAVG",NO_DATA);
    u_rey_zavg = NULL; registerData(u_rey_zavg,"U_REY_ZAVG",NO_DATA);

    if ( mpi_rank == 0 )
      cout << "initializing inflow turbulence with synthetic 3D turbulence" << endl;
    myInflowTurb = new MySyntheticTurbulence;
    myInflowTurb->init(u_turbplane);
  }  
  
  ~MyCliff() {
    delete[] x_no_plane;
  }

  void initialHook() {
    
    // set and initial condition in u and p and any scalars you have requested.
    // make sure the initial condition gets called only when step == 0...
    
    if (step == 0) {
      
      if (mpi_rank == 0)
	cout << "MyCliff::initialHook()" << endl;
      
      FOR_INO {
	u[ino][0] = 1.13;
        u[ino][1] = 0.0;
        u[ino][2] = 0.0;
        p[ino] = 0.0;
      }
    }
    
    int ed_mod = 0;    
    // add dissipation around source term
    FOR_IED {
      const int ino0 = nooed[ied][0];
      const int ino1 = nooed[ied][1];
      assert(ino0 >= 0 && ino1 >= 0);
      double x_ed[3] = {0.5*(x_no[ino0][0] + x_no[ino1][0]),
			0.5*(x_no[ino0][1] + x_no[ino1][1]),
			0.5*(x_no[ino0][2] + x_no[ino1][2])};
      if (x_ed[0] < -0.004) {
	ed_alpha[ied] = 1.0;
	++ed_mod;
      }
    }
    int ed_mod_global;
    MPI_Reduce(&ed_mod, &ed_mod_global, 1, MPI_INT, MPI_SUM, 0, mpi_comm);
    if (mpi_rank == 0) cout << " number of ed_alpha's modified: " << ed_mod_global << endl;

    double * v1 = new double[nno];
    FOR_INO {
      v1[ino] = x_no[ino][0];
    }
    average = new AverageOperator(v1, no_local_volume, nno, 1e-8);
    delete[] v1;
    
    // set the points
    nno_turbplane = 0;
    for (int iter = 0; iter < 2; ++iter) {
      /*
      FOR_FAZONE_BOUNDARY {
	if (fazone->getName() == "x0") {
	  for (int ino_b = fazone->ino_boundary_f; ino_b <= fazone->ino_boundary_l; ++ino_b) {
	    int ino = noono_boundary[ino_b];
	    if (iter == 1) {
	      FOR_I3 x_no_plane[nno_turbplane][i] = x_no[ino][i];
	      ino_turbplane[nno_turbplane] = ino;
	    }
	    ++nno_turbplane;
	  }
	}
      }
      */
      
      FOR_INO {
        if (fabs(x_no[ino][0]) < 1e-6) {
          if (iter == 1) {
            FOR_I3 x_no_plane[nno_turbplane][i] = x_no[ino][i];
            ino_turbplane[nno_turbplane] = ino;
          }
          ++nno_turbplane;
        }
      }
      
      if (iter == 0) {
        u_turbplane = new double[nno_turbplane][3];
        x_no_plane = new double[nno_turbplane][3];
        ino_turbplane = new int[nno_turbplane];
        nno_turbplane = 0;
      }
    }
    myInflowTurb->setPoints(x_no_plane, nno_turbplane);
    
    
    const double Um = 1.13;

    const double xod = 52.0;
    double Tu[3];
    Tu[0] = 0.80*pow(xod,-5.0/7.0);
    Tu[1] = 0.89*Tu[0];
    Tu[2] = 0.89*Tu[0];
    Tu[0] += 0.055;  // fudge factor!
    Tu[1] += 0.025;  // fudge factor!
    Tu[2] += 0.025;  // fudge factor!
    
    um_it[0] = Um;
    um_it[1] = 0.0;
    um_it[2] = 0.0;
    Rd_it[0] = Um*Um*Tu[0]*Tu[0];
    Rd_it[1] = Um*Um*Tu[1]*Tu[1];
    Rd_it[2] = Um*Um*Tu[2]*Tu[2];
    Rod_it[0] = 0.0;
    Rod_it[1] = 0.0;
    Rod_it[2] = 0.0;
    
    if (mpi_rank == 0) {
      cout << " diagonal of reynolds stress tensor: " << Rd_it[0] << " " << Rd_it[1] << " " << Rd_it[2] << endl;
      cout << " turbulence intensities: " << Tu[0] << " " << Tu[1] << " " << Tu[2] << endl;
    }

  }

  void temporalHook() {
    
    if (mpi_rank == 0 && step%check_interval == 0)
      cout << "MyCliff::temporalHook()" << endl;
    
    const double (*u_avg)[3] = getR2("U_AVG");
    const double (*u_rms)[3] = getR2("U_RMS");
    const double (*u_rey)[3] = getR2("U_REY");
    
    FOR_INO {
      FOR_I3 {
	u_avg_zavg[ino][i] = u_avg[ino][i];
	u_rms_zavg[ino][i] = u_rms[ino][i];
	u_rey_zavg[ino][i] = u_rey[ino][i];
      }
    }
    
    average->apply(u_avg_zavg);
    average->apply(u_rms_zavg);
    average->apply(u_rey_zavg);    
  }

  int momentumBcHook(double (*u_bc)[3], FaZone * fazone) {
    if (mpi_rank == 0 && step%check_interval == 0)
      cout << "MyCliff::momentumBcHook()" << endl;
    
    assert(myInflowTurb != NULL);
    myInflowTurb->update(dt);
    
    //myInflowTurb->getData(u_turbplane,nno_turbplane);
    myInflowTurb->getDataOnPoints(u_turbplane, x_no_plane, nno_turbplane);
    addReynoldsStressAndSync();
    
    assert(fazone->getName() == "x0");
    
    int ipoint = 0;
    for (int ino_b = fazone->ino_boundary_f; ino_b <= fazone->ino_boundary_l; ++ino_b) {
      // to get the associated local node...
      int ino = noono_boundary[ino_b];
      assert(ino == ino_turbplane[ipoint]);
      u_bc[ino_b][0] = u_turbplane[ipoint][0];
      u_bc[ino_b][1] = u_turbplane[ipoint][1];
      u_bc[ino_b][2] = u_turbplane[ipoint][2];      
      ++ipoint;
    }
    
    return(INJECT_BC);
  }

  void addReynoldsStressAndSync() {
    
    for (int ipoint = 0; ipoint < nno_turbplane; ++ipoint) {
      // build the Cholesky decomposition of the Reynolds stresses
      double a11 =  sqrt(Rd_it[0]);
      double a21 = 0.0;
      if ( a11 > 0.0 ) a21 = Rod_it[0] / a11;
      double a22 = sqrt( max(Rd_it[1]-a21*a21,0.0) );
      double a31 = 0.0;
      if ( a11> 0.0 ) a31 = Rod_it[1] / a11;
      double a32 = 0.0;
      if ( a22 > 0.0 ) a32 =  (Rod_it[2] - a21*a31) / a22;
      double a33 = sqrt(  max(Rd_it[2]-a31*a31-a32*a32,0.0) );

      // velocity component
      double my_u[3];
      my_u[0]  = a11 * u_turbplane[ipoint][0];
      my_u[1]  = a21 * u_turbplane[ipoint][0] + a22 * u_turbplane[ipoint][1];
      my_u[2]  = a31 * u_turbplane[ipoint][0] + a32 * u_turbplane[ipoint][1] + a33 * u_turbplane[ipoint][2];

      // copy
      FOR_I3 u_turbplane[ipoint][i] = my_u[i];
    }
    
    // sync across processors
    {
      double (*u_no_sync)[3] = new double[nno][3];
      FOR_INO FOR_I3 u_no_sync[ino][i] = 0.0;
      for (int ipoint = 0; ipoint < nno_turbplane; ++ipoint) {
        const int ino = ino_turbplane[ipoint];
        FOR_I3 u_no_sync[ino][i] = u_turbplane[ipoint][i];
      }
      updateNoDataActive(u_no_sync,REPLACE_GHOST_DATA_);
      for (int ipoint = 0; ipoint < nno_turbplane; ++ipoint) {
        const int ino = ino_turbplane[ipoint];
        FOR_I3 u_turbplane[ipoint][i] = u_no_sync[ino][i];
      }
      delete[] u_no_sync;
    }

    // add the mean
    for (int ipoint = 0; ipoint < nno_turbplane; ++ipoint)
      FOR_I3 u_turbplane[ipoint][i] += um_it[i];
  }
  
  void addMomentumSourceHook(double * Au, double (*Au_grad)[3], double (*rhs)[3]) {
    
    if (mpi_rank == 0 && step%check_interval == 0)
      cout << "MyCliff::addMomentumSourceHook()" << endl;
    
    assert(myInflowTurb != NULL);
    myInflowTurb->update(dt);
    
    //myInflowTurb->getData(u_turbplane,nno_turbplane);
    myInflowTurb->getDataOnPoints(u_turbplane, x_no_plane, nno_turbplane);
    addReynoldsStressAndSync();
    
    const double sigma = 1.0/dt;
    for (int ipoint = 0; ipoint < nno_turbplane; ++ipoint) {
      const int ino = ino_turbplane[ipoint];
      const int non_diag = noono_i[ino];
      Au[non_diag] += sigma*no_local_volume[ino];
      
      FOR_I3 {
        rhs[ino][i] += sigma*no_local_volume[ino]*u_turbplane[ipoint][i];
      }
    }

    
  }

};

// ======================================================
// main
// ======================================================

int main(int argc, char * argv[]) {
  
  try {
  
    // initialize the environment: basically 
    // initializes MPI and parameters...
    
    CTI_Init(argc,argv,"cliff.in");
    
    // the solver...
    
    {
      
      // Customized solver...
      
      MyCliff solver;
      
      solver.init();
      
      solver.run();

    }

    // finalize the environment: reports parameter usage, shuts down MPI...
    
    CTI_Finalize();
    
    // rudimentary error management...

  }
  catch (int e) {
    if (e == 0)
      CTI_Finalize();
    else
      CTI_Abort();
  }
  catch(...) {
    //cerr << "unhandled Exception.\n" << endl;
    CTI_Abort();
  }
    
  return(0);

}


