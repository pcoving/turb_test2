// ===========================
// inflow turbulence generator
// ===========================

#include "CTI.hpp"
using namespace CTI;

// fft stuff
# include "FFT.hpp"

# include "ArrayNd.hpp"
# include <complex>
# include <vector>
# include <time.h> 
# include <stdlib.h>

// inflow plane type
// XXXX (already delcared in inflowturb...)
//enum InflowGeom {
//  PLANE_YZ_GEOM,
//  PLANE_XY_GEOM
//};


// ======================================
// Inflow data on a 2D plane for a single 
// component such as u, v, or w
// ======================================
class MyInflowDataDF {

private:

  array< complex<double> > * b_hat;       // filter coeffs in Fourier space
  array< complex<double> > * r_hat;       // random number in Fourier space
  int Nf1, Nf2;                 // filter half width
  int n1, n2;                   // length scale 
  int Nr1, Nr2;                 // array size
  fftw_plan fftw_plan_calcData; // fftw plan for inflow computation
  array<double> * ur;           // extended field
  double tau;                   // timescale in the norml direction  
  array<double> * v_old;        // inflow data from previous step
  int block_mpi_rank;           // mpi rank on which the data is computed
  int N1, N2;                   // box dimension N1 x N2

public:

  array<double> * v;
      
  // ===========
  // constructor
  // =========== 
  MyInflowDataDF(const int block_mpi_rank){
    
    if (mpi_rank == 0)
      cout << "MyInflowDataDF():" << endl;

    // set data
    this->block_mpi_rank = block_mpi_rank;

    // nullify data  
    b_hat = NULL;
    r_hat = NULL;
    ur    = NULL;
    v     = NULL;
    v_old = NULL;

  }

  // =============
  // deconstructor
  // =============
  ~MyInflowDataDF(){

    if (mpi_rank == 0)
      cout << "~MyInflowDataDF():" << endl;
    
    if (mpi_rank == block_mpi_rank)
      fftw_destroy_plan(fftw_plan_calcData);

    if ( b_hat != NULL ) delete b_hat;
    if ( r_hat != NULL ) delete r_hat;
    if ( ur    != NULL ) delete ur;
    if ( v     != NULL ) delete v;
    if ( v_old != NULL ) delete v_old;
    
  }


  // ====================================
  // initialization done on one processor
  // ====================================
  void init(const int N1,      const int N2, 
	    const double& dx1, const double& dx2, 
	    const double& L1,  const double& L2, 
	    const double& tau) {

    // copy
    this->N1  = N1;
    this->N2  = N2;
    this->tau = tau;
    
    // check
    assert(N1>2);
    assert(N2>2);
    
    n1  = (int)(L1/dx1);
    n2  = (int)(L2/dx2);
    
    // set the filter width to be 2*n
    Nf1 = 2*n1;
    Nf2 = 2*n2;
    // compute array size
    Nr1 = N1 + 2*Nf1;
    Nr2 = N2 + 2*Nf2;
    
    // allocate memory
    v     = new array<double>(N1,N2);
    v_old = new array<double>(N1,N2);

     // only the dedicated rank allocates memory and computes the coeffs
    if (mpi_rank == block_mpi_rank) {
        
      // allocate memory
      b_hat = new array< complex<double> >(Nr1,Nr2/2+1);
      r_hat = new array< complex<double> >(Nr1,Nr2/2+1);
      ur    = new array<double>(Nr1,Nr2);
      
      // build the fft plan
      fftw_plan_calcData =						\
	fftw_plan_dft_c2r_2d(Nr1,Nr2,					\
			     reinterpret_cast<fftw_complex*>(&((*r_hat)(0,0))), \
			     &((*ur)(0,0)),FFTW_ESTIMATE);
      
      // build filter coeffs
      initFilter();
    }

    // dump information
    if (mpi_rank == 0 ) {
      cout << "    > computing mpi_rank              = "  << block_mpi_rank << endl;
      cout << "    > turbulence length scales        = (" << L1 << "," << L2 << ")" << endl;
      cout << "    > points per turb length scale    = (" << n1 << "," << n2 << ")" << endl;
      cout << "    > points per half filter width    = (" << Nf1 << "," << Nf2 << ")" << endl;
      cout << "    > time scale in normal direction  = "  << tau << endl;
    }
    
  }


  // ====================
  // calculate inflow data
  // ====================
  void regenerate(const double dt){

    //if ( mpi_rank == 0 )
    //  cout << "regenerate:" << endl;

    double coeff_old = 0.0;
    double coeff_new = 1.0;
    if ( tau != 0.0 ){
      coeff_old = exp(-0.5*M_PI*dt/tau);
      coeff_new = sqrt(1.0-exp(-M_PI*dt/tau));
    }
       
    // copy the array to old one
    for ( int i=0; i<N1; ++i)
      for ( int j=0; j<N2; ++j) 
	(*v_old)(i,j) = (*v)(i,j);
    
    // update the inflow turbulence
    if (mpi_rank == block_mpi_rank){
      // first do it in the plane
      calcInflowDataF(v);
      // now appy the time correlation
      for ( int i=0; i<N1; ++i)
	for ( int j=0; j<N2; ++j) {
	  (*v)(i,j) *= coeff_new;
	  (*v)(i,j) += coeff_old*(*v_old)(i,j); 
	}   
    }
    
  }


  // ==================
  // update inflow data
  // ==================
  void update(){
    // send the data to all processors
    MPI_Bcast(&((*v)(0,0)),v->getSize(),MPI_DOUBLE,block_mpi_rank,mpi_comm);
  }

private:  
 
  // ===================
  // build filter coeffs
  // ===================
  void initFilter(){
    
    // compute the fourier transform of the filter coeffs
    // build the filter coeffs
    // b1
    int Nb1 = 2*Nf1+1;
    double * b1 = new double[Nb1];
    double tmp = 0.0;
    for ( int i=0; i<Nb1; ++i){
      int k = i-Nf1;
      b1[i] = exp(-M_PI*fabs((double)k)/(double)n1);
      //b2[i] = exp(-0.5*M_PI*(double)k * (double)k / ( (double)n1 * (double)n1 ));
      tmp += b1[i]*b1[i];
    }
    for ( int i=0; i<Nb1; ++i)
      b1[i] /= sqrt(tmp);  

    // b2
    int Nb2 = 2*Nf2+1;
    double * b2 = new double[Nb2];
    tmp = 0.0;
    for ( int i=0; i<Nb2; ++i){
      int k = i-Nf2;
      b2[i] = exp(-M_PI*fabs((double)k)/(double)n2);
      //b2[i] = exp(-0.5*M_PI*(double)k * (double)k / ( (double)n2 * (double)n2 ));
      tmp += b2[i]*b2[i];
    }
    for ( int i=0; i<Nb2; ++i)
      b2[i] /= sqrt(tmp);  

    // pad the filter with zeros
    
    // b1
    double * b1r = new double[Nr1];
    for (int i=0; i<Nr1; ++i) b1r[i] = 0.0;
    for ( int i=-Nf1; i<=Nf1; ++i){
      int k = i;
      if ( k>=Nr1 ) k -=Nr1;
      if ( k<0    ) k +=Nr1;
      b1r[k] = b1[i+Nf1];
    }
    delete[] b1;

    // b2
    double * b2r = new double[Nr2];
    for (int i=0; i<Nr2; ++i) b2r[i] = 0.0;
    for ( int i=-Nf2; i<=Nf2; ++i){
      int k = i;
      if ( k>=Nr2 ) k -=Nr2;
      if ( k<0    ) k +=Nr2;
      b2r[k] = b2[i+Nf2];
    }
    delete[] b2;

    // 2D filter is simply b1*b2
    for (int i=0; i<Nr1; ++i)
      for (int j=0; j<Nr2; ++j)
	(*ur)(i,j) = b1r[i]*b2r[j];
    delete[] b1r;
    delete[] b2r;

    // now take the Fourier Transform
    {
      fftw_plan plan =							\
	fftw_plan_dft_r2c_2d(Nr1,Nr2,&((*ur)(0,0)),			\
			     reinterpret_cast<fftw_complex*>(&((*b_hat)(0,0))),FFTW_ESTIMATE);
      fftw_execute(plan);
      fftw_destroy_plan(plan);
    }

  }


  // ===============================================
  // compute the inflow data using Fourier transform
  // ===============================================
  void calcInflowDataF(array<double> * u){

    // compute signal 

    int this_Nr2 = Nr2/2+1;
    complex<double> my_i(0.0,1.0);
    for (int i=0; i<Nr1; ++i) {
      for (int j=0; j<this_Nr2; ++j) {
	double tmp = M_PI*(2.0*(double)(rand())/(double)(RAND_MAX)-1.0);
	(*r_hat)(i,j) = exp(my_i*tmp);
      }
    }
    // remove the mean
    complex<double> my_z(0.0,0.0);
    (*r_hat)(0,0) = my_z;
    
    //convolution in time
    for (int i=0; i<Nr1; ++i) 
      for (int j=0; j<this_Nr2; ++j) 
	(*r_hat)(i,j) *= conj((*b_hat)(i,j));
    
    
    // now take the Fourier Transform
    fftw_execute(fftw_plan_calcData);
    

    // here, we need to set the mean and std to 
    // zero and one, respectively. Note that since
    // after filtering there is no gaurantee that
    // ur has a unity variance
    double mymean = 0.0;
    double myvar  = 0.0;
    for (int i=0; i<Nr1; ++i)
      for (int j=0; j<Nr2; ++j){
	mymean += (*ur)(i,j);
	myvar  += (*ur)(i,j)*(*ur)(i,j);
      }
    double coeff = 1.0 / (double)(Nr1*Nr2);
    mymean *= coeff;
    myvar  *= coeff;
    myvar  -= mymean*mymean;
    double normFac = 1.0/sqrt(myvar);
    for (int i=0; i<Nr1; ++i)
      for (int j=0; j<Nr2; ++j){
	(*ur)(i,j) -= mymean;
	(*ur)(i,j) *= normFac; 
      }  
    
    // copy
    for (int i=Nf1; i<Nr1-Nf1; ++i){
      int k = i - Nf1;	
      for (int j=Nf2; j<Nr2-Nf2; ++j){
	int l = j - Nf2;
	(*u)(k,l) = (*ur)(i,j);
      }
    }

  }

};

// ===================================
// Inflow data for for a 2D plane with 
// three velocity components (u,v,w)
// ===================================
class MyInflowDataDF3D {
  
private:

  MyInflowDataDF * u1;
  MyInflowDataDF * u2;
  MyInflowDataDF * un;
  
public:
  
  string     name;
  InflowGeom Geom;
  int        N1, N2;
  double     dx1, dx2;
  double     x1_min, x1_max;
  double     x2_min, x2_max;

  // array access
  array<double> * v[3];


  MyInflowDataDF3D() {
    if ( mpi_rank == 0 ) 
      cout << "MyInflowDataDF3D():" << endl;

    u1 = NULL;
    u2 = NULL;
    un = NULL;

    FOR_I3 v[i] = NULL;
  }


  ~MyInflowDataDF3D() {
    if ( mpi_rank == 0 ) 
      cout << "~MyInflowDataDF3D():" << endl;

    if ( u1 != NULL) delete u1;
    if ( u2 != NULL) delete u2;
    if ( un != NULL) delete un;

  }

  // ==========
  // initialize
  // ==========
  void init(int block_mpi_rank) {

    // check the rank
    assert(block_mpi_rank>=0);
    assert(block_mpi_rank<=mpi_size);

    // constructor
    if ( block_mpi_rank == mpi_size) block_mpi_rank = 0;  
    u1 = new MyInflowDataDF(block_mpi_rank); ++block_mpi_rank;

    if ( block_mpi_rank == mpi_size) block_mpi_rank = 0;  
    u2 = new MyInflowDataDF(block_mpi_rank); ++block_mpi_rank;

    if ( block_mpi_rank == mpi_size) block_mpi_rank = 0;  
    un = new MyInflowDataDF(block_mpi_rank); ++block_mpi_rank;
    
    // velocity and lenght scale normal to the plane
    double Un;
    double Ln_1, Ln_2, Ln_n;
    double L1_1, L2_1, L1_2, L2_2, L1_n, L2_n;
    
    // rather than looking in input file, just hard-code for now...
        
    if (mpi_rank == 0)
      cout << " > initializing inflow turbulence " << endl;
    
    // these should be in input file...
    {
      Geom = PLANE_YZ_GEOM;
      Un = 1.13;
      const double xod = 52.0;
      const double d = 0.0127;
      const double I = 0.2;
      const double Li_i = d*(I*sqrt(xod));
      const double Li_j = 0.5*Li_i;
      
      L1_1  = Li_i;
      L1_2  = Li_j;
      L1_n  = Li_j;
      L2_1  = Li_j;
      L2_2  = Li_i;
      L2_n  = Li_j;
      Ln_1  = Li_j;
      Ln_2  = Li_j;
      Ln_n  = Li_i;
    }

    // range
    x1_min  = -0.3;
    x1_max  = 0.3;
    x2_min  = -0.3;
    x2_max  = 0.3;
  
    // calculate dx1 and dx2 assuming we have minimum 20 points per integral length scale
    {
      double minL  = min(L1_1, min(L1_2, L1_n));
      //double dx    = minL / 20.0;
      double dx    = minL / 20.0;
      this->N1     = (int)((x1_max-x1_min)/dx) + 1;
      this->dx1    = (x1_max-x1_min) / (double)(N1-1);    
    }
    {    
      double minL  = min(L2_1, min(L2_2, L2_n));
      //double dx    = minL / 20.0;
      double dx    = minL / 20.0;
      this->N2     = (int)((x2_max-x2_min)/dx) + 1;
      this->dx2    = (x2_max-x2_min) / (double)(N2-1);    
    }      
    
    // check
    assert(N1>2);
    assert(N2>2);
    /*
    if ( (N1>2000) || (N2>2000) ) {
      cerr << "Error: block size (N1,N2) = (" << N1 << "," << N2 << ") is larger than the limit (2000,2000)" << endl;
      cerr << "   choose a smaller block dimensitons or a larger turbuelence lenght scale" << endl;
      throw(0);
    }
    */
    // initialze the blocks
    if ( mpi_rank == 0 ) 
      cout << "  > for direction 1:" << endl;
    double tau = Ln_1/Un;  
    u1->init(N1,N2,dx1,dx2,L1_1,L2_1,tau);
  
    // initialze the blocks
    if ( mpi_rank == 0 ) 
      cout << "  > for direction 2:" << endl;
    tau = Ln_2/Un;  
    u2->init(N1,N2,dx1,dx2,L1_2,L2_2,tau);
    
    // initialze the blocks
    if ( mpi_rank == 0 ) 
      cout << "  > for direction n:" << endl;
    tau = Ln_n/Un;  
    un->init(N1,N2,dx1,dx2,L1_n,L2_n,tau);
    
    // set the pointer'   
    if ( Geom == PLANE_YZ_GEOM ) {
      v[0] = un->v;     
      v[1] = u1->v;
      v[2] = u2->v;    
    }
    else if ( Geom == PLANE_XY_GEOM ) {
      v[0] = u1->v;
      v[1] = u2->v;
      v[2] = un->v;         
    }
    else {
      cerr << "Error: unrecognized geometry " << Geom << endl;
      throw(0);
    }

    // dump information
    if (mpi_rank == 0 ){
      cout << "  > x1 range                          = (" << x1_min << "," << x1_max << ")" << endl;
      cout << "  > x2 range                          = (" << x2_min << "," << x2_max << ")" << endl;
      cout << "  > grid spacing                      = (" << dx1 << "," << dx2 << ")" << endl;      
      cout << "  > gid size                          = (" << N1 << "," << N2 << ")" << endl;
      cout << "  > turbulence length scale(1)        = (" << L1_1 << "," << L1_2 << "," << L1_n <<")" << endl;
      cout << "  > turbulence length scale(2)        = (" << L2_1 << "," << L2_2 << "," << L2_n <<")" << endl;
      cout << "  > turbulence length scale(n)        = (" << Ln_1 << "," << Ln_2 << "," << Ln_n <<")" << endl;
      cout << "  > mean velocity in normal direction = "  << Un << endl;
      cout << " ********* end of information for inflow block: " << name << " *********\n" << endl;
    }
       
  }


  // ==========
  // regenerate 
  // ==========
  void regenerate(const double& dt ) {
    un->regenerate(dt);
    u1->regenerate(dt);
    u2->regenerate(dt);  
  }

  // =======
  // update 
  // =======
  void update() {
    un->update();
    u1->update();
    u2->update();
  }


};


// =================
// main inflow class
// =================

typedef struct {
  int nno;
  int (*ijk)[3];
} myijknno;

class MySyntheticTurbulence {

private:

  double (**u)[3];
  MyInflowDataDF3D * data;
  myijknno ijknno;

public:

  MySyntheticTurbulence(){
    if ( mpi_rank == 0 )
      cout << "MySyntheticTurbulence():" << endl;

    u = NULL;
  }
  
  ~MySyntheticTurbulence(){
    if ( mpi_rank == 0 )
      cout << "~MySyntheticTurbulence():" << endl;
  }

  // ==========
  // initialize 
  // ==========
  void init(double (*&u)[3]) {
    
    // copy the data pointer
    this->u = &u;
    
    // initialize
    data = new MyInflowDataDF3D;
    data->init(0);
  }
  
  // ==================
  // update inflow data
  // ==================
  void update(const double& dt){
    
    data->regenerate(dt);
    data->update();
        
  }

  void getData(double (*v)[3],const int n){
    
    assert(ijknno.nno == n);
    
    for ( int ino=0; ino < ijknno.nno; ++ino){
      int k = ijknno.ijk[ino][2];
      
      if ( k >= 0 ) {
        // indecies
        int i = ijknno.ijk[ino][0];
        int j = ijknno.ijk[ino][1];
	
        for ( int l=0; l<3; ++l)
          v[ino][l] = (*(data->v[l]))(i,j);
      }
      // otherwise set it to zero
      else
        FOR_I3 v[ino][i] = 0.0;
    }

  }

  void getDataOnPoints(double (*v)[3], const double (*x)[3], const int nno) {

    // hard coded for z streamwise direction...
    for (int ino = 0; ino < nno; ++ino) {
      assert((x[ino][1] <= data->x1_max) &&
	     (x[ino][1] >= data->x1_min) &&
	     (x[ino][2] <= data->x2_max) &&
	     (x[ino][2] >= data->x2_min));
      
      int ii = floor((x[ino][1] - data->x1_min)/data->dx1 + 0.5);
      assert((ii >= 0) && (ii < data->N1));
      int jj = floor( (x[ino][2] - data->x2_min)/data->dx2 + 0.5);
      assert((jj >= 0) && (jj < data->N2));
      
      FOR_I3 v[ino][i] = (*(data->v[i]))(ii, jj);
    }
    
  }

  void setPoints(const double (*x)[3], const int nno){

    if ( mpi_rank == 0 )
      cout << "setPoints: " << endl;
    
    // allocate memory
    ijknno.nno = nno;
    if ( nno > 0 )
      ijknno.ijk = new int[nno][3];

    // for now if a point belongs to two inflow
    // blocks, just pick the first block on the list
    // we might want to do averaging later ..
    for (int ino=0; ino < nno; ++ino){
      
      // initialize the array
      ijknno.ijk[ino][0] = -1;
      ijknno.ijk[ino][1] = -1;
      ijknno.ijk[ino][2] = -1;

      // set the coordinates
      double x1, x2;
      if ( data->Geom == PLANE_YZ_GEOM ) {
	x1 = x[ino][1];
	x2 = x[ino][2];
      }
      else if ( data->Geom == PLANE_XY_GEOM ) {
	x1 = x[ino][0];
	x2 = x[ino][1];
      }
      else {
	cerr << "Error: unrecognized geometry: " << data->Geom << endl;
	throw(0);
      }

      // if the point is inside the block
      if ( (data->x1_min <= x1)  &&
	   (data->x1_max >= x1)  &&
	   (data->x2_min <= x2)  &&
	   (data->x2_max >= x2)       ){
	// vector
	// (zero for now...)
	ijknno.ijk[ino][2] = 0;
	// i index
	int index = (int)( (x1-data->x1_min)/data->dx1 );
	assert((index>=0)&&(index<data->N1));
	// length
	double dx = (x1-data->x1_min)/data->dx1 - (double)index;
	if ( dx<0.5 )
	  ijknno.ijk[ino][0] = index;
	else
	  ijknno.ijk[ino][0] = index+1;
	assert(ijknno.ijk[ino][0]<data->N1);
	// j index
	index = (int)( (x2-data->x2_min)/data->dx2 );
	assert((index>=0)&&(index<data->N2));
	// length
	dx = (x2-data->x2_min)/data->dx2 - (double)index;
	if ( dx<0.5 )
	  ijknno.ijk[ino][1] = index;
	else
	  ijknno.ijk[ino][1] = index+1;
	assert(ijknno.ijk[ino][1]<data->N2);
      }

    }

  }


};

