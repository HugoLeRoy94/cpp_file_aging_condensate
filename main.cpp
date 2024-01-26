#include "Shared_Objects/Header.h"
#include "Gillespie/Gillespie.h"
using namespace std;
int main(int argc, char* argv[]){
  int t_tot(1000);
  double ell_tot(100.);
  //double distance_anchor(1000.);
  double rho0(pow(10,-2));
  double BindingEnergy(-10);
  int bind(0.);
  double* R;
  Gillespie* S = new Gillespie(ell_tot,0.,BindingEnergy,5.20142*pow(10,-5),5984365,false,5,3);
  cout<<"beginning of the simulation"<<endl;
  for(int n(0);n<t_tot;n++){
  //cout<<n<<endl;
  double time(S->evolve(&bind));
  //S->reset_crosslinkers();
  //if(S->get_r_gillespie_size() != 9)
  //{ 
  //  cout<<S->get_r_gillespie_size()<<endl;
  //  cout<<n<<endl;
  //  exit(0);
  //}
 S->Print_Loop_positions();
  }
  delete S;
  //set<array<double,3>> res;
  //array<double,3> main_ax = {1.,1.,1.};
  //array<double,3> ctr_mass = {1.,2.,3.};
  //double a = 25;
  //double b = 7.5;
//
  //generate_point_in_ellipse(main_ax,ctr_mass,a,b,res,1000);
//
  //ofstream file;
  //file.open("data.txt",ios::trunc);
  //for(auto pts  : res)
  //{
  //  file<<pts[0]<<" "<<pts[1]<<" "<<pts[2]<<endl;;
  //}
  }
