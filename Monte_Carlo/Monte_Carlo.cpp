#include "../Shared_Objects/Header.h"
#include "Monte_Carlo.h"

using namespace std;

Monte_Carlo::Monte_Carlo(double ell_tot, 
                        double BindingEnergy,
                        double k_diff, 
                        int seed, 
                        int Nlinker_max,
                        int dimension) : distrib(1, 10000000)
{
    IF(true) { cout << "MC : creator" << endl; }
    // srand(seed);
    generator.seed(seed);
    // ---------------------------------------------------------------------------
    // -----------------------constant of the simulation--------------------------
    ell = ell_tot;
    rho =0.;
    slide = false;
    LoopLinkWrap::dimension=dimension;
    Nlinker= Nlinker_max;
    Linker::counter = 0;
    binding_energy = BindingEnergy;
    loop_link.create_new_occupied_linker(0.,0.,0.);
    Linker* R0 = loop_link.get_linkers().at({0.,0.,0.});
    Dangling dummy_dangling(R0, 0., ell, rho,slide); // dummy dangling that helps generate crosslinkers but has none initially
    Strand* dummy_strand(loop_link.Create_Strand(dummy_dangling));
    // ---------------------------------------------------------------------------
    //-----------------------------initialize crosslinkers------------------------
    set<array<double,3>> linkers(generate_crosslinkers(0));
  
    for(auto& linker : linkers)
    {
      loop_link.create_new_free_linker(linker.at(0),linker.at(1),linker.at(2));// the linker is free by default
    }
    loop_link.Remove_Strand(dummy_strand);
    // ---------------------------------------------------------------------------
    //-----------------------------initialize dangling----------------------------
    IF(true){cout<< "MC : create dangling" << endl;}
    loop_link.Create_Strand(Dangling(R0, 0., ell, rho,slide));
    IF(true) { cout << "MC : created" << endl; }
}
Monte_Carlo::~Monte_Carlo()
{


}

void Monte_Carlo::Evolve(bool* accepted, int* move_type, double* DE)
{
    // propose a move
        
    // compute the energy difference

    // apply the move if DE < 0

    // refuse the move if  DE>0
}

set<array<double,3>> Monte_Carlo::generate_crosslinkers(bool remake){
// loop over all the strand, and generate crosslinkers within their ellipse.
  set<array<double,3>> res;
  int N_linker_to_make(0.);
  // if we generate a fixed number of linkers
  if(Nlinker>0)
    {
      Linker::counter = loop_link.get_linker_size()-loop_link.get_N_free_linker();// number of occupied linkers      
      N_linker_to_make = Nlinker-Linker::counter; // total number of linkers that has to be added
      //cout<<"counter : "<<Linker::counter<<endl;
      //cout<<"N_linker_max :"<<N_linker_max<<endl;      
      // compute the probability to place N linkers propto the volume of each strands
      double total_volume(0.);
      for(auto& strand : loop_link.get_strands()){total_volume+=strand->get_V();}
      for(auto& strand : loop_link.get_strands())
      {
        double a,b;
        array<double,3> ctr_mass,main_ax;
        // number of linker to add to this specific strand
        int Nlinkers_strand(round(strand->get_V()/total_volume * N_linker_to_make));
        //cout<<strand->get_V()/total_volume * N_linker_to_make<<endl;
        //cout<<Nlinkers_strand<<endl;
        strand->get_volume_limit(main_ax,ctr_mass,a,b);
        generate_point_in_ellipse(main_ax,ctr_mass,a,b,res,Nlinkers_strand);
      }
    }
  else
  {
  for(auto & strand : loop_link.get_strands())
  {
    double a,b;
    array<double,3> ctr_mass,main_ax;
    poisson_distribution<int> distribution(rho * strand->get_V());
    cout<<"volume of the strand "<<strand->get_V()<<endl;
    
    if(remake){N_linker_to_make =N_linker_to_make = max(0.,distribution(generator)-(double)strand->get_occ_r().size());}
    else{N_linker_to_make = max(0.,distribution(generator)-(double)strand->get_r().size()-(double)strand->get_occ_r().size());}
    
    strand->get_volume_limit(main_ax,ctr_mass,a,b);
    generate_point_in_ellipse(main_ax,ctr_mass,a,b,res,N_linker_to_make);
  }
  }
  // transform res depending on the dimension
  set<array<double,3>> dimensional_res;
  if (LoopLinkWrap::dimension == 2){for(auto& xyz: res){dimensional_res.insert({xyz[0],xyz[1],0.});}}
  else if(LoopLinkWrap::dimension==1){for(auto& xyz: res){dimensional_res.insert({xyz[0],0.,0.});}}
  else{return res;}
  return dimensional_res;
}