#include "Header.h"
using namespace std;
Dangling::Dangling() : Strand()
{}
// dsl : distance sur ell : ratio of the two usefull to know how to adjust the concentration of linkers
string Dangling::whoIam() const{return "dangling";}
Dangling::Dangling(Linker* R0,
                  Linker* R1,
                  double ell_0,
                  double ell_in,
                  double rho,
                  bool sliding ) : Strand(R0,R1,ell_0,rho,sliding)
{
  IF(true) { cout << "Dangling : creator" << endl; }
  ell = ell_in;
  radius = sqrt(ell/2.);
  rho0 = rho;
  //V = 4 / 3 * Pi * pow(2 * ell, 1.5);
  V = 4./3.*Pi*pow(radius,3); // linkers are generated into a squared box
  if(R0!=nullptr){
  xg = R0->r().at(0);
  yg = R0->r().at(1);
  zg = R0->r().at(2);
  }
  else{
  xg = R1->r().at(0);
  yg = R1->r().at(1);
  zg = R1->r().at(2);
  }
  IF(true) { cout << "Dangling : constructor over." << endl; }
}
Dangling::Dangling(const Dangling &dangling,
                  Linker* new_left_linker,Linker* new_right_linker) : Strand(dangling,new_left_linker,new_right_linker)
{
  IF(true) { cout << "Dangling : copy constructor" << endl; }
  radius = dangling.radius;
}

Dangling::Dangling(const Dangling &dangling) : Strand(dangling)
{
  IF(true) { cout << "Dangling : copy constructor" << endl; }
  radius = dangling.radius;
}
Strand* Dangling::clone() const{return new Dangling(*this);}
// ---------------------------------------------------------------------------
//-----------------------------------accessor---------------------------------
// ---------------------------------------------------------------------------
double Dangling::get_S(double dl) const { return (ell+dl) * log(4 * Pi); }
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
double Dangling::Omega(double ell) const
{
  return pow(4 * Pi, ell);
}
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
//________________________geometric child_dependend computation_____________
array<double, 3> Dangling::random_in_volume()
{
  uniform_real_distribution<double> distribution(-1, 1); // doubles from -1 to 1
  double x(distribution(generator) * radius);
  double y(distribution(generator) * radius);
  double z(distribution(generator) * radius);
  array<double, 3> res{x + xg, y + yg, z + zg};
  return res;
}

void Dangling::get_volume_limit(array<double,3>& main_ax, 
                                array<double,3>& ctr_mass,
                                double& a, double& b) const
{
  ctr_mass={xg,yg,zg};
  a=radius;
  b=radius;
  main_ax = {1.,1.,1.}; // the main ax does not matter when a = b
}
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
double Dangling::compute_binding_rate(double li, Linker* rlinker) const {
    // Check if both Rleft and Rright are null, if yes, return 0.0
    if (Rleft == nullptr && Rright == nullptr) {
        return 0.0;
    }
    // Determine which R to use
    Linker* currentR = (Rleft != nullptr) ? Rleft : Rright;
    // Pre-calculate some values
    double r_diff = get_square_diff(currentR->r(), rlinker->r());
    double r_li_ratio = r_diff / li;
    // Check the condition
    if (r_diff > li) {
        return 0.0;
    }
    // Compute the result
    // OK, I don't understand why there isn't 1/ell that multiply the rate. Thus, I will remove it, but maybe it was correct.
    // The explanation I gave back then is left intact as well. Understand who can.
    // Maybe I understood: Assuming that each monomer can diffuse once per time step, then a complete change of micro-state takes ell steps.
    //double rate = exp(1.5 * log(1.5 / (Pi * li)) - 1.5 * r_li_ratio) ;/// ell; division by ell for binding per unit length, multiplication by ell to account for higher rate of visit of mu-states.
    double rate = exp(1.5 * log(1.5 / (Pi * li)) - 1.5 * r_li_ratio) /ell;
    return rate;
}

void Dangling::compute_cum_rates(vector<double>& sum_l_cum_rates,
                               vector<vector<double>>& cum_rates) const
{
//reserve the correct amount of memory
sum_l_cum_rates.resize(free_linkers.size());
cum_rates.resize(free_linkers.size());
for(auto& cum_rate : cum_rates){cum_rate.resize((int)ell-1);}
// the entry vectors must be empty
int rindex(0);
Linker* currentR = (Rleft != nullptr) ? Rleft : Rright;
for (auto &rlink : free_linkers)
  {
    // iterate through each linker
    // and compute a cumulative binding rate vector for each
    // length and.
    int ellindex(0);
    double squared_diff(get_square_diff(currentR->r(), rlink->r()));
    for (int ELL = 1; ELL < (int)ell; ELL++)
    {
      // add it to the cumulative vector
      if (ellindex == 0)
      {
        cum_rates[rindex][ellindex] = binding_rate_to_integrate(ELL, squared_diff);//compute_binding_rate((double)ELL,rlink);
      }
      else
      {
        cum_rates[rindex][ellindex] = cum_rates[rindex][ellindex-1]+binding_rate_to_integrate(ELL,squared_diff);//compute_binding_rate((double)ELL,rlink);
      }
      ellindex++;
    }
    // Add the end of this vector, which is the total probability
    // to bind to this specific linker to sum_l_cum_rates.
    if(rindex==0){sum_l_cum_rates[0] = cum_rates[rindex].back();}
    else{sum_l_cum_rates[rindex] = sum_l_cum_rates[rindex-1]+cum_rates[rindex].back();}
    // add the whole vector to to cum_rates:
    rindex++;
  }
}
double Dangling::compute_total_rate(Linker* rlinker) const
{
  Linker* currentR = (Rleft != nullptr) ? Rleft : Rright;
  double DIFF(get_square_diff(currentR->r(),rlinker->r()));
  double integral(0);
  for(int i=1; i<ell;i++){integral+=binding_rate_to_integrate(i, DIFF);}
  return integral;
}

pair<unique_ptr<Strand>,unique_ptr<Strand>> Dangling::bind() const
{
  // return a reference to a the two loop that must be created
  // when binding the current loop to a linker randomly choosen
  // in the vicinity of the current loop. The choice is made
  // in accordance with the binding rate of each specific linker
  Linker* linker_selected;
  double length;
  select_link_length(length,linker_selected);
  //linker_selected->set_bounded();
  if (Rleft == nullptr){
    // if the dangling, is a left dangling bond, then invert the right, and left linker and invert the right and left ell_coordinate.
  unique_ptr<Loop> right_loop = make_unique<Loop>(linker_selected,Rright,length,ell,rho0,slide);
  unique_ptr<Dangling> dangling = make_unique<Dangling>(nullptr,linker_selected,0.,length,rho0,slide);
  return  {move(dangling),move(right_loop)};  
  }
  unique_ptr<Loop> left_loop = make_unique<Loop>(Rleft,linker_selected,ell_coordinate_0,ell_coordinate_0+length,rho0,slide);
  unique_ptr<Dangling> dangling = make_unique<Dangling>(linker_selected,nullptr,ell_coordinate_0+length,ell-length,rho0,slide);
  return  {move(left_loop),move(dangling)};
}
unique_ptr<Strand> Dangling::unbind_from(Strand* left_strand) const
{
  // return a reference to the loop that must be created when
  // unbinding a linker between the current loop (at its right)
  // and "left_loop" that is at its left.
  if(Rright == nullptr and left_strand->get_Rleft()==nullptr){throw logic_error("issue with the unbind from dangling");}
  return make_unique<Dangling>(left_strand->get_Rleft(),
                                Rright,
                               left_strand->get_ell_coordinate_0(),
                               ell + left_strand->get_ell(),
                               rho0,
                               slide);
}
unique_ptr<Strand> Dangling::do_slide(double dl,bool right) const
{
  return make_unique<Dangling>(Rleft,Rright,ell_coordinate_0+dl,ell-dl,rho0,slide);
}
double Dangling::get_ell_coordinate_1() const{return ell_coordinate_0+ell;}