#include "Header.h"

using namespace std;

Linker::Linker(std::array<double,3> r_c,int dim,int id){R = r_c;free=true;dimension=dim;ID=id;}
Linker::~Linker(){}

array<double,3> Linker::r() const{
    if(dimension == 3){
        return R;}
    else if(dimension == 2){
        return {R[0],R[1],0.};
    }
    else if(dimension == 1){
        return {R[0],0.,0.};
    }
    else{throw invalid_argument("invalid dimension value");}
}

int Linker::g_ID() const{return ID;}

void Linker::set_free(){free = true;}

void Linker::set_bounded(){free=false;}

bool Linker::is_free() const{return free;}

void Linker::add_strand(Strand* strand){strands.insert(strand);}

void Linker::remove_strand(Strand* strand){strands.erase(strand);}

set<Strand*,LessLoop> Linker::get_strands() const
{
    return strands;
}
void Linker::print_position(string end)const
{
    cout<<R[0]<<" "<<R[1]<<" "<<R[2]<<end;
}
void Linker::diffuse(array<double,3> r)
{
    // chose a direction to make the move.
    if(r[0]!=0 ||  r[1]!=0. || r[2]!=0){
        R = r;
    }
    normal_distribution<double> distribution(0.,1.);
    double dx(distribution(generator)); //generate a diffusion vector
    double dy(distribution(generator)); //generate a diffusion vector
    double dz(distribution(generator)); //generate a diffusion vector
    double norm = sqrt(dx*dx+dy*dy+dz*dz); // compute its norm
    //double x,y,z; // cartesian coordinate of the linker
    //sph2cart(R[0],R[1],R[2],x,y,z); // convert the spherical to cartesian
    R[0]+=dx/norm; // add the diffusion vector
    R[1]+=dy/norm; // add the diffusion vector
    R[2]+=dz/norm; // add the diffusion vector
    //cart2sph(x,y,z,R[0],R[1],R[2]); // con
}

