#include "Header.h"

using namespace std;

Linker::Linker(std::array<double,3> r_c,int dim,int id){R = r_c;deltaR={0,0,0};free=true;dimension=dim;ID=id;}
Linker::~Linker(){}

array<double,3> Linker::r(bool periodic) const{
    array<double,3> res = R;
    if(periodic){
        for(int i = 0;i<3;i++){
            res[i]-=deltaR[i];
            }
    }
    switch (dimension ){
        case 3 :
            return res;
        case 2:
            return {res[0],res[1],0.};
        case 1 :
            return {res[0],0.,0.};
        default :
            throw invalid_argument("invalid dimension value");
    }    
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
            //cout<<r[0]<<" "<<r[1]<<" "<<r[2]<<endl;
            //cout<<deltaR[0]<<" "<<deltaR[1]<<" "<<deltaR[2]<<endl;
            //cout<<R[0]<<" "<<R[1]<<" "<<R[2]<<endl;
        for(int i = 0; i<3;i++){
            deltaR[i] += r[i] - R[i];
        }
        R = r;
        return;
            //cout<<endl;
            //cout<<r[0]<<" "<<r[1]<<" "<<r[2]<<endl;
            //cout<<deltaR[0]<<" "<<deltaR[1]<<" "<<deltaR[2]<<endl;
            //cout<<R[0]<<" "<<R[1]<<" "<<R[2]<<endl;
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

