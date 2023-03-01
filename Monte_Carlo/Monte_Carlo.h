#ifndef Monte_Carlo_h
#define Monte_Carlo_h
class Monte_Carlo
{
    public:
        Monte_Carlo(double ell_tot, // total length of the polymer
                    double BindingEnergy, // binding energy...
                    double percent_diff_move, // percentage of move dedicated to diffusion
                    int seed, 
                    int Nlinker_max, // total number of linker in the system
                    int dimension); // dimension of the system
        ~Monte_Carlo();
        void Evolve(bool* accepted, int* move_type, double* DE);
    private:
        double ell,binding_energy,percent_diff,rho;
        int see, Nlinker;
        bool slide;
        std::uniform_int_distribution<int> distrib;
        LoopLinkWrap loop_link;


        std::set<std::array<double,3>> generate_crosslinkers(bool remake);
};
#endif