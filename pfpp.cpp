//! @file pfpp.cpp
//! @brief All particle filter code

//! With this definition, we turn off Eigen library debugging code.
#define EIGEN_NO_DEBUG
#include <Eigen/Eigen>

// for Kronecker product:
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
namespace ublas = boost::numeric::ublas;

#include <boost/math/distributions/geometric.hpp>
#include <boost/math/distributions/normal.hpp>
#include <boost/math/distributions/cauchy.hpp>

// For portable, replicable distributions:
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/cauchy_distribution.hpp>


#include <boost/uuid/uuid.hpp>
#include <boost/uuid/sha1.hpp>

#include "prettyprint.hpp"

#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <utility>
#include <algorithm>
#include <array>
#include <vector>
#include <map>
#include <set>
#include <tuple>
#include <chrono>
#include <random>

using namespace std;
using namespace Eigen;


#include "country_vax_efficacy.hpp"// This includes the map of ISO3 codes to vaccine efficacy category

//! This flag should be defined whenever using data from 2018 or after
#define DATA2018

//! This definition will produce replicable random numbers, useful for testing and development.
#define STATICSEED
// You do not want to run with the static seed when producing replicate runs
// for the pf_summarizer.  In fact, presently, the code will exit if you try to do so.


//! If defined, print I and S compartments for particles from best-fit parameter set
#define PRINT_PARTICLES

//! The following will output compartments by ages
#define PRINTAGES
// Specifically, it will generate two files:
// ISOmodel_iap_particles_fxsave.txt
// ISOmodel_sap_particles_fxsave.txt
// Each file will contain space-delimited counts for their respective compartments
// for each age, with one particle per line

// This will output the "long" format which is needed for the pf_summarizer R code.
#define PRINTAGES_BYYEAR

// The following will output compartments, dividing ages above and below 5:
// #define PRINTAGESSPLIT

// The following defines will activate additional logging and asserts:

//#define DEBUG_PARTICLES
//#define DEBUG_PARTICLES_FX
//#define DEBUGWEIGHTS
//#define DEBUG_PARTICLES_DEMOG

#define SUSCEPTIBLE_METHOD_2

constexpr double theta_value = 0.05; 

constexpr int gridsize = 2000; //!< Number of cells in parameter grid, can be altered
constexpr int npart = 1000; //!< Total number of particles to simulate, can be altered.
constexpr int nyrs = 3;  //!< Number of years stored in memory, should not be edited.
constexpr int calendaryrs = 43;//!< Number of years to be simulated, should always correspond to demographic data supplied.

constexpr int nages = 100;//!< Number of age classes to be simulated. Should not be edited.

constexpr int nagesDecremented = nages-1; //!< Calculated at compile time for efficiency.
constexpr int SIAcalendaryrs = calendaryrs + 1;//!< Calculated at compile time for efficiency.

//! The initial version of this project allocated S, I, and N compartments for each particle.
//! However, the N compartments vary with time, but not between particles.  As such, they do not
//! need to be tracked on a per-particle basis.  By default, we do not allocate N for each particle.
//! Prevents allocating redundant storage of N compartment data.
#define DROPN
#ifdef DROPN
constexpr int npars = 2;//!< compartments: SI
#else
constexpr int npars = 3;//!< compartments: SIN
#endif


#ifdef STATICSEED
seed_seq my_seed_seq {42};
#endif

#ifndef STATICSEED
// Generate a random seed, as per discussion in: https://stackoverflow.com/a/34490647
random_device r;
seed_seq my_seed_seq{r(), r(), r(), r(), r(), r()};
#endif

// Create Mersenne Twister generator
mt19937_64 generator{my_seed_seq};

// For portable Normal distributions:
boost::random::mt19937 rng{my_seed_seq};

std::string get_sha1(const std::string& p_arg)
{
    boost::uuids::detail::sha1 sha1;
    sha1.process_bytes(p_arg.data(), p_arg.size());
    unsigned hash[5] = {0};
    sha1.get_digest(hash);

    // Back to string
    char buf[41] = {0};

    for (int i = 0; i < 5; i++)
    {
        std::sprintf(buf + (i << 3), "%08x", hash[i]);
    }

    return std::string(buf);
}


std::string get_timestamp(){

    //get current time in format of time_t
    time_t t = time(NULL);

    //convert time_t to char*
    //char* charTime = ctime(&t);
    //display current time
    //cout << "Now is " << charTime << endl;

    //convert time_t to tm
    //tm* my_time = localtime(&t);

    std::stringstream buffer;

    //buffer << std::put_time(localtime(&t), "%y%m%d%H%M"); // 2207011148
    buffer << std::put_time(localtime(&t), "%y%m%d"); // 220701

cout << buffer.str() << endl;
return(buffer.str());
}

//! This struct holds particle information
struct ptuple
{
    int index;//!< Index of this particle in vector
    double pf;//!< Sample weight
    double cf;//!< Run probability
    int count;//!< Number of times this particle is selected to propagate
};


//! A struct to hold and pass command line arguments.
struct cli_inputs
{
    int action;//!< Action: Grid Search = "1", Single Parameter Set = "5"
    string country_model;//!< ISO3 code with model as suffix, e.g.: "NGAmcv2"
    string demofile;//!< Demographic input file
    string siafile;//!< Input file for SIAs
    int modelnumber;//!< 0 = "normal"; 1 = "discounting mcv2"; 2 = "discounting SIAs"; 3 = "discounting MCV2 and SIAs"
    double sigma;//!< Noise term
    string iestfile;//!< Optional, if supplied contains a single parameter set for a single run
    bool fixed_prhigh;//!< If true, set probability of reporting to "high" for all years
    double prhigh;//!< Optional, if set will determine the "high" reporting rate
    double cauchy_weight;//!< Term which scales contribution to observation weight from Cauchy distribution
    string timestamp;//!< Run timestamp
    string fullcl;//!< String which contains full command line 
    cli_inputs()
    {
        action=0;
        country_model=string("");
        demofile=string("");
        siafile=string("");
        modelnumber=0;
        sigma=0.0;
        iestfile=string("");
        fixed_prhigh = false;
        prhigh = -20;
        cauchy_weight = 0.0;
    }
};

//! Simply packages parameters and their observed negative log-likelihood for output
struct Iest_struct
{
    double nLL;//!< Negative Log Likelihood
    double S0pct;//!< Noise
    VectorXd paramsBest;//!< b0, b1, pr, prhi
};

//! Store values returned from fx function
struct fx_return_struct
{
    double nLL;//!< Negative Log Likelihood
    double lik;//!<
    double nLLt;//!<
    double depletedyearsum;//!<
    double penalty;//!<
    double S0pct;//!< Noise
    VectorXd paramsBest;//!< b0, b1, pr, prhi
};


//! Writes best nLL and associated parameters to disk
void writeIest(Iest_struct Iest,//!< Iest_struct to write to disk
               const string& countryISO3,//!< country IS03 + model
               int single,//!< Flag indicating whether this Iest_struct resulted from a single parameter set run
               const string& ttt//<! sha hash + timestamp
              )
{
    std::stringstream fc;

    fc << Iest.nLL << "\n";
    fc << Iest.S0pct << "\n";
    fc << Iest.paramsBest << "\n";

    std::stringstream a_fn_ss;

    a_fn_ss << "Iest_" ;
    if(single) {
        // If this not the result from a grid search, indicate that in filename:
        a_fn_ss << "single_";
    }
    a_fn_ss << countryISO3 << ttt << ".txt";
    std::string a_fn_str = a_fn_ss.str();
    std::ofstream a_ofst(a_fn_str.c_str());
    //df.iest_output_sha = get_sha1(fc.str());
    if (a_ofst.is_open())
    {
        a_ofst << fc.str();
    }
    else
    {
        std::cerr << "Could not open Iest file for writing" << std::endl;
    }
}


//! Read Iest file
// File will contain six lines, and look like this:
// 7941.98 // <-- nLL
// 0.100042 // <-- S0pct
//   -4.00622 // <-- Beta0
//      88.44 // <-- Beta1
// 0.00625985 // probabilty of reporting
//  0.0618812 // probabilty of reporting (high)
Iest_struct readIest(const string& iestfile)
{
    Iest_struct Iest;
    VectorXd pb(4);
    Iest.paramsBest = pb;
    std::ifstream intest(iestfile);
    std::string line;

    std::getline(intest,line);
    double dnLL =  stod(line);
    Iest.nLL = dnLL;

    std::getline(intest,line);

    Iest.S0pct = stod(line);

    std::getline(intest,line);
    Iest.paramsBest(0) = stod(line);

    std::getline(intest,line);
    Iest.paramsBest(1) = stod(line);

    std::getline(intest,line);
    Iest.paramsBest(2) = stod(line);

    std::getline(intest,line);
    Iest.paramsBest(3) = stod(line);

    return(Iest);
}

//! This struct holds all data from the demographic input file
struct demographicDataFrame
{
    vector<int> C;//!< Reported cases
    vector<int> N;//!< Population size
    vector<int> B;//!< Births
    vector<int> highyears;//!< Years with higher rates of reporting
    vector<double> mcv1;//!< MCV1 coverage
    vector<double> mcv2;//!< MCV2 coverage
    vector<double> SIA;//!< SIA coverage
    double SIA_matrix [nages][SIAcalendaryrs];//!< Age-specific SIA coverage
    string countryISO3;// ISO3 code for country
    string countryModel;//!< IS03 code + model
    int modelnumber;
    int action;
    // shasums and timestamps?
    string timestamp;//!<  timestamp
    string sha;//!< SHA + timestamp
    string shatimestamp;//!< SHA + timestamp
    string cli_args;
 string dm_sha;//!< post-modification demographic data 
 string sia_sha;//!< post-modification sias 
 string pm_dm_sha;//!< post-modification demographic data 
 string pm_sia_sha;//!< post-modification sias 
string s_particles_sha;//<! SHA for 
string i_particles_sha;//<! SHA for 
#ifdef PRINTAGES
string sap_particles_sha;//<! SHA for 
string iap_particles_sha;//<! SHA for 
#endif
string runs_showing_penalty_sha;//<! SHA for
string Iest_sha; //<! SHA for 
float  nLL;
float  S0pct;
float  p_b0;
float  p_b1;
float  p_pr;
float  p_pr_high;
};


//! Print demographic data to console
void printDemographicDataFrameToConsole(const string& hstring,//!< String prepending to each line of output
                                        const demographicDataFrame &df)//!< Demographic data frame
{
    int yrs = df.C.size();
    for(int t=0; t<yrs; t++)
    {
        cout << hstring << " " << t << " C: " << df.C[t] << " N: " << df.N[t]<< "  mcv1: " << df.mcv1[t] << " mcv2: " << df.mcv2[t] << " SIA: " << df.SIA[t]  << endl;
    }
}

void writeDemographicDataFrameToInfoFile(const demographicDataFrame &df, const cli_inputs &cli){
    
    stringstream ss;
ss << "command: " << cli.fullcl << endl;

 ss << "modelNumber: " << df.modelnumber << endl;
    ss << "action: " << cli.action << endl;
ss << "countryISO3: " << df.countryISO3 << endl;
 ss << "countryModel: " << df.countryModel << endl;
    ss << "demofile: " << cli.demofile << endl;
    ss << "siafile: " << cli.siafile << endl;
    ss << "sigma: " << cli.sigma << endl;
    ss << "cauchy_weight: " << cli.cauchy_weight << endl;
 ss << "theta_value: " << theta_value << endl;
 ss << "gridsize: " << gridsize << endl;
 ss << "npart: " << npart << endl;


 ss << "calendaryrs: " << calendaryrs << endl;
 ss << "ages: " << nages << endl;

#ifdef STATICSEED
    ss << "STATICSEED: TRUE" << "\n";
#else
    ss << "STATICSEED: FALSE" << "\n";
#endif

#ifdef SUSCEPTIBLE_METHOD_1 
    ss << "SUSCEPTIBLE_METHOD: 1" << "\n";
#endif


#ifdef SUSCEPTIBLE_METHOD_2 
    ss << "SUSCEPTIBLE_METHOD: 2" << "\n";
#endif
 // give shas more explanatory names
 // add sha for Iest input if action #5  Or not needed?
 ss << "demofile_siafile_sha: " << df.sha << endl;
 ss << "sha_timestamp: " << df.shatimestamp << endl;
// // add shas for other outputs, if not too slow
 //
 
 if (cli.action == 5){
    ss << "input_iestfile: " << cli.iestfile << endl;
 }

 ss << "demofile_sha: " << df.dm_sha << endl;
 ss << "siafile_sha: " << df.sia_sha << endl;
 ss << "modified_demofile_sha: " << df.pm_dm_sha << endl;
 ss << "modified_siafile_sha: " << df.pm_sia_sha << endl;
 ss << "i_particles_sha: " << df.i_particles_sha << endl;
 ss << "s_particles_sha: " << df.s_particles_sha << endl;

#ifdef PRINTAGES
 ss << "iap_particles_sha: " << df.iap_particles_sha << endl;
 ss << "sap_particles_sha: " << df.sap_particles_sha << endl;
#endif
 
 if (cli.action == 1){
    ss << "runs_with_penalty: " << df.runs_showing_penalty_sha << endl;
ss << "nLL: " << df.nLL << endl;
ss << "S0pct: " << df.S0pct << endl;
ss << "b0: " << df.p_b0 << endl;
ss << "b1: " << df. p_b1 << endl;
ss << "pr: " << df.p_pr << endl;
ss << "pr_high: " << df.p_pr_high << endl;
 }

    std::stringstream fnss;
    fnss << df.countryModel;

 if (cli.action == 5){
        // If this not the result from a grid search, indicate that in filename:
        fnss << "_single";
    }
       fnss << "_info" << df.shatimestamp  << ".txt";
    std::string fns= fnss.str();
    std::ofstream a_xofstd(fns.c_str());
    if (a_xofstd.is_open())
    {
        a_xofstd << ss.str();
    }
    else
    {
        std::cerr << "Could not open " << fns << std::endl;
    }

}

//! Print age-class specific SIA coverages to console
//void printDataFrameSIAs(const demographicDataFrame& df)
void printDataFrameSIAs(demographicDataFrame& df)
{
    stringstream sss;
    for(int v = 1; v<calendaryrs; v++)
    {
        sss << "\"V" << v << "\",";
    }

    sss << endl;
    for(int j=0; j< nages; j++)
    {
        sss << "\"" << (j+1) << "\"," ;
        for(unsigned i=0; i< df.SIA.size(); i++)
        {
            //sss << i << "| " << df.SIA[i] << "| ";
            sss << df.SIA_matrix[j][i+1] << ",";
        }
        sss << "\n";
    }

    std::stringstream a_fn_ss_dbg;
    a_fn_ss_dbg << "sia" << df.countryModel << "_post-modification" << df.shatimestamp << ".txt";
    std::string a_fn_strd = a_fn_ss_dbg.str();
    std::ofstream a_ofstd(a_fn_strd.c_str());
    if (a_ofstd.is_open())
    {
        a_ofstd << sss.str();
    }
    else
    {
        std::cerr << "Could not open " << a_fn_strd << std::endl;
    }

    df.pm_sia_sha = get_sha1(sss.str());


    stringstream ssx;
    ssx <<"\"C\" \"N\" \"B\" \"mcv1.true\" \"mcv2.true\" \"SIA\" \"high.years\"" << endl;
    int yrs = df.C.size();
    for(int t=0; t<yrs; t++)
    {
        ssx << df.C[t] << " " << df.N[t]<< " " << df.B[t] << " " << df.mcv1[t] << " " << df.mcv2[t] << " " << df.SIA[t]  << " " << df.highyears[t] << endl;
    }


    std::stringstream a_fn_ss_xdbg;
    a_fn_ss_xdbg << df.countryModel << "_post-modification" << df.shatimestamp << ".txt";
    std::string a_fn_xstrd = a_fn_ss_xdbg.str();
    std::ofstream a_xofstd(a_fn_xstrd.c_str());
    if (a_xofstd.is_open())
    {
        a_xofstd << ssx.str();
    }
    else
    {
        std::cerr << "Could not open " << a_fn_xstrd << std::endl;
    }

df.pm_dm_sha = get_sha1(ssx.str());

}



//! Read demographic input file into data structure.
// Input should have the following format:
// "C" "N" "B" "mcv1.true" "mcv2.true" "SIA" "high.years"
// 943 39619603 444411 0 1 0 1
// 1140 38242022 429178 0.213 0 1 1
// 1214 36929085 415035 0.212 0 1 1
// [...]
// 33 16307382 182324 0.867 0.5785955069 0 1
// 72 15849443 177567 0.887 0.5748005761 0 1
// 635 15403406 172314 0.895 0.5844050531 0 1
std::string readDemographicDataFrame(std::istream& str,//!< Path to demographic input file.
                              demographicDataFrame& df//!< Reference to demographic data frame.
                             )
{
    std::string line;
    std::stringstream ffile;
    while(std::getline(str,line))
    {
        cout << line << endl;
        ffile << line << endl;
        std::stringstream          lineStream(line);
        std::string                cell;
        std::vector<std::string>   result;

        while(std::getline(lineStream,cell,' '))
        {
            result.push_back(cell);
        }

        if(result[0] != "\"C\"")
        {
            int cases = -1;
            if(result[0] != "NA")
            {
                cases = stoi(result[0]);
            }
            df.C.push_back(cases);
            df.N.push_back(stoi(result[1]));
            df.B.push_back(stoi(result[2]));
            df.mcv1.push_back(stod(result[3]));
            df.mcv2.push_back(stod(result[4]));
            df.SIA.push_back(stod(result[5]));
            int hy = 0;
            if(result[6] != "NA")
            {
                hy = stoi(result[6]);
            }
            df.highyears.push_back(hy);
        }
        // handle sha
//        unsigned char ta[] = lineStream.str().c_str();
//        char* c = const_cast<char*>(lineStream.str().c_str());
        //cmakeSHA1boostV1(ta);
    }
        //cout << get_sha1((ffile.str())
        return(ffile.str());
}


//! Read SIA input file into data structure.
// Input should have the following format:
// "V1","V2","V3","V4","V5",[...],"V36","V37","V38","V39","V40"
// "1",0,0,0,0,0,[...],0.187,0.142375,0.210375,0,0.18275
// "2",0,0,0,0,0,[...],0.4092,0.6231,0.9207,0,0.7998
// [...]
// "99",0,0,0,0,0,0,[...],0,0,0,0,0
// "100",0,0,0,0,0,0,[...],0,0,0,0,0
std::string readAgeSpecificSIAs(std::istream& str, demographicDataFrame& df )
{
    std::string line;
    int x = 0;
    int lx = 0;
    std::stringstream ffile;
    while(std::getline(str,line))
    {
        //cout << line << endl;
        ffile << line << endl;
        std::stringstream          lineStream(line);
        std::string                cell;
        std::vector<std::string>   result;

        cout << line << endl;
        //  while(std::getline(lineStream,cell,' '))
        while(std::getline(lineStream,cell,','))
        {
            result.push_back(cell);
        }

#ifdef DATA2018
        if(result[0] != "\"V1\"")
        {
#else
        if(result[1] != "\"V1\"")
        {
#endif
            for(unsigned i=0; i<(calendaryrs+1); i++)
            {
                if(result[i] != "0")
                {
                    double sv = 0.0;
                    if(i>0)
                    {
                        sv = stod(result[i]);
                    }
                    if(sv>0.0)
                    {
                        df.SIA_matrix[x][i] = sv;
                    }
                }
            }
            ++x;
        }
        ++lx;
    }

        return(ffile.str());
}


//! Calculate likelihood (i.e. weight) from the observation model
//! @param[out] VectorXd A vector of particle weights
VectorXd obswtIt(const int t,//!< Timestep
                 const MatrixXd& statesIt,//!< "I" compartments for all particles
                 const MatrixXd& params,//!< Parameter set
                 const vector<int>& observation, //!< The observed cases column ("C"), from the country input file
                 int yradj,//!< If this is a high-reporting year, value will be 1
                 const double cauchy_weight//!< Mixture weight on Cauchy contribution
                )
{
    // obs.wt <- function(t,states, params, observation, particle, yr.adj){
    // # calculate likelihood (i.e. weight) from the observation model
    // #    if(is.na(yr.adj)){
    //     p.obs <- min(.99, params[3])
    //
    double pobs = 0.99;

    // If this is a "high-reporting" year, scale probability of  observation
    if(yradj==1)
    {
        pobs = params(2);
    }
    else
    {
        pobs = params(3);
    }


    // wt <-
    // log((1-cauchy.weight)*dnorm( observation, mean = p.obs*It, sd = sqrt(pmax(It,.1) * p.obs * (1 - p.obs)),log=F) +
    // cauchy.weight*dcauchy( observation, location = p.obs*It, scale = sqrt(pmax(It,.1) * p.obs * (1 - p.obs)),log=F)
    // +1e-108)

    VectorXd myIt(npart);
    myIt = statesIt.colwise().sum();

    VectorXd mymean(npart);
    for(unsigned i=0; i<npart; i++)
    {
        mymean[i] = myIt(i) * pobs;
    }

    VectorXd mysd(npart);
    double pobsm1 = 1.0 - pobs;
    double pobs_pobsm1 = pobs * pobsm1;
    for(unsigned i=0; i<npart; i++)
    {
        double myIt_i = 0.1;
        if(myIt(i) > myIt_i)
        {
            myIt_i = myIt(i);
        }

        mysd[i] = sqrt( myIt_i * pobs_pobsm1); //sqrt( statesIt[t][i] * pobs * (1 - pobs));
    }

    VectorXd pwt(npart);
    VectorXd cwt(npart);
    VectorXd wt(npart);
    pwt.setZero();
    wt.setZero();

    VectorXd forsort_wt(npart);
    forsort_wt.setZero();

    for(int i=0; i<npart; i++)
    {
        if(mymean[i]>0)
        {
            boost::math::normal_distribution<double> dnorm(
                std::lround(mymean[i]),// mean
                mysd[i]// standard deviation
            );
            pwt[i] = pdf(dnorm, observation[t]);

            boost::math::cauchy_distribution<double> dcauchy(
                std::lround(mymean[i]),// location
                mysd[i]// scale
            );

            cwt[i] = pdf(dcauchy, observation[t]);
        }
        else
        {
            pwt[i] =  0.0;
            cwt[i] =  0.0;
        }
    }

    // Calculate final weight for each particle
    for(int i=0; i<npart; ++i)
    {
        wt[i] = log((1.0-cauchy_weight)*pwt[i] + cauchy_weight*cwt[i]+1e-108);
    }

    // Return vector of particle weights
    return(wt);

}


//! Advance the simulation model one time step.
VectorXd simstep(
    const unsigned t, //!< Number of years from 1980
    const unsigned ft,//!< Index of current year's data in states
    const unsigned bft,//!< Index to prior year's data in states
    double (*states)[nyrs][npars][nages][npart],//!< Pointer to heap-allocated states
    const MatrixXd& params,//!< Parameters to be evaluated
    const double sigma,//!< Noise term
    const double theta,//!< Theta term
    const demographicDataFrame df,//!< Demographic dataframe
    const double efficacy//!<
)
{
    // # Advance the simulation model one time step
    // # dimensions are times, states, ages, particles

    // St <- states[t-1,1,,] # last years
    // It <- states[t-1,2,,]
    // Nt <- states[t,3,,] # current year

    MatrixXd St(nages,npart);
    for(int age=0; age<nages; age++)
    {
        for(int p=0; p<npart; p++)
        {
            St(age,p) = states[0][bft][0][age][p];
#ifdef DEBUG_PARTICLES_DEMOG
            if(isnan( states[0][bft][0][age][p]))
            {
                cout << "isnan t:" << t << " age: " << age << " p:" << p << endl;
                assert(0);
            }
#endif
        }
    }

    // It <- states[t-1,2,,]
    MatrixXd It(nages,npart);
    for(int age=0; age<nages; age++)
    {
        for(int p=0; p<npart; p++)
        {
            It(age,p) = states[0][bft][1][age][p];
        }
    }

    if(t==1)
    {
        It.setZero();
    }

#ifndef DROPN
    // Nt <- states[t,3,,] # current year
    MatrixXd Nt(nages,npart);
    for(int age=0; age<nages; age++)
    {
        for(int p=0; p<npart; p++)
        {
            //Nt(age,p) = *states[ft][2][age][p];
            Nt(age,p) = states[0][ft][2][age][p];

        }
    }
#endif


    // St <- S.update(t, St, births, n.ages, mcv1, mcv2, SIA, sigma)
    //  Supdate(t, St, sigma, df, efficacy,particle);
    // *************** SUPDATE ******************************

    double ub = df.B[t]*(1.0-df.mcv1[t]*efficacy);

    // St <- rbind(births[t]*(1-mcv1[t]*.85), St[-n.ages,])  # Rbind(unvac births,  all but oldest group)
    for(int i=(nages-1); i>0; i--)
    {
        for(int j=0; j<npart; j++)
        {
            St(i,j) = St(i-1,j);
        }
    }
    for(int i=0; i<npart; i++)
    {
        St(0,i)=ub;
    }
    // St[3,] <- St[3,]  * (1 - mcv2[t]*.95) # second dose efficacy is 95%
    double second_dose = (1.0 - df.mcv2[t]*0.95);
    for(int j=0; j<npart; j++)
    {
        //St(2,j) = St(2,j) * (1.0 - df.mcv2[t]*0.95);
        St(2,j) = St(2,j) * second_dose;
    }

    if(df.SIA[t] > 0.0)   // This is the SIA in the ISO3.txt file; is there a SIA this year?
    {
        // If there is a SIA this year, we access the age-specific SIA coverages:
        VectorXd csia(nages);
        VectorXd mcsia(nages);
        for(int j=0; j<nages; j++)
        {
            csia(j) = df.SIA_matrix[j][t+1];
            mcsia(j) = 1.0 - df.SIA_matrix[j][t+1];
        }

        //     sia <- SIA[,t]
        //     St <- (St) *(1 - sia)
        for(int age=0; age<nages; age++)
        {
            if(csia[age]>0.0)   // If SIA value for this age from siaISO.txt file
            {
                for(int aparticle=0; aparticle<npart; aparticle++)
                {
                    double nv = St(age,aparticle) * mcsia(age);
                    if(nv<0.0)
                    {
                        nv=0.0;
                    }
                    St(age,aparticle) = nv;
                }
            }
        }
    }
    // End of S update

    // St.mat <- St    # Rbind(unvac births,  all but oldest group)    # MATT ADDED ON 11 MAR
    MatrixXd Stmat(nages,npart);
    Stmat = St;

    // St.mat[1,] <- floor(St.mat[1,]*0.25)  # MATT ADDED Discounting first age class to account for maternal immunity
    for(int p=0; p<npart; p++)
    {
        Stmat(0,p) = floor(  Stmat(0,p) * 0.25 );
    }

    // S.prop <- apply(St.mat,2,sum) / Nt[1,]  # all susceptibles; not just age group's

    VectorXd Spropp(npart);
    Spropp.setZero();

    for(int p=0; p<npart; p++)
    {
        std::vector<double> fs(nages, 0.0);
        for(int a=0; a<nages; a++)
        {
            fs[a] =  Stmat(a,p);
        }
        std::sort (fs.begin(),fs.end());

        Spropp[p] = std::accumulate(fs.begin(), fs.end(), 0.0); //Stmat(a,p);
#ifdef DEBUG_PARTICLES_DEMOG
        if(isnan(Spropp[p]))
        {
            cout << "ISNAN: time: " << t << " Spropp:" << Spropp[p] << " particle: " << p << endl;
            assert(0);
        }
#endif
    }

    vector<double> Sprop;
    Sprop.reserve(npart);

#ifdef DROPN
    double dfnt = df.N[t];
#else
    VectorXd Ntrow = Nt.row(0.0);
#endif
    for(int p=0; p<npart; p++)
    {
#ifndef DROPN
        Sprop.push_back(Spropp[p]/Ntrow[p]);
#else
        Sprop.push_back(Spropp[p]/dfnt);
#endif

    }


    // *************** IUPDATE ******************************

    double b0 = params(0);
    double b1 = params(1);

    It.setZero();

    double forexp[npart];

    //     noise <- rnorm(particle, mean = 0, sd = sigma)
    //     p.inf <- rep(exp(b0+b1*S.prop+noise)/(1 + exp(b0+b1*S.prop+noise)),25)

    boost::random::normal_distribution<double> distribution(0,sigma);

    vector<double> store_sigmanoise;
    store_sigmanoise.reserve(npart);

    for(int i=0; i<npart; i++)
    {
        double sigmanoise =  distribution(rng);

        forexp[i] =  b0 + b1 * Sprop[i] + sigmanoise;

        store_sigmanoise.push_back(sigmanoise);

#ifdef DEBUG_PARTICLES_DEMOG
        if(isnan(forexp[i]))
        {
            cout << "ISNAN: time: " << t << " " << Sprop[i] << " " << b0 << " " << b1 << " particle: " << i
                 << " sigmanoise: " << sigmanoise << " forexp " << forexp[i] << endl;
            assert(0);
        }
#endif

    }

    VectorXd xpinf(npart);
    for(int i=0; i<npart; i++)
    {
        double ef = exp(forexp[i]);
        xpinf[i] = ef/(1.0+ef);
    }


    //  p.inf <- pmax(p.inf + rnorm(length(p.inf),0,pmax(theta - S.prop*(theta/.05),0)),0)  # adds noise to pinf when less than 0.05
    for(int i=0; i<npart; i++)
    {
        double rnorm_sd = theta - (Sprop[i] * (theta/0.05));
        if(rnorm_sd < 0.0)
        {
            rnorm_sd = 0.0;
        }


        boost::random::normal_distribution<double> noise_distribution(0.0,rnorm_sd);
        double potentialnoise_preabs = noise_distribution(rng);

        double potentialnoise = abs(potentialnoise_preabs);

#ifdef DEBUG_PARTICLES_DEMOG
        if(isnan(potentialnoise))
        {
            cout << "ISNAN:  pinf " << xpinf[i] << " sd: " << rnorm_sd <<  " potentialnoise: " << potentialnoise << endl;
            assert(0);
        }

        if(isnan(xpinf[i]))
        {
            cout << __LINE__ << " " << i << " pinf " << xpinf[i] << " sd: " << rnorm_sd <<  " potentialnoise: " << potentialnoise << endl;
            assert(0);
        }
#endif

        // intent is for potential noise to decrease possible number of infecteds
        double potentialpinf = xpinf[i] - potentialnoise;

        if(potentialpinf>0.0)
        {
            xpinf[i] = potentialpinf;
        }
        else
        {
            xpinf[i] = 0.0;
        }

    }


#ifdef DEBUG_PARTICLES_DEMOG
    stringstream pdbg_particles;
    pdbg_particles << "S N SN pinf I forexp Sprop sigmanoise" << endl;
    for(int particle=0; particle<npart; particle++)
    {
        double sfullpar = 0.0;
        double ifullpar = 0.0;
        double nfullpar = 0.0;
        for(int par=0; par<npars; par++)
        {
            float fullpar = 0.0;
            for(int age=0; age<nages; age++)
            {
                fullpar += states[0][ft][par][age][particle];
            }
            if(par==0)
            {
                sfullpar = fullpar;
            }
            else if (par == 1)
            {
                ifullpar = fullpar;
            }
            else
            {
                nfullpar = fullpar;
            }
        }
        pdbg_particles << sfullpar << " " << nfullpar << " " << (sfullpar/nfullpar) << " " << xpinf[particle] << " " << ifullpar << " " << forexp[particle] << " " << Sprop[particle] << " " << store_sigmanoise[particle] << endl;
    }

    pdbg_particles << endl;
    std::stringstream a_fn_ss_pdbg;
    a_fn_ss_pdbg << df.countryModel << "_pdbg_particles_fx_" << t << df.shatimestamp << ".txt";
    std::string a_fn_strd = a_fn_ss_pdbg.str();
    std::ofstream a_ofstd(a_fn_strd.c_str());
    if (a_ofstd.is_open())
    {
        a_ofstd << pdbg_particles.str();
    }
    else
    {
        std::cerr << "Could not open " << a_fn_strd << std::endl;
    }
#endif


    //  # It is now a matrix of size n.ages X particles

    MatrixXd mxpinf(nages,npart);
    MatrixXd m1mxpinf(nages,npart);
    for(int j=0; j<npart; j++)
    {
        double pval = xpinf[j];
        double pvalm1 = 1.0 - pval;
        for(int i=0; i<nages; i++)
        {
            mxpinf(i,j) =  pval;//xpinf[j];
            m1mxpinf(i,j) =  pvalm1;//1.0 - xpinf[j];
        }
    }

    //     St.mat <- St
    //     St.mat[1,] <- floor(St.mat[1,]*0.25)

    Stmat = St;
    for(int p=0; p<npart; p++)
    {
        //Discounting first age class to account for maternal immunity
        Stmat(0,p) = floor(  Stmat(0,p) * 0.25 );
    }

    //    mean.mat =p.inf*St.mat
    MatrixXd meanmat(nages,npart);
    meanmat = mxpinf.cwiseProduct(Stmat);

    // var.mat =(1-p.inf)*p.inf*St.mat
    MatrixXd varmat(nages,npart);

    // meanmat, varmat are nagesXparticles
    varmat = m1mxpinf.cwiseProduct(meanmat);

    // It <- matrix(floor(
    //		      rnorm(particle*n.ages,
    //                          mean = mean.mat,
    //                          sd = sqrt(var.mat)
    //			    )
    // ),n.ages,particle)//

    //      rnorm(particle*n.ages,mean = mean.mat,sd = sqrt(var.mat) # -> 250
    //			    )
    //     It <- pmax(It,0)

    MatrixXd sqrtvarmat(nages,npart);
    for(int p=0; p<npart; p++)
    {
        for(int age=0; age<nages; age++)
        {
            sqrtvarmat(age,p) = sqrt(varmat(age,p));
            if(isnan(sqrtvarmat(age,p)))
            {
                cout << __LINE__ << " " << age << " " << p << " " << meanmat(age,p) << " " << m1mxpinf(age,p) << " " << varmat(age,p) << " " << sqrtvarmat(age,p) << endl;
                assert(0);
            }
        }
    }

    // R col-first:
    for(int p=0; p<npart; p++)
    {
        for(int age=0; age<nages; age++)
        {

            boost::random::normal_distribution<double> idistribution(meanmat(age,p),sqrtvarmat(age,p));
            double pp = floor(idistribution(rng));

            if(pp<0)
            {
                pp=0;
            }

            It(age,p) = pp;

#ifdef DEBUG_PARTICLES_DEMOG
            if(isnan(It(age,p)))
            {
                cout << __LINE__ << " " << age << " " << p << " " << meanmat(age,p) << " " << sqrtvarmat(age,p) << " " << varmat(age,p) << endl;
                assert(0);
            }
#endif
        }
    }
    // *************** END IUPDATE ******************************


    MatrixXd StmItSt(nages,npart);
    StmItSt = St - It;

    for(int a=0; a<nages; a++)
    {
        for(int p=0; p<npart; p++)
        {
            StmItSt(a,p) = floor(StmItSt(a,p));
            // #NOTE -- ensure that each susceptible age class has at least 1 susceptible
            if(StmItSt(a,p)<1.0)
            {
                StmItSt(a,p) = 1.0;
            }
        }
    }

    for(int a=0; a<nages; a++)
    {
        for(int p=0; p<npart; p++)
        {
            states[0][ft][0][a][p] =  StmItSt(a,p); // Assign S - I back into S compartment in states
        }
    }


    for(int a=0; a<nages; a++)
    {
        for(int p=0; p<npart; p++)
        {
            states[0][ft][1][a][p] = It(a,p);// Assign I compartment back into states
        }
    }


    return(xpinf);
}



//! Runs full simulation and particle filter for one set of parameters, with additional outputs and optional debugging code
fx_return_struct fx(
    double (*states)[nyrs][npars][nages][npart],//!< Pointer to heap-allocated states storage
    VectorXd params, //!< Parameters to be evaluated
    const double *S0,//!< Pointer to an array of intial S0 values
    //const demographicDataFrame df,//!< Demographic dataframe
    demographicDataFrame& df,//!< Demographic dataframe
    const string& countryISO3, //!< countryISO3 ISO3 code
    const double sigma, //!< Noise term
    const double cauchy_weight,//!<Weight with which to apply Cauchy-distributed noise
    bool save//<! Determines whether particles will be written to file
)
{
    vector<int> data = df.C;
    vector<int> yradj = df.highyears;

    double efficacy = 0.85;
    int sw = mcv1_eff[countryISO3];
    if(sw == 1)
    {
        efficacy = 0.93;
    }


#ifdef SUSCEPTIBLE_METHOD_1
    // age.dist <- dgeom(c(0:(n.ages-1)),.25)
    // age.dist[n.ages] = pgeom(n.ages, .25,lower.tail = FALSE)
    ublas::vector<double> agedist(nages);
    ublas::vector<double> S0N(npart);
    double const lambda = 0.25;//1.0;
    auto d = boost::math::geometric_distribution<double> {lambda};

    for(int i=0; i<nages; ++i)
    {
        agedist(i) = pdf(d, i);
    }

    agedist(nagesDecremented) = (1.0 - cdf(d,nagesDecremented));
    agedist[99] = 2.405402e-13; //boost::geometric_distribution and pgeom differ, we use pgeom value

    //// write agedist to file
    // ofstream agedistfile;
    // agedistfile.open("agedist.csv");
    // for(int i=0; i<nages; ++i)
    // {
    //    agedistfile << agedist(i) << endl;
    //    
    //    //cout << agedist(i) << endl;

    //    //assert(agedist(i) > 0);
    //    //assert(agedist(i) < 1);
    //    //assert(!isnan(agedist(i)));
    //    //assert(!isinf(agedist(i)));
    // }
    // agedistfile.close();

    for(int i=0; i<npart; ++i)
    {
        S0N(i) = S0[i] * df.N[0] ;
    }

    // states[1, 1,,] <- matrix((S0*N[1])%x%age.dist, nrow = n.ages, byrow = F)	# initialize with 3% susceptible
    ublas::matrix<double> kron(nages,npart);
    kron = outer_prod(agedist,S0N);

    for(int p=0; p<nages; p++)
    {
        for(int q=0; q<npart; q++)
        {
            states[0][0][0][p][q] = kron(p,q);
        }
    }
#endif

#ifdef SUSCEPTIBLE_METHOD_2


    // r <- (S0*N[1]) / (B[1] + S0*N[1])
    double r[npart];

    double N = df.N[0];
    double B = df.B[0];

                                                                                                                                                                                                                                                                                                                                                                                   
 // # calculate the age distribution for the initial susceptible population -- B[1] is birth cohort in first year; N[1] is popualtion size in first year; mcv1[1] is vaccination coverage; .83 is assumed vaccine efficacy                                                                                                                                                                      
 // B0 <- B[1] * (1-mcv1[1]*.83)     # This is the assumed firt unvaccinated birth cohort. For coutnries with no vaccination, this is B[1]; for countries with high coverage in 1980, this is the fraction unimmunized by vaccination                                                                                                                                                           
 // r <- (S0*N[1]) / (B0 + S0*N[1])    # r is the decay rate of a geometric series that guarantees that the sum of all susceptibles is equal to S0*N                                                                                                                                                                                                                                            
 // # generate a unique decay rate for each of the random draws of S0; should be of length particle                                                                                                                                                                                                                                                                                             
                                                                                         


    double B0 = B * (1.0 - (df.mcv1[0] * 0.83)); 

    for(int p=0; p<npart; p++)
    {
        //r[p] = (S0[p] * N) / (B + S0[p] * N);
        r[p] = (S0[p] * N) / (B0 + S0[p] * N);
    }
    ////// write r to file
    //std::stringstream rfss;
    //rfss << "r_" << df.countryModel << "_" << df.shatimestamp << ".csv";
    //ofstream rfile;
    //rfile.open(rfss.str());
    //for(int p=0; p<npart; p++)
    //{
    //    rfile << r[p] << endl;
    //}
    //rfile.close();


    //// sum the values of r
    //double rsum = 0.0;
    //for(int p=0; p<npart; p++)
    //{
    //    rsum += r[p];
    //}
    ////cout << "rsum: " << rsum << endl;


    for(int p=0; p<nages; p++)
    {
        for(int q=0; q<npart; q++)
        {
        // states[1, 1,,ii] <- B[1] * r[ii]^(1:n.ages)
            states[0][0][0][p][q] = B0 * (pow(r[q],(p+1)));
        }
    }
#endif

    stringstream s_particles;
    //stringstream s_particles_full;
    for(int particle=0; particle<npart; particle++)
    {
        double sparticle = 0;
        for(int age=0; age<nages; age++)
        {
            sparticle += states[0][0][0][age][particle]; // We are summing all age classes here
        //// write the value to the stringstream using full precision to avoid rounding errors
        ////  s_particles_full << std::fixed << std::setprecision(16) << states[0][0][0][age][particle] << " ";
           // s_particles_full << states[0][0][0][age][particle] << " ";
        }
        s_particles << sparticle << " ";
        //s_particles_full << "\n";
    }
    s_particles << "\n";

    //// write s_particles to file
    //ofstream s_particles_file;
    //std::stringstream fnss;
    //fnss << df.countryModel << "_sap_" << df.shatimestamp << ".txt";
    //s_particles_file.open(fnss.str());
    //s_particles_file << s_particles.str();
    //s_particles_file.close();


    //std::stringstream sfnss;
    //ofstream s_particles_file;
    //sfnss << df.countryModel << "_sapFULL_" << df.shatimestamp << ".txt";
    //s_particles_file.open(sfnss.str());
    //s_particles_file << s_particles_full.str();
    //s_particles_file.close();
    ////exit(0);

    MatrixXd msit(calendaryrs, npart);
    MatrixXd msst(calendaryrs, npart);

    MatrixXd observationwt(calendaryrs, npart);
    MatrixXd observationlik(calendaryrs, npart);
    for(int p=0; p<npart; p++)
    {
        observationwt(0,p) = 0.0;
        observationlik(0,p) = 0.0;
    }

    int ft = 1;  // time index
    int bft = 0; // time index to previous time step

    stringstream i_particles;

#ifdef PRINTAGES
    stringstream i_age_particles;
    stringstream s_age_particles;
#endif
#ifdef PRINTAGESSPLIT
    stringstream i_age_split_particles;
    stringstream s_age_split_particles;
#endif
    if(save) {
        // For R comparison:
        i_particles << "NA";
        for(int h=1; h<npart; h++)
        {
            i_particles << " NA";
        }
        i_particles << "\n";

    }// end save

    for(int t=1; t<calendaryrs; t++)   // TIME STARTS HERE AT 1; TIMESTEP 0 ALREADY HAPPENED
    {

#ifndef DROPN
        for(int a=0; a<nages; a++)
        {
            for(int p=0; p<npart; p++)
            {
                states[0][ft][2][a][p] = df.N[t];
            }
        }
#endif

        double theta = theta_value; // 0.05;
        VectorXd xpi = simstep(t,ft,bft,states,params,sigma,theta,df,efficacy);

#ifndef DROPN
        // states[i, 3,,] <- N[i]
        for(int p=0; p<nages; p++)
        {
            for(int q=0; q<npart; q++)
            {
                states[0][ft][2][p][q] = df.N[t];
            }
        }
#endif


        MatrixXd ItStates(nages,npart);
        ItStates.setZero();
        for(unsigned i=0; i<npart; i++)
        {
            for(unsigned j=0; j<nages; j++)
            {
                ItStates(j,i) = states[0][ft][1][j][i]; // only updates current T
            }
        }

        VectorXd vobservationwt = obswtIt(t, ItStates, params, data, yradj[t], cauchy_weight);
        observationwt.row(t) = vobservationwt;

        // Find maximum observation weight
        double maxob = 0.0;
        for(int g=0; g<npart; g++)
        {
            if(g==0)
            {
                maxob = vobservationwt[g];
            }
            if(vobservationwt[g]>maxob)
            {
                maxob = vobservationwt[g];
            }
        }

        //smp.wt <- (exp(observation.wt[i,] - max(observation.wt[i,], na.rm = TRUE))) / sum(exp(observation.wt[i,] - max(observation.wt[i,], na.rm = TRUE)), na.rm= TRUE)
        VectorXd psmpwt(npart);
        VectorXd smpwt(npart);

        for(int q=0; q<npart; q++)
        {
            psmpwt(q) = exp(vobservationwt(q) - maxob);
        }

        double sum = psmpwt.sum();
        smpwt = psmpwt/sum;

        VectorXd newind(npart);
        VectorXd newindcount(npart);
        newindcount.setZero();

        vector<ptuple> ptable;
        ptable.reserve(npart);

        for(int r=0; r<npart; r++)
        {
            ptuple p;
            p.index = r;
            p.pf = smpwt[r];
            p.cf = 0.0;
            p.count = 0;
            ptable.push_back(p);
        }
        // sort:
        std::sort(ptable.begin(), ptable.end(),
                  [](const ptuple & a, const ptuple & b) -> bool
        {
            return a.pf < b.pf;
        });

        double runprob = 0;
        for(int r=0; r<npart; r++)
        {
            ptable[r].cf = runprob + ptable[r].pf;
            runprob += ptable[r].pf;
        }

        std::uniform_real_distribution<> dis(0.0, 1.0);
        for(int q=0; q<npart; q++)
        {
            double sprob = dis(generator);

            int fndr = 0;
            while(ptable[fndr].cf<sprob)
            {
                fndr++;
            }

            int ap = ptable[fndr].index;

            newind[q] = ap;
            newindcount[ap] += 1;
        }

#ifdef DEBUGWEIGHTS
        for(int r=0; r<npart; r++)
        {
            cout << "PTABLE " << r << " " << ptable[r].index << " " << ptable[r].pf << " " << ptable[r].cf << endl;
        }
        vector<double> pcase_counts(npart);
        for(int y=0; y<npart; y++)
        {
            pcase_counts.push_back(0.0);
        }
        for(int particle=0; particle<npart; particle++)
        {
            double fulli = 0;
            for(int age=0; age<nages; age++)
            {
                fulli += states[0][ft][1][age][particle];
            }
            pcase_counts[particle] = fulli;
        }
        for(int y=0; y<npart; y++)
        {
            cout << y << " CC: " << pcase_counts[y] << endl;
        }
#endif

#ifdef DEBUG_PARTICLES
        // find p for each particle:
        vector<double> ppartable;
        ppartable.reserve(npart);
        for(int r=0; r<npart; r++)
        {
            ppartable[ptable[r].index] = ptable[r].pf;
        }

        vector<double> pcartable;
        pcartable.reserve(npart);
        for(int r=0; r<npart; r++)
        {
            pcartable[ptable[r].index] = ptable[r].cf;
        }

        stringstream dbg_particles;
        dbg_particles << "S N SN pinf I newp fprob smp pf cf" << endl;
        for(int particle=0; particle<npart; particle++)
        {
            double sfullpar = 0.0;
            double ifullpar = 0.0;
            double singlenpar = 0.0;
            for(int par=0; par<npars; par++)
            {
                double fullpar = 0.0;
                for(int age=0; age<nages; age++)
                {
                    fullpar += states[0][ft][par][age][particle];
#ifndef DROPN 
                    if (par == 2)
                    {
                        singlenpar = states[0][ft][par][age][particle];
                        // cout << states[0][ft][par][age][particle] << endl;
                    }
#endif
                }
                if(par==0)
                {
                    sfullpar = fullpar;
                }
                else if (par == 1)
                {
                    ifullpar = fullpar;
                }
            }

            dbg_particles << sfullpar << " " << singlenpar << " " << (sfullpar/singlenpar) << " " << xpi[particle] << " " << ifullpar << " " << newindcount[particle] << " " ;
            dbg_particles << ppartable[particle] << " " <<  smpwt[particle] << " " << ptable[particle].pf << " " << ptable[particle].cf ;
            dbg_particles << endl;
        }

        std::stringstream a_fn_ss_dbg;
        a_fn_ss_dbg << df.countryModel << "_dbg_particles_fx_save" << t << df.shatimestamp << ".txt";
        std::string a_fn_strd = a_fn_ss_dbg.str();
        std::ofstream a_ofstd(a_fn_strd.c_str());
        if (a_ofstd.is_open())
        {
            a_ofstd << dbg_particles.str();
        }
        else
        {
            std::cerr << "Could not open " << a_fn_strd << std::endl;
        }
#endif


        double newstates[npars][nages][npart];// sample particles using probability of sample weight
        int app = 0;
        for(int par=0; par<npars; par++)
        {
            for(int age=0; age<nages; age++)
            {
                for(int q=0; q<npart; q++)
                {
                    app = newind[q];
                    newstates[par][age][q] = states[0][ft][par][age][app];// REPLACEMENT

                    if(isinf(states[0][ft][par][age][app]))
                    {
                        cout << __LINE__ << ft << " " << par << " " << age << " " << app << endl;
                        assert(0);
                    }
                    if(isnan(states[0][ft][par][age][app]))
                    {
                        cout << __LINE__ << " " << ft << " " << par << " " << age << " " << app << endl;
                        assert(0);
                    }
                }
            }
        }

        if(df.C[t] != -1 )  // If Cases are not NA, actually update according to weights; otherwise simply propagate all particles
        {
            for(int par=0; par<npars; par++)
            {
                for(int age=0; age<nages; age++)
                {
                    for(int particle=0; particle<npart; particle++)
                    {
                        states[0][ft][par][age][particle] =  newstates[par][age][particle];
                    }
                }
            }
        }


        //Record particle compartments for output:
        for(int particle=0; particle<npart; particle++)
        {
            unsigned ciparticle = 0;
            unsigned sparticle = 0;
#ifdef PRINTAGESSPLIT
            double agesum_i = 0.0;
            double agesum_s = 0.0;
#endif
            for(int age=0; age<nages; age++)
            {
                ciparticle += states[0][ft][1][age][particle];
                sparticle += states[0][ft][0][age][particle]; // We are summing all age classes here
#ifdef PRINTAGES

                if(save) {
                    i_age_particles <<  states[0][ft][1][age][particle] << ' ';
                    s_age_particles <<  states[0][ft][0][age][particle] << ' ';
                }
#endif

#ifdef PRINTAGESSPLIT
                agesum_i += states[0][ft][1][age][particle];
                agesum_s += states[0][ft][0][age][particle];
                if(age == 4)
                {
                    i_age_split_particles << agesum_i << ' ';
                    s_age_split_particles << agesum_s << ' ';
                    agesum_i = 0.0;
                    agesum_s = 0.0;
                }
                if(age == 99)
                {
                    i_age_split_particles << agesum_i << ' ';
                    s_age_split_particles << agesum_s << ' ';
                    agesum_i = 0.0;
                    agesum_s = 0.0;
                }
#endif
            }

            if(save) {
                i_particles << ciparticle << " ";
                s_particles << sparticle << " ";
#ifdef PRINTAGES
#ifdef PRINTAGES_BYYEAR
                i_age_particles << "\n";
                s_age_particles << "\n";
#endif
#endif

//                i_particles << "\n";
//                s_particles << "\n";


#ifdef PRINTAGESSPLIT
                i_age_split_particles << "\n";
                s_age_split_particles << "\n";
#endif
            }
        }

#ifdef PRINTAGES
#ifndef PRINTAGES_BYYEAR
        if(save) {
            i_age_particles << "\n";
            s_age_particles << "\n";
                i_particles << "\n";
                s_particles << "\n";
        }
#endif
#endif



#ifdef PRINTAGES
#ifdef PRINTAGES_BYYEAR
        if(save) {
                i_particles << "\n";
                s_particles << "\n";
        }
#endif
#endif


        // Step time indicies here:
        if(ft==1)
        {
            ft=2;
            bft=1;
        }
        else
        {
            ft = 1;
            bft= 2;
        }


    } // end of calendaryrs

    if(save) {
        std::stringstream a_fn_ss;
        a_fn_ss << df.countryModel << "_i_particles_fxsave" << df.shatimestamp << ".txt";
        std::string a_fn_str = a_fn_ss.str();
        std::ofstream a_ofst(a_fn_str.c_str());
        //df.i_particles_sha = get_sha1(i_particles.str());
        if (a_ofst.is_open())
        {
            a_ofst << i_particles.str();
        }
        else
        {
            std::cerr << "Could not open " << a_fn_str << std::endl;
        }

        std::stringstream a_fn_sss;
        a_fn_sss << df.countryModel << "_s_particles_fxsave" << df.shatimestamp << ".txt";
        std::string a_fn_strs = a_fn_sss.str();
        std::ofstream a_ofsts(a_fn_strs.c_str());
        if (a_ofsts.is_open())
        {
            a_ofsts << s_particles.str();
        }
        else
        {
            std::cerr << "Could not open " << a_fn_strs << std::endl;
        }

        df.i_particles_sha = get_sha1(i_particles.str());
        df.s_particles_sha = get_sha1(s_particles.str());

#ifdef PRINTAGES
        std::stringstream iap_fn_ss;
        iap_fn_ss << df.countryModel << "_iap_particles_fxsave" << df.shatimestamp << ".txt";
        std::string iap_fn_str = iap_fn_ss.str();
        std::ofstream iap_ofst(iap_fn_str.c_str());
        if (iap_ofst.is_open())
        {
            iap_ofst << i_age_particles.str();
        }
        else
        {
            std::cerr << "Could not open " << iap_fn_str << std::endl;
        }

        //if(df.action == 5){
        std::stringstream sap_fn_ss;
        sap_fn_ss << df.countryModel << "_sap_particles_fxsave" << df.shatimestamp << ".txt";
        std::string sap_fn_str = sap_fn_ss.str();
        std::ofstream sap_ofst(sap_fn_str.c_str());
        if (sap_ofst.is_open())
        {
            sap_ofst << s_age_particles.str();
        }
        else
        {
            std::cerr << "Could not open " << sap_fn_str << std::endl;
        }
        
        df.iap_particles_sha = get_sha1(i_age_particles.str());
        df.sap_particles_sha = get_sha1(s_age_particles.str());
        //}
#endif


#ifdef PRINTAGESSPLIT
        std::stringstream spiap_fn_ss;
        spiap_fn_ss << df.countryModel << "_iap_split_particles_fxsave" << df.shatimestamp << ".txt";
        std::string spiap_fn_str = spiap_fn_ss.str();
        std::ofstream spiap_ofst(spiap_fn_str.c_str());
        if (spiap_ofst.is_open())
        {
            spiap_ofst << i_age_split_particles.str();
        }
        else
        {
            std::cerr << "Could not open " << spiap_fn_str << std::endl;
        }

        std::stringstream spsap_fn_ss;
        spsap_fn_ss << df.countryModel << "_sap_split_particles_fxsave" << df.shatimestamp << ".txt";
        std::string spsap_fn_str = spsap_fn_ss.str();
        std::ofstream spsap_ofst(spsap_fn_str.c_str());
        if (spsap_ofst.is_open())
        {
            spsap_ofst << s_age_split_particles.str();
        }
        else
        {
            std::cerr << "Could not open " << spsap_fn_str << std::endl;
        }
#endif

    }// end save

    VectorXd forsum(calendaryrs);
    VectorXd forsumll(calendaryrs);
    for(int t=0; t<calendaryrs; t++)
    {
        double rsum = 0.0;
        double rsumlik = 0.0;
        std::vector<double> nps(npart, 0.0);
        std::vector<double> npsll(npart, 0.0);
        for(int particle=0; particle<npart; particle++)
        {
            nps[particle] = observationwt(t,particle);
            npsll[particle] = observationlik(t,particle);
        }
        std::sort (nps.begin(),nps.end());
        rsum = std::accumulate(nps.begin(), nps.end(), 0.0);
        forsum(t) = rsum/static_cast<double>(npart);

        std::sort (npsll.begin(),npsll.end());
        rsumlik = std::accumulate(npsll.begin(), npsll.end(), 0.0);
        forsumll(t) = rsumlik/static_cast<double>(npart);
    }

// construct and apply penalty:
    double depletedyearsum = 0.0;
    for(unsigned cy=0; cy<calendaryrs; ++cy)
    {
        if (forsum(cy) < -248.5)
        {
            depletedyearsum += 1.0;
        }
    }
    double penalty = (depletedyearsum*depletedyearsum);
    double nLLt = -(forsum.sum());
    double nLL = -(forsum.sum()) + penalty;
    double lik = forsumll.sum();
// End of penalty application

    VectorXd states00sum(npart);
    states00sum.setZero();

    for(int a=0; a<nages; a++)
    {
        for(int p=0; p<npart; p++)
        {
            states00sum(p) += states[0][0][0][a][p];
        }
    }

    double S0pct = states00sum.mean() / df.N[0];


    if(save) {
        std::cout <<"["<<__FILE__<<"|"<<__FUNCTION__<<"|"<<__LINE__<<"]"<< " S0pct: " << S0pct << " " << std::endl;
    }

    fx_return_struct fxr;

    fxr.nLL=nLL;
    fxr.lik=lik;
    fxr.nLLt=nLLt;
    fxr.depletedyearsum=depletedyearsum;
    fxr.penalty=penalty;
    // Iest_struct fields follow:
    fxr.S0pct=S0pct;
    fxr.penalty=penalty;
    fxr.paramsBest = params;



    //return Iest;
    return fxr;
}










//! Loop over the parameter values for grid search, storing the likelihoods recieved from the particle filter
// * @param[out] A vector containing parameters, negative log likelihoods, penalties for particle collapse
vector< array<double, gridsize> > gridSearch(
    double (*states)[nyrs][npars][nages][npart], //<!
    const MatrixXd& params, //<! All parameter sets to be sampled in this grid search
    double *S0, //<! Pointer to an array of inital S0 values
    demographicDataFrame df, //<! Demographic DataFrame
    const string& countryISO3, //<! ISO3 code
    double sigma, //<! Noise term
    const double cauchy_weight //<! Scale Cauchy distribution contribution to particle observation weights
)
{
    // Create required arrays for function output:
    array<double,gridsize> nLLo;
    array<double,gridsize> lik;
    array<double,gridsize> nllt;
    array<double,gridsize> depletedyearsum;
    array<double,gridsize> penalty;

    // Loop over states and initialize:
    for(int j=0; j<gridsize; j++)
    {
        int toplevel = 1;
        for (int tl = 0; tl<toplevel; tl++)
        {
            for(int y=0; y<nyrs; y++)
            {
                for(int p=0; p<npars; p++)
                {
                    for(int a=0; a<nages; a++)
                    {
                        for(int x=0; x<npart; x++)
                        {
                            states[tl][y][p][a][x] = 0.0;
                        }
                    }
                }
            }
        }

        RowVectorXd io = params.row(j); // feed specific parameter set here


        bool save = false;
        fx_return_struct fxcr = fx( states,
                                    io,
                                    S0,
                                    df,
                                    countryISO3,
                                    sigma,
                                    cauchy_weight,
                                    save
                                  );


        nLLo[j] = fxcr.nLL;
        lik[j] = fxcr.lik;
        nllt[j] = fxcr.nLLt;
        depletedyearsum[j] = fxcr.depletedyearsum;
        penalty[j] = fxcr.penalty;

        std::cout << j << " " << io << " " << fxcr.nLL << std::endl;

    }

    vector< array<double, gridsize> > fxr;
    fxr.push_back(nLLo);
    fxr.push_back(lik);
    fxr.push_back(nllt);
    fxr.push_back(depletedyearsum);
    fxr.push_back(penalty);

    return(fxr);

}



//! Generates mcv2 effective coverage under assumption that second dose is given only to those children who received the first dose
double mcv2_dep(
    const double mcv1, //!< Coverage with first dose
    const double mcv2, //!< Coverage with second dose
    const double ve//!< vaccine efficacy of first dose (countryISO3 specific based on age of 1st dose)
)
{
    return ((1.0-ve)*mcv2)/(1-mcv1*ve);
}


//! Generates SIA effective coverage under assumption that campaigns preferentially vaccinate children who received the first routine dose
double sia_dep(const double mcv1,//!< Coverage with first dose
               const double sia,//!< Coverage of SIA
               const double ve//!< vaccine efficacy of first dose (countryISO3 specific based on age of 1st dose)
              )
{
    double sia_new = 0.0;
    if(sia>0.0)
    {
        double p = (mcv1 * (1.0 - ve));// fraction not immunized by routine
        if(sia<mcv1)
        {
            sia_new = sia/mcv1 * p + (1.0-p) * 0;// SIA only reaches those not immunized by routine
        }
        else
        {
            sia_new = p * 1.0 + (1.0-p) * (sia-mcv1)/(1.0-mcv1);// SIA reaches into those missed by routine
        }

    }
    return sia_new;
}

//! Discount MCV2 efficacy for models which require it.
void scale_mcv2(demographicDataFrame &frame)
{
    double efficacy = 0.85;
    int sw = mcv1_eff[frame.countryISO3];
    if(sw == 1)
    {
        efficacy = 0.93;
    }
    int yrs = frame.N.size();
    for(int t=0; t<yrs; t++)
    {
        double cmcv2 = frame.mcv2[t];
        frame.mcv2[t] = mcv2_dep(frame.mcv1[t],cmcv2,efficacy);
    }
}

//! Discount SIA efficacy for models which require it.
void scale_sia(demographicDataFrame &frame)
{
    double efficacy = 0.85;
    int sw = mcv1_eff[frame.countryISO3];
    if(sw == 1)
    {
        efficacy = 0.93;
    }
    int yrs = frame.N.size();
    for(int t=0; t<yrs; t++)
    {
        cout << t << " " << yrs << endl;
        double sia_orig = frame.SIA[t];
        if(sia_orig > 0.0)
        {
            cout << "sia_orig: " << sia_orig << endl;
            frame.SIA[t] = sia_dep(frame.mcv1[t],sia_orig,efficacy);
            for(int a=0; a<nages; a++)
            {
                sia_orig = frame.SIA_matrix[a][t+1];
                if(sia_orig > 0.0)
                {
                    cout << "---------------" << endl;
                    cout << "AT: " << t << " " << a << endl;
                    cout << "---------------" << endl;
                }
                frame.SIA_matrix[a][t+1] = sia_dep(frame.mcv1[t],sia_orig,efficacy);
            }
        }
    }
}


//! Create and search initial parameter grid, then create and search constrained parameter grid.
//! @param[out] An Iest_struct containing best-fit parameters and nLL
Iest_struct initiateGridSearch(
    const cli_inputs& cli,//!< Command-line input struct
    //const demographicDataFrame& df,//!< Demographic DataFrame
    demographicDataFrame& df,//!< Demographic DataFrame
    const string& countryISO3,//!< countryISO3 ISO3 code. "AFG"
    const string& countryModel//!< countryISO3 ISO3 code + model.  "ETHnormal"
)
{
    cout << "allocating states" << endl;
    auto states = new double[1][nyrs][npars][nages][npart]();
    cout << "allocated states" << endl;

    int toplevel = 1;
    for (int tl = 0; tl<toplevel; tl++)
    {
        for(int y=0; y<nyrs; y++)
        {
            for(int p=0; p<npars; p++)
            {
                for(int a=0; a<nages; a++)
                {
                    for(int x=0; x<npart; x++)
                    {
                        states[tl][y][p][a][x] = 0.0;
                    }
                }
            }
        }
    }


    // #initial values for the number of susceptibles (expressed as proportion of total population)
    std::uniform_real_distribution<double> distributionS0(0.001,0.2);
    double S0[npart]= {};
    for (int i=0; i<npart; ++i)
    {
        S0[i] =  distributionS0(generator);
    }

    //std::ifstream file("S01k.txt");
    //std::string line;
    //int index = 0;

    //if (file.is_open()) {
    //    while (getline(file, line) && index < 1000) {
    //        try {
    //            S0[index] = std::stod(line);
    //            index++;
    //        } catch (const std::invalid_argument& e) {
    //            std::cerr << "Invalid number found in the file: " << line << std::endl;
    //        }
    //    }
    //    file.close();
    //} else {
    //    std::cerr << "Unable to open file" << std::endl;
    //}

    //// Write S0 to file
    //ofstream S0file;
    //std::stringstream fnss;
    //fnss << countryModel << "_S0_" << df.shatimestamp << ".txt";
    //S0file.open(fnss.str());
    //for (int i=0; i<npart; ++i)
    //{
    //    S0file << S0[i] << "\n";
    //}
    //S0file.close();

    // Generate beta0 values for grid cells:
    std::uniform_real_distribution<double> distribution_beta0(-10.0,-4.0);
    double b0inits[gridsize]= {};
    for (int i=0; i<gridsize; ++i)
    {
        b0inits[i] =  distribution_beta0(generator);
    }

    std::uniform_real_distribution<double> distribution_pobs(0.001,0.1);
    std::uniform_real_distribution<double> distributione(0.001,0.2);
    MatrixXd params_inits4par(gridsize,4);
    for (int i=0; i<gridsize; ++i)
    {
        double b0100 = -100.0 * b0inits[i];
        double b1upper = 200.0;

        if(b0100<b1upper)
        {
            b1upper = b0100;
        }
        std::uniform_real_distribution<double> distribution_beta1((4*(-b0inits[i])),b1upper);
        params_inits4par(i,0)=b0inits[i];
        // Generate beta1 values for this grid cell:
        params_inits4par(i,1)=distribution_beta1(generator);
        // Generate probability of reporting value for this grid cell:
        params_inits4par(i,2)=distribution_pobs(generator);
        // Generate probability of "high" reporting values for this grid cell:
        params_inits4par(i,3)= params_inits4par(i,2) + distributione(generator);
        if(cli.fixed_prhigh == true)
        {
            assert(cli.fixed_prhigh == false);
            params_inits4par(i,3)=cli.prhigh;
        }
    }


    vector< array<double,gridsize> > run1;
    run1 = gridSearch(states,
                      params_inits4par,
                      S0,
                      df,
                      countryISO3,
                      cli.sigma,
                      cli.cauchy_weight);


    // ## # 2. SECOND RUN

    // # get the range of the best runs from run 1
    // rg <- apply(run1$params[order(run1$nLL)[1:100],],2,range, na.rm = TRUE)

    vector<tuple<int, double>> v;
    for (int i=0; i<gridsize; ++i)
    {
        v.push_back(std::make_tuple(i,run1[0][i]));
    }

    sort(v.begin(),v.end(),
         [](const tuple<int,double>& a,
            const tuple<int,double>& b) -> bool
    {
        return std::get<1>(a) < std::get<1>(b);
    });

    unsigned bestparamsize = 100;// How many "best" param sets to bound beta search?
    MatrixXd bestparams(bestparamsize,4);// Number of "best" sets by number of parameters

    for (unsigned i=0; i<bestparamsize; ++i)
    {
        bestparams.row(i) = params_inits4par.row(get<0>(v[i]));
        cout << "BESTPARAMS " << i << " " << bestparams.row(i) << endl;
    }

    MatrixXd rg(2,4);
    rg.row(0) = bestparams.colwise().minCoeff();
    rg.row(1) = bestparams.colwise().maxCoeff();

    cout << "rg.row(0): " << rg.row(0) << endl;
    cout << "rg.row(1): " << rg.row(1) << endl;

    //  # limit grid to the best runs from run 1
    // b0inits <- runif(gridsize,rg[1,1],rg[2,1])
    // #include "b0initsRUN2.h"
    //     b0inits = b0inits2;

    // params.initsRED <- matrix(c(b0inits,
    //                             runif(gridsize, 4*(-b0inits),rg[2,2]),
    //                             runif(gridsize, rg[1,3],rg[2,3])),
    // nr = gridsize,nc = 3)


    std::uniform_real_distribution<double> distributionb2(rg(0,0),rg(1,0));
    for (int i=0; i<gridsize; ++i)
    {
        b0inits[i] = distributionb2(generator);
    }

    std::uniform_real_distribution<double> distributionc2(rg(0,2),rg(1,2));
    std::uniform_real_distribution<double> distribution_prhigh(rg(0,3),rg(1,3));

    MatrixXd params_initsRED(gridsize,4);
    for (int i=0; i<gridsize; ++i)
    {
        params_initsRED(i,0)=b0inits[i];
        std::uniform_real_distribution<double> distributiond2((4*(-b0inits[i])),rg(1,1));
        params_initsRED(i,1)=distributiond2(generator);
        params_initsRED(i,2)=distributionc2(generator);
        double prhighsc = distribution_prhigh(generator);
        while(prhighsc < params_initsRED(i,2))
        {
            prhighsc = distribution_prhigh(generator);
        }
        params_initsRED(i,3)= prhighsc;

        if(cli.fixed_prhigh == true)
        {
            params_initsRED(i,3)=cli.prhigh;
        }
    }


    std::cout <<"["<<__FILE__<<"|"<<__FUNCTION__<<"|"<<__LINE__<<"]"<<  "SECOND RUN " << std::endl;

    vector< array<double,gridsize> > run2;
    run2 = gridSearch(
               states,
               params_initsRED,
               S0,
               df,
               countryISO3,
               cli.sigma,
               cli.cauchy_weight);

    // # 3. GENERATE STATE CI
    // best <- run2$params[which.min(run2$nLL),]

    int bestIndex = 0;
    double minnLL = run2[0][0];
    for(unsigned i=0; i<gridsize; i++)
    {
        if(run2[0][i]<minnLL)
        {
            bestIndex = i;
            minnLL = run2[0][i];
        }
    }

    std::cout <<"["<<__FILE__<<"|"<<__FUNCTION__<<"|"<<__LINE__<<"]"<< " bestIndex: " << bestIndex << " " << std::endl;
    std::cout <<"["<<__FILE__<<"|"<<__FUNCTION__<<"|"<<__LINE__<<"]"<< " minnLL: " << minnLL << " " << std::endl;

    VectorXd paramsBest;
    paramsBest = params_initsRED.row(bestIndex);

    std::cout <<"["<<__FILE__<<"|"<<__FUNCTION__<<"|"<<__LINE__<<"]"<< " paramsBest: " << paramsBest << " " << std::endl;

    // Iest<-fx(best, S0, df = df,n.ages, age.dist, SIA, particle = npart, high.years = df$high.years, directory=save.here, save = TRUE,country=countryISO3, sigma=s)



    bool save = true;
    fx_return_struct fxr = fx(states,
                              paramsBest,
                              S0,
                              df,
                              countryISO3,
                              cli.sigma,
                              cli.cauchy_weight,
                              save);

    Iest_struct bestIest;

    bestIest.nLL = fxr.nLL;
    bestIest.S0pct = fxr.S0pct;
    bestIest.paramsBest = fxr.paramsBest;


    std::stringstream run1c;
    run1c << "b0 b1 pr prhigh prepenalty_nLL depletedyears penalty nLL grid\n";
    for(unsigned i=0; i<gridsize; i++)
    {
        run1c << params_inits4par(i,0) << " "<< params_inits4par(i,1) << " "<< params_inits4par(i,2) << " " << params_inits4par(i,3) << " " << run1[2][i] << " " << run1[3][i] << " " << run1[4][i] << " " << run1[0][i] << " 1" <<   "\n";
    }

    for(unsigned i=0; i<gridsize; i++)
    {
        run1c << params_initsRED(i,0) << " "<< params_initsRED(i,1) << " "<< params_initsRED(i,2) << " " <<  params_initsRED(i,3) << " " << run2[2][i] << " " << run2[3][i] << " " << run2[4][i] << " " << run2[0][i] << " 2" <<  "\n";
    }

    std::stringstream fnss2;
    fnss2 << countryModel << "_runs_showing_penalty" << df.shatimestamp << ".txt";
    std::ofstream filerun2(fnss2.str());
    if (filerun2.is_open())
    {
        filerun2 <<  run1c.str();
    }
    filerun2.close();
    df.runs_showing_penalty_sha = get_sha1(run1c.str());

    delete[] states;
    return(bestIest);
}

//! Run a single parameter set, without setting up a full grid search.
void runSingleParameterSet(
    const cli_inputs& cli,
    //const demographicDataFrame& df,
    demographicDataFrame& df,
    const string& countryISO3,
    const string& countryModel
)
{
    auto states = new double[1][nyrs][npars][nages][npart]();


    // initial values for the number of susceptibles
    std::uniform_real_distribution<double> distributionS0(0.001,0.2);// S0 <- runif(npart, .001, .2)
    double S0[npart]= {};
    for (int i=0; i<npart; ++i)
    {
        S0[i] =  distributionS0(generator);
    }

   // std::ifstream file("S01k.txt");
   // std::string line;
   // int index = 0;

   // if (file.is_open()) {
   //     while (getline(file, line) && index < 1000) {
   //         try {
   //             S0[index] = std::stod(line);
   //             index++;
   //         } catch (const std::invalid_argument& e) {
   //             std::cerr << "Invalid number found in the file: " << line << std::endl;
   //         }
   //     }
   //     file.close();
   // } else {
   //     std::cerr << "Unable to open file" << std::endl;
   // }


   // // Write S0 to file
   // ofstream S0file;
   // std::stringstream fnss;
   // fnss << countryModel << "_S0_" << df.shatimestamp << ".txt";
   // S0file.open(fnss.str());
   // for (int i=0; i<npart; ++i)
   // {
   //     S0file << S0[i] << "\n";
   // }
   // S0file.close();

    Iest_struct iestBest = readIest(cli.iestfile);

    VectorXd paramsBest;
    paramsBest = iestBest.paramsBest;

    std::cout <<"["<<__FILE__<<"|"<<__FUNCTION__<<"|"<<__LINE__<<"]"<< " paramsBest: " << paramsBest << " " << std::endl;


    bool save = true;

    fx_return_struct fxr = fx(states,
                              paramsBest,
                              S0,
                              df,
                              countryISO3,
                              cli.sigma,
                              cli.cauchy_weight,
                              save);

    Iest_struct bestIest;

    bestIest.nLL = fxr.nLL;
    bestIest.S0pct = fxr.S0pct;
    bestIest.paramsBest = fxr.paramsBest;

    cout << "nLL: " << bestIest.nLL << endl;
    cout << "S0: "  << bestIest.S0pct << endl;
    cout <<  "params " << bestIest.paramsBest << endl;
    writeIest(bestIest,cli.country_model,1,df.shatimestamp);


    delete[] states;
}




int main(int argc, char** argv)
{


#ifdef STATICSEED
    // Logging 10 random numbers; only needed for checking state of
    // RNG when comparing runs with static RNG
    std::uniform_real_distribution<double> distributionTest(0.001,0.1);
    for (int i=0; i<10; ++i)
    {
        cout <<  distributionTest(generator) << endl;
    }
#endif

    cli_inputs cli = cli_inputs();

    // Read, convert, store command-line arguments:
    cli.action = atoi(argv[1]);
    cli.country_model = string(argv[2]);
    cli.demofile = string(argv[3]);
    cli.siafile = string(argv[4]);
    cli.modelnumber = atoi(argv[5]);
    cli.sigma = atof(argv[6]);
    cli.cauchy_weight = atof(argv[7]);

    cli.timestamp = get_timestamp();

    std::vector<std::string> args;
    std::copy(argv + 1, argv + argc, std::back_inserter(args));
    std::stringstream arg_ss;
    for(unsigned i=0;i<args.size();i++){
        arg_ss << args[i] << " ";
    }

cout << arg_ss.str() << endl;

    cli.fullcl = arg_ss.str();

    if(cli.action == 5)
    {
        cli.iestfile = string(argv[8]);
    }

    cout << "action: " << cli.action << endl;
    cout << "demofile: " << cli.demofile << endl;
    cout << "siafile: " << cli.siafile << endl;
    cout << "modelnumber: " << cli.modelnumber << endl;
    cout << "sigma: " << cli.sigma << endl;
    cout << "cauchy_weight: " << cli.cauchy_weight << endl;

    demographicDataFrame df=demographicDataFrame();
    df.timestamp = cli.timestamp;
    df.modelnumber = cli.modelnumber;

    for(int age=0; age<nages; age++)
    {
        for(int year=0; year<SIAcalendaryrs; year++)
        {
            df.SIA_matrix[age][year] = 0.0;
        }
    }

    string countryModel = cli.country_model;
    df.countryModel = countryModel;
    string countryISO3 = countryModel.substr(0,3);// Must be country ISO3 code
    df.countryISO3 = countryISO3;
    df.action = cli.action;

    cout << "countryModel: " << countryModel << endl;


    std::stringstream inpss;
    std::stringstream dmss;
    std::stringstream siass;
    // Read demographic data:
    std::ifstream in(cli.demofile);
    if (!in)
    {
        std::cout << "could not open demographic data input file " << cli.demofile << std::endl;
        return(1);
    }
    else
    {
        dmss << readDemographicDataFrame(in,df);
        inpss << dmss.str();
    }
    std::cout <<"["<<__FILE__<<"|"<<__FUNCTION__<<"|"<<__LINE__<<"]" << std::endl;

    // Read SIA data:
    ifstream siain(cli.siafile);
    if (!siain)
    {
        std::cout << "could not open SIA input file " << cli.siafile << std::endl;
        return(1);
    }
    else
    {
siass <<         readAgeSpecificSIAs(siain,df);
inpss <<   siass.str();
    }

    // calculate df + sia_df SHA
    df.dm_sha = get_sha1(dmss.str());
    df.sia_sha = get_sha1(siass.str());
    df.sha = get_sha1(inpss.str());
    std::stringstream sh_ts;
    sh_ts << "_"; 
    sh_ts << df.sha.substr(0, 9);
    sh_ts << "_"; 
    sh_ts << df.timestamp;
    df.shatimestamp = sh_ts.str(); 


    // If model requires discounting MCV2:
    if(atoi(argv[5]) == 1 || atoi(argv[5]) == 3)
    {
        scale_mcv2(df);
    }

    // If model requires discounting SIAs:
    if(atoi(argv[5]) == 2 || atoi(argv[5]) == 3)
    {
        scale_sia(df);
    }

    printDemographicDataFrameToConsole("pre-run",df);
    printDataFrameSIAs(df);

    std::cout <<"["<<__FILE__<<"|"<<__FUNCTION__<<"|"<<__LINE__<<"]"<<  "action:" << cli.action << std::endl;

    // Command should look like the following:
    // ./pfpp 1 NGAnormal Data_for_fitting2019/NGA.txt Data_for_fitting2019/siaNGA.txt 0 0.5 0.01
    if(cli.action == 1)
    {
        Iest_struct Iest = initiateGridSearch(cli,df,countryISO3,countryModel);
        writeIest(Iest,countryModel,0,df.shatimestamp);

        df.nLL =  Iest.nLL;
        df.S0pct =  Iest.S0pct;
        df.p_b0 = Iest.paramsBest[0];
        df.p_b1 = Iest.paramsBest[1];
        df.p_pr = Iest.paramsBest[2];
        df.p_pr_high = Iest.paramsBest[3];
    }

    // This action runs a single parameter set:
    // ./pfpp 5 NGAnormal_5 Data_for_fitting2019/NGA.txt Data_for_fitting2019/siaNGA.txt 0 0.5 0.01 NGAnormal_5.txt
    if(cli.action == 5)
    {
#ifdef STATICSEED
        cout << "WARNING: As this action is intented to produce random replicates, you likely do not" << endl;
        cout << "wish to have STATICSEED enabled." << endl;
        exit(1);
#endif

#ifndef PRINTAGES_BYYEAR
        cout << "WARNING: As this action is intented to produce random replicates, you likely wish to" << endl;
        cout << "have PRINTAGES_BYYEAR defined, to produce the long format outputs." << endl;
        exit(1);
#endif
        runSingleParameterSet(cli,df,countryISO3,countryModel);
    }


    writeDemographicDataFrameToInfoFile(df,cli);

    cout << "EXIT SUCCESS" << endl;
    return(0);
}
