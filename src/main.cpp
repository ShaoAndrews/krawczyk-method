#include <iostream>
#include <getopt.h>
#include <string>
#include <fstream>
#include <Eigen/Dense>
#include "quadtree.h"
using namespace std;
using namespace Eigen;
template <typename T>
T load_csv(const std::string &path)
{
    std::ifstream in;
    in.open(path);
    std::string line;
    std::vector<double> values;
    uint rows = 0;
    while (std::getline(in, line))
    {
        std::stringstream ss(line);
        std::string cell;
        while (std::getline(ss, cell, ','))
        {
            double val = std::stod(cell);
            values.push_back(val);
        }
        ++rows;
    }
    return Eigen::Map<const Eigen::Matrix<
        typename T::Scalar,
        T::RowsAtCompileTime,
        T::ColsAtCompileTime,
        Eigen::RowMajor>>(values.data(), rows, values.size() / rows);
}

void usage()
{
    std::cout << "Options:" << endl;
    std::cout << " -a,  --f1=f1.csv           A csv file of polynomial about f1" << endl;
    std::cout << " -b,  --f2=f2.csv           A csv file of polynomial about f2" << endl;
    std::cout << " -o,  --out=OUT             The path of output file" << endl;
    std::cout << "                            default value is './curves@f1@f2.csv'" << endl;
    std::cout << "                            The fist line records the time, and " << endl;
    std::cout << "                            the second line records the number of roots" << endl;
    std::cout << " -t,  --tolerance           Set default precision 1e-6" << endl;
    std::cout << " -h,  --help(no value)      Print the message and exit" << endl;
}

int main(int argc, char *argv[])
{
    string f1_path;
    string f2_path;
    double precision = 1e-6;
    string outpath = "./result.json";
    static struct option long_options[] = {
        {"f1", required_argument, 0, 'a'},
        {"f2", required_argument, 0, 'b'},
        {"out", optional_argument, 0, 'o'},
        {"help", no_argument, 0, 'h'},
        {"tolerance", optional_argument, 0, 't'},
        {0, 0, 0, 0}};
    int getopt_ret, option_index;
    bool status = true;
    while (1)
    {
        getopt_ret = getopt_long(argc, argv, "a:b:o::ht::", long_options, &option_index);
        if (getopt_ret == -1)
            break;
        switch (getopt_ret)
        {
        case 0:
            status = false;
            break;
        case 'a':
            f1_path = optarg;
            break;
        case 'b':
            f2_path = optarg;
            break;
        case 'o':
            outpath = optarg;
            break;
        case 'h':
            usage();
            status = false;
            break;
        case 't':
            precision = atof(optarg);
            break;
        case '?':
            status = false;
            break;
        default:
            status = false;
            break;
        }
    }
    if (status)
    {
        if (access((f1_path).c_str(), F_OK) != 0)
        {
            std::cout << "The path '" << f1_path << "' does not exist, pleasure check first" << std::endl;
            return -1;
        }
        if (access((f2_path).c_str(), F_OK) != 0)
        {
            std::cout << "The path '" << f2_path << "' does not exist, pleasure check first" << std::endl;
            return -1;
        }
        Eigen::MatrixXd f1_mat = load_csv<Eigen::MatrixXd>(f1_path);
        Eigen::MatrixXd f2_mat = load_csv<Eigen::MatrixXd>(f2_path);
        IntervalVec<itv> range(itv(0, 1), itv(0, 1));
        quadtree * solver = new quadtree(f1_mat, f2_mat, range);
        solver->search_procedure(range);
        delete solver;
        solver = nullptr;
    }
}
