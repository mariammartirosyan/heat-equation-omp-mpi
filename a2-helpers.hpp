/**
 *   This file contains a set of helper function. 
 *   Comment: You probably do need to change any of these.
*/
#include <vector>
#include <string>
#include <sstream>
#include <iomanip>
#include <fstream>
#include <iostream>
#include <cmath>

double** allocate(int height, int width, const double& val = 0) {    
    double** ptr = new double*[height]; 
    double* mem = new double[height*width]{ val }; 

    for (unsigned i = 0; i < height; ++i, mem += width)
        ptr[i] = mem;

    return ptr;
}

void deallocate(double** data) {
   delete [] data[0];  
   delete [] data;     
}

/**
 * Matrix data structure
 * 2D => (y-dim, x-dim)
 * Init:
 *   Mat mat(800, 600);

 * Access element at position (34, 53):
 *  mat(34,53)
 * 
 * Access height and width with: 
 * mat.height and mat.width
 * 
*/

struct Mat {
    double** data = nullptr;
    
    int height;
    int width;
    
    
    Mat (int height, int width, const double& val = 0) 
        : height(height), width(width), data(nullptr)
    {   
        if ( height > 0 && width >  0)
            data = allocate(height, width, 0);
    }

    ~Mat() {
        if (data) {
            delete [] data[0];  
            delete [] data; 
        }
    }

    double& operator()(int h, int w) {
        return data[h][w];
    }

    void swap(Mat& right)
    {
        std::swap(this->data, right.data);
        std::swap(this->width, right.width);
        std::swap(this->height, right.height);
    }

    double* operator[](unsigned row)
    {
        return data[row];
    }

    const double* operator[](unsigned row) const
    {
        return data[row];
    }

    /**
     * Checks if all values of two input matrices are equal.
     * 
     * @param[in] a
     * @param[in] b
     * @param[in] epsilon
     * 
     * Return true if two matrices of type Mat have the same values.
    */
    bool compare(Mat& m, double eps=std::pow(10, -8)) {
        if ( this->height != m.height || this->width != m.width ) 
            return false;
        
        for (int j=0;j<this->height;++j) {
            for (int i=0;i<this->width;++i) {
                if ( std::abs( this->operator()(j,i) -  m(j,i)) > eps ) {
                    return false; 
                }
            }
        }

        return true;
    }

    /**
     * Print the Matrix to std out 
     * Could be useful for debugging
    */
    void print() {
        std::cout << this->to_string() << std::endl;
    }

    /**
     * Convert matrix to string 
     * inefficient, but useful for debugging
    */
    std::string to_string() {
        std::stringstream ss;
        for (int j = 0; j < height; ++j) {
            for (int i = 0; i < width; ++i) {
                ss << std::fixed << std::setfill(' ') << std::right << std::setw(11) << std::setprecision(6)  << this->operator()(j,i) << " ";
                
            }
            ss << std::endl;
        }
        
        return ss.str();
    }


    /**
     * Save data to disk in the PPM image format
     * Note: N and M are required arguments
     * 
    */
    void save_to_disk(std::string filename) {
        std::ofstream ofs(filename, std::ofstream::out);
        ofs << this->height << "\n"
            << this->width << "\n";
        
        for (int i = 0; i < this->height; ++i)
        {
            for (int j = 0; j < this->width; ++j)
            {
                ofs << std::fixed << std::setfill(' ') << std::right << std::setw(11) << std::setprecision(6) <<  this->operator()(i,j) << " ";
            }
            ofs << "\n";
        }
        ofs.close();
    }

};

/**
 * Process command line arguments
 * Note: N and M are required arguments
 * 
*/
void process_input(int argc, char **argv, int& N, int& M, int& max_iterations, double& epsilon, bool& verify, bool& print_config) {
    int input_check = 0;

    // process input arguments
    for (int i = 0; i < argc; ++i) {
        if ( std::string(argv[i]).compare("--n") == 0 ) {
            N=atoi(argv[++i]);
            input_check++;
        }
        if ( std::string(argv[i]).compare("--m") == 0 ) {
            M=atoi(argv[++i]);
            input_check++;
        }
        if ( std::string(argv[i]).compare("--no-verify") == 0 ) {
            verify = false;
        }
        if ( std::string(argv[i]).compare("--print-config") == 0 ) {
            print_config = true;
        }
        if ( std::string(argv[i]).compare("--max-iterations") == 0 ) {
            max_iterations=atoi(argv[++i]);
        }
        if ( std::string(argv[i]).compare("--epsilon") == 0 ) {
            epsilon=atof(argv[++i]);
        }
    }
    if ( input_check < 2 ) {
        std::cout << "Usage: --n <columns> --m <rows> [ --max-iterations <int> --epsilon <double> --verify ]" << std::endl;
        exit(-1);
    }
}

/**
 * Jacobi iterative solver for Heat Equation - sequential version
 * @param[inout] U
 * @param[in] max_iterations
 * @param[in] epsilon
 * @param[inout] iteration_count
*/

void heat2d_sequential(Mat& U, int max_iterations, double epsilon, int& iteration_count) {
    int i, j;
    double diffnorm;

    int M = U.height, N = U.width;

    // allocate another 2D array
    Mat W(M,N); 

    // Init & Boundary
    for (i = 0; i < M; ++i) {
        for (j = 0; j < N; ++j) {
            W[i][j] = U[i][j] = 0.0;
        }

        W[i][0] = U[i][0] = 0.05; // left side
        W[i][N-1] = U[i][N-1] = 0.1; // right side
    }

    for (j = 0; j < N; ++j) {
        W[0][j] = U[0][j] = 0.02; // top 
        W[M - 1][j] = U[M - 1][j] = 0.2; // bottom 
    }
    // End init

    int icount = 0;
    do
    {
        icount++;
        diffnorm = 0.0;
        
        /* Compute new values (but not on boundary) */
        for (i = 1; i < M - 1; ++i){
            for (j = 1; j < N - 1; ++j)
            {
                W(i,j) = (U(i,j + 1) + U(i,j - 1) + U(i + 1,j) + U(i - 1,j)) * 0.25;
                diffnorm += (W(i,j) - U(i,j))*(W(i,j) - U(i,j));
            }
        }

        // Only transfer the interior points
        for (i = 1; i < M - 1; ++i)
            for (j = 1; j < N - 1; ++j)
                U(i,j) = W(i,j);

        diffnorm = sqrt(diffnorm);
    } while (epsilon <= diffnorm  && icount < max_iterations);

    iteration_count = icount; // output
}

