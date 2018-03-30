#ifndef universe_hpp
#define universe_hpp

#include <vector>
#include "vertex.hpp"
#include "triangle.hpp"
#include "pool.hpp"
#include "bag.hpp"

#define N_TRIANGLES 20000
#define N_VERTICES 10000

class Universe {
    public:
        Pool<Triangle, N_TRIANGLES> triangles;
        Pool<Vertex, N_VERTICES> vertices;

        //Bag<Vertex, N_VERTICES> verticesDel;  //  vertices with coordination number 4. These are candidates for the delete move
        //Bag<Vertex, N_VERTICES> verticesFlip;   //  vertices with more than two upwards pointing links. These can be used in a flip move

        std::vector<int> sliceSizes;
        int nSlices;

        Universe(int n_slices);

        void initialize();

}; 
#endif