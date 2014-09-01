/*
 *
 *   2014
 *   Author:       Giulia Picciau - DIBRIS, Università degli studi di Genova
 *   Supervisors:  Leila De Floriani - DIBRIS, Università degli studi di Genova
 *                 Patricio Simari - Department of Electrical Engineering and Computer Science, The Catholic University of America
 *
 *   Title:          Fast and scalable mesh superfacets
 *   Submission to Eurographics Symposium on Geometry Processing 2014
 *
 *
 **/
#ifndef NORMALS_H
#define NORMALS_H

#include "Mesh.h"
#include "Vertex3D.h"
#include "Triangle.h"

/**
 *
 * @brief The Normals class stores information about the triangle normals
 *
 */
class Normals
{
public:
    Normals();

    /// Constructor from 3 points
    Normals(Vertex3D, Vertex3D, Vertex3D);

    /// Scalar product
    float dotProd(Normals N);

    /// Vector Normalization
    void Normalize();

    inline float getNx(){
        return this->nx;
    }

    inline float getNy(){
        return this->ny;
    }

    inline float getNz(){
        return this->nz;
    }

    inline void setNX(float a){
        this->nx=a;
    }

    inline void setNY(float a){
        this->ny=a;
    }

    inline void setNZ(float a){
        this->nz=a;
    }

private:
    /// Values for the normal along the 3 axes
    float nx;
    float ny;
    float nz;
};

#endif // NORMALS_H
