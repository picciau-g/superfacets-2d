/*
 *
 *   2014
 *   Author:       Giulia Picciau - DIBRIS, Università degli studi di Genova
 *   Supervisors:  Leila De Floriani - DIBRIS, Università degli studi di Genova
 *                 Patricio Simari - Department of Electrical Engineering and Computer Science, The Catholic University of America
 *
 *   Title:          Fast and scalable mesh superfacets
 *   Submission to Pacific Graphics 2014
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

    /**
     * @brief getNx
     * @return x component of the normal vector
     */
    inline float getNx(){
        return this->nx;
    }

    /**
     * @brief getNy
     * @return y component of the normal vector
     */
    inline float getNy(){
        return this->ny;
    }

    /**
     * @brief getNz
     * @return z component of the normal vector
     */
    inline float getNz(){
        return this->nz;
    }

    /**
     * @brief setNX set x component of the normal vector
     * @param a value of the x component
     */
    inline void setNX(float a){
        this->nx=a;
    }

    /**
     * @brief setNY set y component of the normal vector
     * @param a value of the y component
     */
    inline void setNY(float a){
        this->ny=a;
    }

    /**
     * @brief setNZ set z component of the normal vector
     * @param a value of the z component
     */
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
