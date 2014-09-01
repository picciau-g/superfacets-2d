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

#include "normals.h"

Normals::Normals()
{
}

/**
 * @brief Normals::Normals calculates normals to a plane (represented with 3 points)
 * @param a first point
 * @param b second point
 * @param c third point
 */
Normals::Normals(Vertex3D a, Vertex3D b, Vertex3D c){

    float vec1[3], vec2[3];

    vec1[0]=b.getX()-a.getX();
    vec1[1]=b.getY()-a.getY();
    vec1[2]=b.getZ()-a.getZ();

    vec2[0]=c.getX()-a.getX();
    vec2[1]=c.getY()-a.getY();
    vec2[2]=c.getZ()-a.getZ();

    this->nx=(vec1[1]*vec2[2])-(vec1[2]*vec2[1]);
    this->ny=(vec1[2]*vec2[0])-(vec1[0]*vec2[2]);
    this->nz=(vec1[0]*vec2[1])-(vec1[1]*vec2[0]);

    //To have them in range [0,1]
    Normalize();
}

/**
 * @brief Normals::Normalize to have normal components in the interval [0,1]
 */
void Normals::Normalize(){
    float normFact=sqrt(nx*nx + ny*ny + nz*nz);

    nx /= normFact;
    ny /= normFact;
    nz /= normFact;
}

/**
 * @brief Normals::dotProd scalar product
 * @param N normal vector
 * @return
 */
float Normals::dotProd(Normals N){

    return this->nx*N.nx + this->ny*N.ny + this->nz*N.nz;
}
