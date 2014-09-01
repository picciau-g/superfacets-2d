#ifndef _VERTEX3D_H
#define	_VERTEX3D_H

#include <math.h>
#include "Vertex2D.h"

/// An inner-class, extending Vertex2D, representing a vertex in a triangle mesh
class Vertex3D : public Vertex2D
{
public:
    ///A constructor method
    Vertex3D();
    ///A constructor method
    Vertex3D(const Vertex3D& orig);
    ///A constructor method
    /*!
     * \param x a float argument, representing the x coordinate
     * \param y a float argument, representing the y coordinate
     * \param z a float argument, representing the z coordinate
     */
    Vertex3D(double x, double y, double z);
    ///A destructor method
    virtual ~Vertex3D();
    ///
    friend bool operator== (const Vertex3D& p, const Vertex3D &q);
    ///
    friend bool operator!= (const Vertex3D& p, const Vertex3D &q);
    ///A public method that returns the z coordinate
    /*!
     * \return a double value, representing the z coordinate
     */
    double getZ();
    ///A public method that sets the x coordinate
    /*!
     * \param x a double argument, represents the value of the x coordinate to set
     */
    void setZ(double x);

    //this function returns the norm of vector vec-v
    double norma(Vertex3D& v){return(sqrt(((v.getX()-x)*(v.getX()-x))+((v.getY()-y)*(v.getY()-y))+((v.getZ()-z)*(v.getZ()-z))));}

    //this function returns the scalar products between vectors v1-vec and v2-vec
    double prodscal(Vertex3D& v1,Vertex3D& v2){return(((v1.getX()-x)*(v2.getX()-x))+((v1.getY()-y)*(v2.getY()-y))+((v1.getZ()-z)*(v2.getZ()-z)));}

    //norm of this vector
    double norma() {return sqrt((x)*(x)+(y)*(y)+(z)*(z));}

    //scalar products between vectors vec and v1
    double prodscal(Vertex3D& v1){return(((v1.getX())*(x))+((v1.getY())*(y))+((v1.getZ())*(z)));}

    inline double distance(Vertex3D& v)
    {
        double xdist = this->x-v.getX();
        double ydist = this->y-v.getY();
        double zdist = this->z-v.getZ();

        return sqrt(xdist*xdist + ydist*ydist + zdist*zdist);
    }

    inline void setSaliency(double sal)
    {
        this->vertexSaliency = sal;
    }


    inline double saliency()
    {
        return this->vertexSaliency;
    }

    inline bool operator<(Vertex3D& v){

        if(this->x < v.getX()) return true;
        else if(this->x == v.getX() && this->y < v.getY()) return true;
        else if(this->x == v.getX() && this->y == v.getY() && this->z < v.getZ()) return true;
        return false;
    }

    inline bool operator>(Vertex3D& v){

        if(this->x > v.getX()) return true;
        else if(this->x == v.getX() && this->y > v.getY()) return true;
        else if(this->x == v.getX() && this->y == v.getY() && this->z > v.getZ()) return true;
        return false;
    }

private:    
    ///A protected variable representing the z coordinate of the point
    double z;
    double vertexSaliency;

};

#endif	/* _VERTEX3D_H */

