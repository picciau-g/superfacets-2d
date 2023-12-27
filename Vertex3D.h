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
    Vertex3D(const Vertex3D& pOrig);
    ///A constructor method
    /*!
     * \param x a float argument, representing the x coordinate
     * \param y a float argument, representing the y coordinate
     * \param z a float argument, representing the z coordinate
     */
    Vertex3D(double pX, double pY, double pZ);
    ///A destructor method
    virtual ~Vertex3D();
    ///
    friend bool operator== (const Vertex3D& pP, const Vertex3D &pQ);
    ///
    friend bool operator!= (const Vertex3D& pP, const Vertex3D &pQ);
    ///A public method that returns the z coordinate
    /*!
     * \return a double value, representing the z coordinate
     */
    double getZ();
    double getZ() const;
    ///A public method that sets the x coordinate
    /*!
     * \param x a double argument, represents the value of the x coordinate to set
     */
    void setZ(double pZ);

    /**
     * @brief norma
     * @param pV; A 3D vertex
     * @return The norm of the vector representing the vertex
     */
    double Norm();
    double Norm(const Vertex3D& pOth);

    /**
     * @brief prodscal
     * @param v1 one vertex
     * @param v2 another vertex
     * @return Scalar product between this-v1 and this-v2
     */
    double prodscal(const Vertex3D& v1,const Vertex3D& v2){return(((v1.getX()-m_X)*(v2.getX()-m_X))+((v1.getY()-m_Y)*(v2.getY()-m_Y))+((v1.getZ()-m_Z)*(v2.getZ()-m_Z)));}

    //scalar products between vectors vec and v1
    double prodscal(const Vertex3D& v1);

    inline double distance(const Vertex3D& v)
    {
        double xdist = this->m_X-v.getX();
        double ydist = this->m_Y-v.getY();
        double zdist = this->m_Z-v.getZ();

        return sqrt(xdist*xdist + ydist*ydist + zdist*zdist);
    }

    /**
     * @brief setSaliency
     * @param sal saliency value to be set
     */
    inline void setSaliency(double sal)
    {
        this->vertexSaliency = sal;
    }


    inline double saliency()
    {
        return this->vertexSaliency;
    }

    /**
     * @brief operator <
     * @param v A vertex
     * @return
     */
    inline bool operator<(const Vertex3D& v){

        if(this->m_X < v.getX()) return true;
        else if(this->m_X == v.getX() && this->m_Y < v.getY()) return true;
        else if(this->m_X == v.getX() && this->m_Y == v.getY() && this->m_Z < v.getZ()) return true;
        return false;
    }

    /**
     * @brief operator >
     * @param v A vertex
     * @return
     */
    inline bool operator>(const Vertex3D& v){

        if(this->m_X > v.getX()) return true;
        else if(this->m_X == v.getX() && this->m_Y > v.getY()) return true;
        else if(this->m_X == v.getX() && this->m_Y == v.getY() && this->m_Z > v.getZ()) return true;
        return false;
    }

private:
    ///A protected variable representing the m_Z coordinate of the point
    double m_Z;
    double vertexSaliency;

};

#endif	/* _VERTEX3D_H */

