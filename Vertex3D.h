#ifndef _VERTEX3D_H
#define	_VERTEX3D_H

#include <glm/glm.hpp>>
#include <glm/gtx/norm.hpp>

/// An inner-class, extending Vertex2D, representing a vertex in a triangle mesh
class Vertex3D
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

    Vertex3D(const glm::vec3& pCoords);

    ///A destructor method
    virtual ~Vertex3D();
    ///
    friend bool operator== (const Vertex3D& pP, const Vertex3D &pQ);
    ///
    friend bool operator!= (const Vertex3D& pP, const Vertex3D &pQ);


    //manage the VTstar relation (vertex to a face)
    inline void SetVTStar(int pV) {m_VTstar = pV;}
    inline int VTstar() {return m_VTstar;}
    inline void VTstar(int pT) {m_VTstar = pT;}

    ///A public method that returns the x coordinate
    /*!
     * \return a double value, representing the x coordinate
     */
    double getX();
    double getX() const;

    ///A public method that returns the y coordinate
    /*!
     * \return a double value, representing the y coordinate
     */
    double getY();
    double getY() const;


    ///A public method that returns the z coordinate
    /*!
     * \return a double value, representing the z coordinate
     */
    double getZ();
    double getZ() const;

    ///A public method that sets the x coordinate
    /*!
     * \param pX a double argument, represents the value of the x coordinate to set
     */
    void setX(double pX);

    ///A public method that sets the x coordinate
    /*!
     * \param pY a double argument, represents the value of the y coordinate to set
     */
    void setY(double pY);

    ///A public method that sets the x coordinate
    /*!
     * \param pZ a double argument, represents the value of the z coordinate to set
     */
    void setZ(double pZ);

    glm::vec3 GetCoordinates() const;

    /**
     * @brief norma
     * @param pV; A 3D vertex
     * @return The norm of the vector representing the vertex
     */
    double SquaredMagnitude();
    double SquaredDistance(const Vertex3D& pOth);

    /**
     * @brief prodscal
     * @param v1 one vertex
     * @param v2 another vertex
     * @return Scalar product between this-v1 and this-v2
     */
    double prodscal(const Vertex3D& v1,const Vertex3D& v2);

    //scalar products between vectors vec and v1
    double prodscal(const Vertex3D& v1);


    /**
     * @brief distance
     * @param v the Vertex3D we're calculating the distance to
     * @return the distance between this and v
     */
    inline double distance(const Vertex3D& v)
    {
        return glm::distance(m_Coordinates, v.GetCoordinates());
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

        if(m_Coordinates.x < v.getX()) return true;
        else if(m_Coordinates.x == v.getX() && m_Coordinates.y < v.getY()) return true;
        else if(m_Coordinates.x == v.getX() && m_Coordinates.y == v.getY() && m_Coordinates.z < v.getZ()) return true;
        return false;
    }

    /**
     * @brief operator >
     * @param v A vertex
     * @return
     */
    inline bool operator>(const Vertex3D& v){

        if(m_Coordinates.x > v.getX()) return true;
        else if(m_Coordinates.x == v.getX() && m_Coordinates.y > v.getY()) return true;
        else if(m_Coordinates.x == v.getX() && m_Coordinates.y == v.getY() && m_Coordinates.z > v.getZ()) return true;
        return false;
    }

private:
    ///A protected variable representing the m_Z coordinate of the point
    glm::vec3 m_Coordinates;
    double vertexSaliency;
    int m_VTstar;

};

#endif	/* _VERTEX3D_H */

