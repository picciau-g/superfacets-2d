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

#include "Vertex3D.h"

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
    float DotProd(Normals N);

    /// Vector Normalization
    void Normalize();

    /**
     * @brief getNx
     * @return x component of the normal vector
     */
    inline float getNx()
    {
        return m_Normals.x;
    }

    /**
     * @brief getNy
     * @return y component of the normal vector
     */
    inline float getNy()
    {
        return m_Normals.y;
    }

    /**
     * @brief getNz
     * @return z component of the normal vector
     */
    inline float getNz()
    {
        return m_Normals.z;
    }

    /**
     * @brief setNX set x component of the normal vector
     * @param a value of the x component
     */
    inline void setNX(float pA)
    {
        m_Normals.x = pA;
    }

    /**
     * @brief setNY set y component of the normal vector
     * @param a value of the y component
     */
    inline void setNY(float pA)
    {
        m_Normals.y = pA;
    }

    /**
     * @brief setNZ set z component of the normal vector
     * @param a value of the z component
     */
    inline void setNZ(float pA)
    {
        m_Normals.z = pA;
    }

    inline glm::vec3 GetNormal() const
    {
        return m_Normals;
    }

private:
    /// Values for the normal along the 3 axes
    glm::vec3 m_Normals;
};

#endif // NORMALS_H
