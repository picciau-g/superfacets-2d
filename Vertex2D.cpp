#include "Vertex2D.h"

/**
 * @brief Vertex2D::Vertex2D class constructor
 */
Vertex2D::Vertex2D()
:
    m_Coordinates(glm::vec2(0.0, 0.0))
{

}

/**
 * @brief Vertex2D::Vertex2D class constructor
 * @param orig origin
 */
Vertex2D::Vertex2D(const Vertex2D& pOrig)
:
    m_Coordinates(pOrig.m_Coordinates)
{

}

/**
 * @brief Vertex2D::Vertex2D class constructor
 * @param pX pX coordinate
 * @param pY pY coordinate
 */
Vertex2D::Vertex2D(double pX, double pY)
:
    m_Coordinates(glm::vec2(pX, pY))
{

}

Vertex2D::~Vertex2D()
{
}

double Vertex2D::getX()
{
    return m_Coordinates.x;
}

double Vertex2D::getY()
{
    return m_Coordinates.y;
}

double Vertex2D::getX() const
{
    return m_Coordinates.x;
}

double Vertex2D::getY() const
{
    return m_Coordinates.y;
}


void Vertex2D::setX(double pX)
{
    m_Coordinates.x = pX;
}

void Vertex2D::setY(double pY)
{
    m_Coordinates.y = pY;
}


bool operator== (const Vertex2D &p, const Vertex2D &q)
{
    return p.m_Coordinates == q.m_Coordinates;//((p.m_X == q.m_X) && (p.m_Y == q.m_Y));
}

bool operator !=(const Vertex2D& p, const Vertex2D& q)
{
        return !(p == q);
}
