#include "Vertex2D.h"

/**
 * @brief Vertex2D::Vertex2D class constructor
 */
Vertex2D::Vertex2D()
:
    m_X(0.0)
    , m_Y(0.0)
    , m_VTstar(-1)
{

}

/**
 * @brief Vertex2D::Vertex2D class constructor
 * @param orig origin
 */
Vertex2D::Vertex2D(const Vertex2D& pOrig)
:
    m_X(pOrig.m_X)
    , m_Y(pOrig.m_Y)
    , m_VTstar(pOrig.m_VTstar)
{

}

/**
 * @brief Vertex2D::Vertex2D class constructor
 * @param pX pX coordinate
 * @param pY pY coordinate
 */
Vertex2D::Vertex2D(double pX, double pY)
:
    m_X(pX)
    , m_Y(pY)
    , m_VTstar(-1)
{

}

Vertex2D::~Vertex2D()
{
}

double Vertex2D::getX()
{
    return m_X;
}

double Vertex2D::getY()
{
    return m_Y;
}

double Vertex2D::getX() const
{
    return m_X;
}

double Vertex2D::getY() const
{
    return m_Y;
}


void Vertex2D::setX(double pX)
{
    this->m_X = pX;
}

void Vertex2D::setY(double pY)
{
    this->m_Y = pY;
}

int Vertex2D::VTstar()
{
    return m_VTstar;
}

void Vertex2D::VTstar(int vtstar)
{
    m_VTstar = vtstar;
}

bool operator== (const Vertex2D &p, const Vertex2D &q)
{
    return ((p.m_X == q.m_X) && (p.m_Y == q.m_Y));
}

bool operator !=(const Vertex2D& p, const Vertex2D& q)
{
        return !(p == q);
}
