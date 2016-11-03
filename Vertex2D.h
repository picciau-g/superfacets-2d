#ifndef VERTEX2D_H
#define VERTEX2D_H

/**
 * @brief The Vertex2D class represent a point in 2D space
 */

class Vertex2D
{
public:
    ///A constructor method
    Vertex2D();
    ///A constructor method
    Vertex2D(const Vertex2D& orig);
    ///A constructor method
    /*!
     * \param x a float argument, representing the x coordinate
     * \param y a float argument, representing the y coordinate
     */
    Vertex2D(double x, double y);
    ///A destructor method
    virtual ~Vertex2D();
    ///
    friend bool operator== (const Vertex2D& p, const Vertex2D &q);
    ///
    friend bool operator!= (const Vertex2D& p, const Vertex2D &q);
    ///A public method that returns the x coordinate
    /*!
     * \return a double value, representing the x coordinate
     */
    double getX();
    ///A public method that returns the y coordinate
    /*!
     * \return a double value, representing the y coordinate
     */
    double getY();
    ///A public method that sets the x coordinate
    /*!
     * \param x a double argument, represents the value of the x coordinate to set
     */
    void setX(double x);
    ///A public method that sets the x coordinate
    /*!
     * \param x a double argument, represents the value of the x coordinate to set
     */
    void setY(double x);
    ///
    int VTstar();
    ///
    void VTstar(int vtstar);
protected:
    ///A protected variable representing the x coordinate of the point
    double x;
    ///A protected variable representing the y coordinate of the point
    double y;
    ///
    int vtstar;
};

#endif // VERTEX2D_H
