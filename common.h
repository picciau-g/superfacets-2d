#ifndef COMMON_H
#define COMMON_H

#include "normals.h"
#include "Timer.h"
#include "Reader.h"
#include <unordered_map>
#include <queue>
#include <unordered_set>
#include <numeric>
#include <QString>
#include <QStringList>
#include <cfloat>
#include <fstream>



/// types for indices of triangles and for adjacencies between them
typedef unsigned long int faceind;
typedef unsigned long long int edgekey;

/// types for grid initialization
typedef signed char coordind;
typedef int gridkey;

/*!
 * \brief The pointDist struct: used for the priority queue in the Dijkstra-based expansion
 */
struct pointDist
{
    int indexP;
    double distanceP;
};

/*!
 * \brief The compare struct: redefinition of the > operator to be able to use a priority queue
 */
struct compare
{

    bool operator ()(pointDist p1, pointDist p2)
    {

        return p1.distanceP > p2.distanceP;
    }
};

/*!
     * \brief getKey
     * \param a index of the first face
     * \param b index of the second face
     * \return the key value (used by the structures storing the distances between faces)
     */
edgekey getKey(faceind a, faceind b)
{

    if(a<=b)
        return (edgekey(a) << 32 | edgekey(b));
    else
        return (edgekey(b) << 32 | edgekey(a));
}


gridkey getGridKey(coordind x, coordind y, coordind z)
{
    //return (gridkey(x) << 16 | gridkey(y) << 8 | gridkey(z));
    if(x <= y && x <= z)
    {
        if(y<=z)
            return (gridkey(fabs(x)) << 16 | gridkey(fabs(y)) << 8 | gridkey(fabs(z)));
        else
            return (gridkey(fabs(x)) << 16 | gridkey(fabs(z)) << 8 | gridkey(fabs(y)));
    }
    else
    {
        if(y <= z)
        {
            if(x <= z)
                return (gridkey(fabs(y)) << 16 | gridkey(fabs(x)) << 8 | gridkey(fabs(z)));
            return (gridkey(fabs(y)) << 16 | gridkey(fabs(z)) << 8 | gridkey(fabs(x)));
        }
        if(x <= y)
            return (gridkey(fabs(z)) << 16 | gridkey(fabs(x)) << 8 | gridkey(fabs(y)));
        return (gridkey(fabs(z)) << 16 | gridkey(fabs(y)) << 8 | gridkey(fabs(x)));
    }

}

int minArrayIndex(float array[3])
{
    if(array[0]<array[1] && array[0]<array[2])
        return 0;
    else if(array[1]<array[2])
        return 1;
    else
        return 2;
}

#endif // COMMON_H
