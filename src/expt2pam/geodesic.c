#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "geodesic.h"

//
//  geodesic
//  geometry builder
//  Copyright (c) 2013 Robby Kraft
//  MIT open source license
//

////////////////////////////////////////////////////////////////////
//                        PLATONIC SOLIDS                         //
//           centered at origin          radius = 1               //
//                                                                //
//         TETRA               OCTA              ICOSA            //
//             one point aligned along the +X axis                //
//                                                                //
//                   HEXA                DODECA                   //
//          aligned as duals to the octa / icosahedron            //
//                   faces in place of points                     //
//         (one face normal aligned along the +X axis)            //
//                                                                //
//                     + + + + + + + + + + +                      //
//                                                                //
//     even possible to cross-reference indices across duals      //
//                                                                //
//  eg:                                                           //
//   octa's 6 point indices correlate to hexa's 6 face indices    //
//       and octa's point[0] is along hexa's face[0] normal       //
//                                                                //
//      this means: a solid's face normal is that same index      //
//                  in its dual's point array                     //
//                                                                //
////////////////////////////////////////////////////////////////////

// for higher precision replace occurrences of "float" with "double" or "long double"



const double _tetrahedron_points[TETRAHEDRON_POINT_COUNT * 3] = {
        1.0, 0.0, 0.0,
        -0.333333333333333, -0.942809041582063, 0.0,
        -0.333333333333333, 0.471404520791032, 0.816496580927726,
        -0.333333333333333, 0.471404520791032, -0.816496580927726};
const unsigned short _tetrahedron_lines[TETRAHEDRON_LINE_COUNT * 2] = {
        2, 3, 2, 0, 2, 1, 3, 0, 3, 1, 0, 1};
const unsigned short _tetrahedron_faces[TETRAHEDRON_FACE_COUNT * 3] = {
        2, 1, 3,
        2, 3, 0,
        2, 0, 1,
        3, 1, 0};
const double _tetrahedron_dual_points[TETRAHEDRON_POINT_COUNT * 3] = {
        -1.0, 0.0, 0.0,
        0.333333333333333, 0.942809041582063, 0.0,
        0.333333333333333, -0.471404520791032, 0.816496580927726,
        0.333333333333333, -0.471404520791032, -0.816496580927726};
const unsigned short _tetrahedron_dual_lines[TETRAHEDRON_LINE_COUNT * 2] = {
        2, 3, 2, 0, 2, 1, 3, 0, 3, 1, 0, 1};
const unsigned short _tetrahedron_dual_faces[TETRAHEDRON_FACE_COUNT * 3] = {
        2, 1, 3,
        2, 3, 0,
        2, 0, 1,
        3, 1, 0};
const double _octahedron_points[OCTAHEDRON_POINT_COUNT * 3] = {
        1.0, 0.0, 0.0,
        0.0, 1.0, 0.0,
        0.0, 0.0, 1.0,
        -1.0, 0.0, 0.0,
        0.0, -1.0, 0.0,
        0.0, 0.0, -1.0};
const unsigned short _octahedron_lines[OCTAHEDRON_LINE_COUNT * 2] = {
        1, 0, 1, 2, 1, 5, 1, 3, 3, 2, 2, 0, 0, 5, 5, 3, 4, 2, 4, 3, 4, 5, 4, 0};
const unsigned short _octahedron_faces[OCTAHEDRON_FACE_COUNT * 3] = {
        1, 0, 2,
        1, 5, 0,
        4, 0, 5,
        4, 2, 0,
        1, 2, 3,
        1, 3, 5,
        4, 5, 3,
        4, 3, 2};
const double _hexahedron_points[HEXAHEDRON_POINT_COUNT * 3] = {
        0.57735026918963, 0.57735026918963, 0.57735026918963,
        0.57735026918963, 0.57735026918963, -0.57735026918963,
        0.57735026918963, -0.57735026918963, -0.57735026918963,
        0.57735026918963, -0.57735026918963, 0.57735026918963,
        -0.57735026918963, 0.57735026918963, 0.57735026918963,
        -0.57735026918963, 0.57735026918963, -0.57735026918963,
        -0.57735026918963, -0.57735026918963, -0.57735026918963,
        -0.57735026918963, -0.57735026918963, 0.57735026918963};
const unsigned short _hexahedron_lines[HEXAHEDRON_LINE_COUNT * 2] = {
        0, 1, 1, 2, 2, 3, 3, 0, 0, 4, 1, 5, 2, 6, 3, 7, 4, 5, 5, 6, 6, 7, 7, 4};
const unsigned short _hexahedron_faces[HEXAHEDRON_FACE_COUNT * 4] = {
// these are being stored as squares, not triangles
        0, 3, 2, 1,
        4, 0, 1, 5,
        0, 3, 7, 4,
        7, 4, 5, 6,
        3, 7, 6, 2,
        1, 5, 6, 2};
const unsigned short _hexahedron_triangle_faces[HEXAHEDRON_TRIANGLE_FACE_COUNT * 3] = {
        0, 2, 3,
        2, 0, 1,
        4, 1, 0,
        1, 4, 5,
        0, 3, 7,
        7, 4, 0,
        7, 5, 4,
        5, 7, 6,
        3, 6, 7,
        6, 3, 2,
        1, 5, 6,
        6, 2, 1};
const double _icosahedron_points[ICOSAHEDRON_POINT_COUNT * 3] = {
        0.447213595499958, -0.276393202250021, 0.850650808352040,
        -0.447213595499958, 0.276393202250021, 0.850650808352040,
        -0.447213595499958, 0.276393202250021, -0.850650808352040,
        0.447213595499958, -0.276393202250021, -0.850650808352040,
        -0.447213595499958, -0.723606797749979, 0.525731112119134,
        0.447213595499958, 0.723606797749979, 0.525731112119134,
        0.447213595499958, 0.723606797749979, -0.525731112119134,
        -0.447213595499958, -0.723606797749979, -0.525731112119134,
        0.447213595499958, -0.894427190999916, 0.0,
        1.0, 0.0, 0.0,
        -0.447213595499958, 0.894427190999916, 0.0,
        -1.0, 0.0, 0.0};
const unsigned short _icosahedron_lines[ICOSAHEDRON_LINE_COUNT * 2] = {
        0, 8, 0, 9, 0, 1, 0, 4, 0, 5, 8, 3, 8, 9, 8, 7, 8, 4, 9, 3,
        9, 6, 9, 5, 7, 4, 7, 3, 7, 2, 7, 11, 2, 10, 2, 11, 2, 3, 2, 6,
        10, 11, 10, 5, 10, 6, 10, 1, 11, 1, 11, 4, 4, 1, 5, 1, 5, 6, 6, 3};
const unsigned short _icosahedron_faces[ICOSAHEDRON_FACE_COUNT * 3] = {
        8, 7, 4,
        7, 8, 3,    // pair 2
        8, 4, 0,
        8, 0, 9,
        9, 3, 8,
        9, 0, 5,
        9, 5, 6,
        9, 6, 3,
        3, 2, 7,   // pair 2
        3, 6, 2,
        0, 4, 1,
        0, 1, 5,
        11, 4, 7,
        11, 7, 2,
        11, 2, 10,
        1, 11, 10,    // pair 1
        11, 1, 4,
        10, 6, 5,
        10, 5, 1,     // pair 1
        10, 2, 6};
const double _dodecahedron_points[DODECAHEDRON_POINT_COUNT * 3] = {
        -0.794655, 0.491123, 0.356822,
        -0.187593, 0.794654, -0.577350,
        -0.794655, 0.491123, -0.356822,
        -0.794654, -0.187593, -0.577350,
        0.187592, 0.303531, -0.934172,
        0.187592, 0.982247, 0.000000,
        0.187593, -0.794654, -0.577350,
        -0.187592, -0.303531, -0.934172,
        -0.187592, -0.982247, 0.000000,
        -0.794654, -0.607062, 0.000000,
        0.794655, -0.491123, -0.356822,
        0.187593, -0.794654, 0.577350,
        -0.187592, -0.303531, 0.934172,
        -0.794654, -0.187593, 0.577350,
        0.794655, -0.491123, 0.356822,
        -0.187593, 0.794654, 0.577350,
        0.187592, 0.303531, 0.934172,
        0.794654, 0.187593, 0.577350,
        0.794654, 0.607062, 0.000000,
        0.794654, 0.187593, -0.577350};
const unsigned short _dodecahedron_lines[DODECAHEDRON_LINE_COUNT * 2] = {
        19, 18,
        18, 17,
        17, 14,
        14, 10,
        10, 19,
        14, 11,
        11, 8,
        8, 6,
        6, 10,
        19, 4,
        4, 1,
        1, 5,
        5, 18,
        9, 3,
        3, 7,
        7, 6,
        8, 9,
        13, 9,
        3, 2,
        2, 0,
        0, 13,
        12, 11,
        12, 13,
        1, 2,
        5, 15,
        15, 0,
        15, 16,
        16, 17,
        16, 12,
        4, 7};
const unsigned short _dodecahedron_triangle_faces[DODECAHEDRON_TRIANGLE_FACE_COUNT * 3] = {
        5, 2, 1,
        15, 0, 5,
        2, 5, 0,

        9, 2, 13,
        2, 0, 13,
        3, 2, 9,

        10, 17, 19,
        18, 19, 17,
        14, 17, 10,

        19, 1, 4,
        6, 9, 8,
        1, 18, 5,

        11, 14, 8,
        18, 1, 19,
        8, 9, 11,

        10, 8, 14,
        8, 10, 6,
        6, 7, 9,

        9, 7, 3,
        13, 12, 9,
        12, 11, 9,

        18, 15, 5,
        17, 15, 18,
        15, 17, 16,

        11, 12, 14,
        14, 12, 17,
        17, 12, 16,

        12, 13, 16,
        13, 0, 16,
        0, 15, 16,

        7, 2, 3,
        2, 7, 1,
        1, 7, 4,

        4, 7, 19,
        19, 7, 10,
        6, 10, 7};

/**********************************************************************/
/**                         geodesic polyhedra                       **/
/**********************************************************************/


// SPHERE
void _divide_geodesic_faces(geodesicSphere *g, int v);
void _remove_duplicate_points_lines(geodesicSphere *g);
void _spherize_points(flo_t *points, unsigned int numPoints);
void _apply_geodesic_sphere_normals(geodesicSphere *g);


void tetrahedronMesh(flo_t **po, unsigned int *numPoints,
                     unsigned short **li, unsigned int *numLines,
                     unsigned short **fa, unsigned int *numFaces)
{
    *numPoints = TETRAHEDRON_POINT_COUNT;
    *numLines = TETRAHEDRON_LINE_COUNT;
    *numFaces = TETRAHEDRON_FACE_COUNT;
    flo_t *points = malloc(sizeof(flo_t) * (*numPoints) * 3);
    unsigned short *lines = malloc(sizeof(unsigned short) * (*numLines) * 2);
    unsigned short *faces = malloc(sizeof(unsigned short) * (*numFaces) * 3);
    for (int i = 0; i < (*numPoints) * 3; i++) {
        points[i] = _tetrahedron_points[i];
    }
    for (int i = 0; i < (*numLines) * 2; i++) {
        lines[i] = _tetrahedron_lines[i];
    }
    for (int i = 0; i < (*numFaces) * 3; i++) {
        faces[i] = _tetrahedron_faces[i];
    }
    *po = points;
    *li = lines;
    *fa = faces;
}


void octahedronMesh(flo_t **po, unsigned int *numPoints,
                    unsigned short **li, unsigned int *numLines,
                    unsigned short **fa, unsigned int *numFaces)
{
    *numPoints = OCTAHEDRON_POINT_COUNT;
    *numLines = OCTAHEDRON_LINE_COUNT;
    *numFaces = OCTAHEDRON_FACE_COUNT;
    flo_t *points = malloc(sizeof(flo_t) * (*numPoints) * 3);
    unsigned short *lines = malloc(sizeof(unsigned short) * (*numLines) * 2);
    unsigned short *faces = malloc(sizeof(unsigned short) * (*numFaces) * 3);
    for (int i = 0; i < (*numPoints) * 3; i++) {
        points[i] = _octahedron_points[i];
    }
    for (int i = 0; i < (*numLines) * 2; i++) {
        lines[i] = _octahedron_lines[i];
    }
    for (int i = 0; i < (*numFaces) * 3; i++) {
        faces[i] = _octahedron_faces[i];
    }
    *po = points;
    *li = lines;
    *fa = faces;
}


void icosahedronMesh(flo_t **po, unsigned int *numPoints,
                     unsigned short **li, unsigned int *numLines,
                     unsigned short **fa, unsigned int *numFaces)
{
    *numPoints = ICOSAHEDRON_POINT_COUNT;
    *numLines = ICOSAHEDRON_LINE_COUNT;
    *numFaces = ICOSAHEDRON_FACE_COUNT;
    flo_t *points = malloc(sizeof(flo_t) * (*numPoints) * 3);
    unsigned short *lines = malloc(sizeof(unsigned short) * (*numLines) * 2);
    unsigned short *faces = malloc(sizeof(unsigned short) * (*numFaces) * 3);
    for (int i = 0; i < (*numPoints) * 3; i++) {
        points[i] = _icosahedron_points[i];
    }
    for (int i = 0; i < (*numLines) * 2; i++) {
        lines[i] = _icosahedron_lines[i];
    }
    for (int i = 0; i < (*numFaces) * 3; i++) {
        faces[i] = _icosahedron_faces[i];
    }
    *po = points;
    *li = lines;
    *fa = faces;
}


geodesicSphere tetrahedronSphere(unsigned int v)
{
    geodesicSphere g;
    g.frequency = v;
    tetrahedronMesh(&g.points, &g.numPoints, &g.lines, &g.numLines, &g.faces, &g.numFaces);
    _divide_geodesic_faces(&g, v);
    _spherize_points(g.points, g.numPoints);
    _apply_geodesic_sphere_normals(&g);
    return g;
}


geodesicSphere octahedronSphere(unsigned int v)
{
    geodesicSphere g;
    g.frequency = v;
    octahedronMesh(&g.points, &g.numPoints, &g.lines, &g.numLines, &g.faces, &g.numFaces);
    _divide_geodesic_faces(&g, v);
    _spherize_points(g.points, g.numPoints);
    _apply_geodesic_sphere_normals(&g);
    return g;
}


geodesicSphere icosahedronSphere(unsigned int v)
{
    geodesicSphere g;
    g.frequency = v;
    icosahedronMesh(&g.points, &g.numPoints, &g.lines, &g.numLines, &g.faces, &g.numFaces);
    _divide_geodesic_faces(&g, v);

    g.pointsNotSpherized = malloc(sizeof(flo_t) * g.numPoints * 3);
    memcpy(g.pointsNotSpherized, g.points, sizeof(flo_t) * g.numPoints * 3);

    _spherize_points(g.points, g.numPoints);

    g.pointsDeltaSpherized = malloc(sizeof(flo_t) * g.numPoints * 3);
    for (int i = 0; i < g.numPoints * 3; i++) {
        g.pointsDeltaSpherized[i] = g.points[i] - g.pointsNotSpherized[i];
    }

    _apply_geodesic_sphere_normals(&g);
    return g;
}


void deleteGeodesicSphere(geodesicSphere *g)
{
    // be careful with this one:
    // an initially unallocated geodesic will still register
    // TRUE on the if()s and call free() and crash
    g->frequency = 0;
    g->numPoints = 0;
    g->numLines = 0;
    g->numFaces = 0;
    if (g->points) {
        free(g->points);
        g->points = NULL;
    }
    if (g->lines) {
        free(g->lines);
        g->lines = NULL;
    }
    if (g->faces) {
        free(g->faces);
        g->faces = NULL;
    }
    if (g->pointNormals) {
        free(g->pointNormals);
        g->pointNormals = NULL;
    }
    if (g->lineNormals) {
        free(g->lineNormals);
        g->lineNormals = NULL;
    }
    if (g->faceNormals) {
        free(g->faceNormals);
        g->faceNormals = NULL;
    }
}


void _apply_geodesic_sphere_normals(geodesicSphere *g)
{
    // shortcuts are made possible due to
    // - all points lying on the surface of a sphere
    // - centered at the origin
    if (g->numPoints) {
//        flo_t length = 1.0;  // shortcut
        g->pointNormals = malloc(sizeof(flo_t) * g->numPoints * 3);
        for (int i = 0; i < g->numPoints; i++) {
//            length = sqrtf( pow(g->points[0+3*i],2) + pow(g->points[1+3*i],2) + pow(g->points[2+3*i],2) ); // shortcut: radius of the sphere
            g->pointNormals[0 + 3 * i] = g->points[0 + 3 * i];//  / length;// * (i/12.0);
            g->pointNormals[1 + 3 * i] = g->points[1 + 3 * i];//  / length;// * (i/12.0);
            g->pointNormals[2 + 3 * i] = g->points[2 + 3 * i];//  / length;// * (i/12.0);
        }
    }
    if (g->numLines) {
        g->lineNormals = malloc(sizeof(flo_t) * g->numLines * 3);
        for (int i = 0; i < g->numLines; i++) {
            g->lineNormals[i * 3 + 0] = (g->pointNormals[g->lines[i * 2 + 0] * 3 + 0] +
                                         g->pointNormals[g->lines[i * 2 + 1] * 3 + 0]) / 2.0;
            g->lineNormals[i * 3 + 1] = (g->pointNormals[g->lines[i * 2 + 0] * 3 + 1] +
                                         g->pointNormals[g->lines[i * 2 + 1] * 3 + 1]) / 2.0;
            g->lineNormals[i * 3 + 2] = (g->pointNormals[g->lines[i * 2 + 0] * 3 + 2] +
                                         g->pointNormals[g->lines[i * 2 + 1] * 3 + 2]) / 2.0;
        }
    }
    if (g->numFaces) {
        g->faceNormals = malloc(sizeof(flo_t) * g->numFaces * 3);
        for (int i = 0; i < g->numFaces; i++) {
            g->faceNormals[i * 3 + 0] = (g->pointNormals[g->faces[i * 3 + 0] * 3 + 0] +
                                         g->pointNormals[g->faces[i * 3 + 1] * 3 + 0] +
                                         g->pointNormals[g->faces[i * 3 + 2] * 3 + 0]) / 3.0;
            g->faceNormals[i * 3 + 1] = (g->pointNormals[g->faces[i * 3 + 0] * 3 + 1] +
                                         g->pointNormals[g->faces[i * 3 + 1] * 3 + 1] +
                                         g->pointNormals[g->faces[i * 3 + 2] * 3 + 1]) / 3.0;
            g->faceNormals[i * 3 + 2] = (g->pointNormals[g->faces[i * 3 + 0] * 3 + 2] +
                                         g->pointNormals[g->faces[i * 3 + 1] * 3 + 2] +
                                         g->pointNormals[g->faces[i * 3 + 2] * 3 + 2]) / 3.0;
        }
    }
}

// NEW POINTS / FACE
// V0: 1           +2 =
// V1: 3 per face  +3 =
// V2: 6           +4 =
// V3: 10          +5 =
// V4: 15          +6 =
// V5: 21

// NEW LINES / FACE
// V0: 0    +3=   (3*1)
// V1: 3    +6=   (3*2)
// V2: 9    +9=   (3*3)
// V3: 18   +12=  (3*4)
// V4: 30   +15=  (3*5)

// NEW FACES / FACE
// V0:
// V1: 1  +3 =  (1*1)
// V2: 4  +5 =  (2*2)
// V3: 9  +7 =  (3*3)
// V4: 16 +9 =  (4*4)
// V5: 25       (5*5)

//         o   A                               0
//        / \         clockwise winding       2 1
//       /   \                               5 4 3
//      /     \                             9 8 7 6
//  C  o_______o  B
//
//           side length  /  frequency (v)
//       \  = AB     short for AB/v
//        / = AC
//
//        /\         /\        /\           /\
//       /\/\                   /\         /\/\
//      /\/\/\
//     /\/\/\/\
//                 first    only 1 up      add 2
//       goal       row     triangle     up and down
//
//        /\                /\              /\
//       /\/\              /\/\            /\/\
//      /\/\/\            /\/\/\          /\/\/\
//           /\              /\/\          /\/\/\
//
//     after the       add 2 faces
//     first face        each step
//    of each row     an up and a down
//          ______
//       a /\    / b    INDE0 SHORTCUTS (row = row # from A)
//        /  \  /        p->c  i-1
//      p/____\/c        p->a  i-row
//                       p->b  i-row-1

// fills data into geodesic.points, and geodesic.pointsNotSpherized
void _divide_geodesic_faces(geodesicSphere *g, int v)
{
    if (v > 1) {
        // calculate new points per face
        // V0: 1           +2 =
        // V1: 3 per face  +3 =
        // V2: 6           +4 =
        // V3: 10          +5 =
        // V4: 15          +6 =
        // V5: 21
        int pointsPerFace = 3;
        for (int i = 1; i < v; i++) {
            pointsPerFace += (i + 2);
        }
        // calculate new lines per face
        // V0: 0    +3=   (3*1)
        // V1: 3    +6=   (3*2)
        // V2: 9    +9=   (3*3)
        // V3: 18   +12=  (3*4)
        // V4: 30   +15=  (3*5)
        int linesPerFace = 3;
        for (int i = 1; i < v; i++) {
            linesPerFace += 3 * (i + 1);
        }
        // new Points, Faces arrays, and their sizes
        flo_t newPointsArray[g->numFaces * pointsPerFace * 3 + g->numPoints];
        // legacy data
        g->parentFace = malloc(sizeof(unsigned short) * g->numFaces * pointsPerFace + g->numPoints);

        unsigned short newFacesArray[
                v * (v + 1) * g->numFaces * 3 * 3];   // data overflow problem. TODO: correctly approximate array size
        unsigned short newLinesArray[linesPerFace * g->numFaces * 2];
        // incrementers for the new arrays as we increment and add to them
        unsigned int newPI = 0;
        unsigned int newLI = 0;
        unsigned int newFI = 0;
        // original points in their original indices
        for (int i = 0; i < g->numPoints; i++) {
            newPointsArray[i * 3 + 0] = g->points[i * 3 + 0];
            newPointsArray[i * 3 + 1] = g->points[i * 3 + 1];
            newPointsArray[i * 3 + 2] = g->points[i * 3 + 2];
            // legacy
            g->parentFace[newPI] = -1; // edge vertices aren't a member of only one face
            newPI++;
        }
        // bring along the parent polyhedra's faces too
        // makes for interesting non-convex geometry
        // othewise, leave this commented out
//        for(int i = 0; i < numFaces; i++){
//            newFacesArray[i*3+0] = faces[i*3+0];
//            newFacesArray[i*3+1] = faces[i*3+1];
//            newFacesArray[i*3+2] = faces[i*3+2];
//            newFI++;
//        }
        int i, j, k;
        // TODO this can probably remain an int, and just use the variable v
        float segments = v;
        // the 3 vertices of the parent triangle we will subdivide
        int faceEdgeA, faceEdgeB, faceEdgeC;
        flo_t *edgePointA, *edgePointB, *edgePointC;
        // vectors: line segments AB and BC divided by the frequency number
        flo_t dAB[3], dBC[3];
        // increment through the original set of faces
        for (i = 0; i < g->numFaces; i++) {
            // save the original major 3 vertices
            faceEdgeA = g->faces[i * 3 + 0];
            faceEdgeB = g->faces[i * 3 + 1];
            faceEdgeC = g->faces[i * 3 + 2];
            edgePointA = &g->points[faceEdgeA * 3];
            edgePointB = &g->points[faceEdgeB * 3];
            edgePointC = &g->points[faceEdgeC * 3];
            // calculate the vector quantity from line segment A to B divided by frequency
            // same with B to C
            dAB[0] = (edgePointB[0] - edgePointA[0]) / segments;
            dAB[1] = (edgePointB[1] - edgePointA[1]) / segments;
            dAB[2] = (edgePointB[2] - edgePointA[2]) / segments;
            dBC[0] = (edgePointC[0] - edgePointB[0]) / segments;
            dBC[1] = (edgePointC[1] - edgePointB[1]) / segments;
            dBC[2] = (edgePointC[2] - edgePointB[2]) / segments;
            // starting at point A, begin generating points one row at a time
            // incrementing towards line segment BC
            // iterate 1, 12, 123, 1234, 12345, 123456...
            for (j = 0; j <= v; j++) {
                for (k = 0; k <= j; k++) {
                    // skip the 3 original vertices
                    if (!((j == 0 && k == 0) || (j == v & k == 0) ||
                          (j == v && k == v))) {  //ignore 3 points of the triangle
                        // LEGACY
                        g->parentFace[newPI] = i;
                        // POINTS
                        newPointsArray[newPI * 3 + 0] = edgePointA[0] + j * dAB[0] + k * dBC[0];
                        newPointsArray[newPI * 3 + 1] = edgePointA[1] + j * dAB[1] + k * dBC[1];
                        newPointsArray[newPI * 3 + 2] = edgePointA[2] + j * dAB[2] + k * dBC[2];
                        newPI++;
                    }
                    // FACES and LINES
                    if (k != 0) {
                        // build a vertical pointing triangle face
                        int faceP1 = (newPI - 1);
                        int faceP2 = (newPI - 1) - j - 1;
                        int faceP3 = (newPI - 1) - 1;
                        if (j == v) {
                            faceP2++;
                        }  // last row->parent row is offset by one because of skipping one of the original triangle points
                        if ((j == v && k ==
                                       v)) {  // i have no idea why the last triangle on the last row differs from the other triangles on the last row
                            faceP2++;
                            faceP3++;
                        }
                        // why we save original face edges:
                        // instead of generating new edge vertices,
                        // preserve 3 original edge vertices
                        if (j == 1) faceP2 = faceEdgeA;  // (original pointA)
                        if (j == v && k == 1) faceP3 = faceEdgeB;
                        if (j == v && k == v) faceP1 = faceEdgeC;

                        newFacesArray[newFI * 3 + 0] = faceP1;
                        newFacesArray[newFI * 3 + 1] = faceP2;
                        newFacesArray[newFI * 3 + 2] = faceP3;
                        newFI++;

                        // LINES
                        // from vertically pointed triangles
                        newLinesArray[newLI * 2 + 0] = faceP1;
                        newLinesArray[newLI * 2 + 1] = faceP2;
                        newLI++;
                        newLinesArray[newLI * 2 + 0] = faceP2;
                        newLinesArray[newLI * 2 + 1] = faceP3;
                        newLI++;
                        newLinesArray[newLI * 2 + 0] = faceP3;
                        newLinesArray[newLI * 2 + 1] = faceP1;
                        newLI++;
                        // LINES END

                        //also build a downward pointing triangle face
                        if (k != j) {
                            faceP1 = (newPI - 1);
                            faceP2 = (newPI - 1) - j + 1 - 1;
                            faceP3 = (newPI - 1) - j - 1;
                            if (j == v) {
                                faceP2++;
                                faceP3++;
                            }
                            newFacesArray[newFI * 3 + 0] = faceP1;
                            newFacesArray[newFI * 3 + 1] = faceP2;
                            newFacesArray[newFI * 3 + 2] = faceP3;
                            newFI++;
                        }
                    }
                }
            }
        }

        g->numPoints = newPI;
        free(g->points);
        g->points = malloc(sizeof(flo_t) * g->numPoints * 3);
        memcpy(g->points, newPointsArray, sizeof(flo_t) * g->numPoints * 3);
        // for(int i = 0; i < g->numPoints*3; i++)
        //     g->points[i] = newPointsArray[i];

        g->numLines = newLI;
        free(g->lines);
        g->lines = malloc(sizeof(unsigned short) * g->numLines * 2);
        memcpy(g->lines, newLinesArray, sizeof(unsigned short) * g->numLines * 2);
        // for(int i = 0; i < g->numLines*2; i++)
        //     g->lines[i] = newLinesArray[i];

        g->numFaces = newFI;
        free(g->faces);
        g->faces = malloc(sizeof(unsigned short) * g->numFaces * 3);
        memcpy(g->faces, newFacesArray, sizeof(unsigned short) * g->numFaces * 3);
        // for(int i = 0; i < g->numFaces*3; i++)
        // g->faces[i] = newFacesArray[i];


        // what is left missing, is there are duplicate points along the original
        // face edge lines, due to subdividing each original triangle face
        // without being aware of which faces are its neighbors.

        // an un-elegant fix is to heuristically merge points
        // that have the same coordinates into one point
        // and update pointers in lines[] and faces[] arrays
        _remove_duplicate_points_lines(g);
    }
}


void _spherize_points(flo_t *points, unsigned int numPoints)
{
    int i;
    flo_t difference, distance;
    flo_t maxdistance = 1.0;//sqrt( ((1 + sqrt(5)) / 2 ) + 2 );
    for (i = 0; i < numPoints; i++) {
        distance = sqrt(pow(points[i * 3 + 0], 2) +
                        pow(points[i * 3 + 1], 2) +
                        pow(points[i * 3 + 2], 2));
        difference = maxdistance / distance;
        points[i * 3 + 0] *= difference;
        points[i * 3 + 1] *= difference;
        points[i * 3 + 2] *= difference;
    }
}



// sample 128 precision: 1.189731495357231765085759326628007
// sample 64 precision: 1.7976931348623157
// sample 32 precision: 3.4028234
#if _float_tprecision == 128
#define ELBOW .0000000000000001
#elif _float_tprecision == 64
#define ELBOW .00000000001
#elif _float_tprecision == 32
#define ELBOW .00001
#else
#define ELBOW .00001
#endif

// subdividing faces without face-neighbor data allows more freedom
// for the algorithm to work on many objects, but requires more work:
// for each set of joined faces duplicate points will be generated along their shared line

void _remove_duplicate_points_lines(geodesicSphere *g)
{

    // make array of size numPoints which looks like this:
    // -1  -1  -1  -1   3  -1  -1  -1  -1  -1  5   5  -1  -1  -1
    // mostly -1s, except at duplicate points, store the number of the index of the first instance of duplication
    int duplicateIndexes[g->numPoints];
    for (int i = 0; i < g->numPoints; i++) {
        duplicateIndexes[i] = -1;
    }
    for (int i = 0; i < g->numPoints - 1; i++) {
        for (int j = i + 1; j < g->numPoints; j++) {
            if (g->points[0 + i * 3] - ELBOW < g->points[0 + j * 3] &&
                g->points[0 + i * 3] + ELBOW > g->points[0 + j * 3] &&
                g->points[1 + i * 3] - ELBOW < g->points[1 + j * 3] &&
                g->points[1 + i * 3] + ELBOW > g->points[1 + j * 3] &&
                g->points[2 + i * 3] - ELBOW < g->points[2 + j * 3] &&
                g->points[2 + i * 3] + ELBOW > g->points[2 + j * 3]) {
                duplicateIndexes[j] = i;
            }
        }
    }

    // replaces all pointers to duplicated indexes with their first instance
    // FACES
    for (int f = 0; f < g->numFaces * 3; f++) {
        if (duplicateIndexes[g->faces[f]] != -1) {
            g->faces[f] = duplicateIndexes[g->faces[f]];
        }
    }
    // LINES
    unsigned short *lineWasAssociatedWithADuplicate = calloc(g->numLines, sizeof(unsigned short));
    for (int l = 0; l < g->numLines * 2; l++) {
        if (duplicateIndexes[g->lines[l]] != -1) {
            g->lines[l] = duplicateIndexes[g->lines[l]];
            lineWasAssociatedWithADuplicate[l / 2] = 1;  // this is going to help us with searching our line duplicates
        }
    }

    //
    //   DUPLICATE LINES
    //
    // now we have all we need to handle duplicate account of lines
    // build an array of -1s, except where a duplicate lies-
    // it will contain the index of it's first occurrence
    int duplicateLineIndexes[g->numLines];
    for (int i = 0; i < g->numLines; i++) {
        duplicateLineIndexes[i] = -1;
    }
    unsigned int duplicateCount = 0;
    for (int i = 0; i < g->numLines; i++) {
        for (int j = i + 1; j < g->numLines; j++) {
            // loop in a loop, bad news
            // use the following to cut down on calls
            if (lineWasAssociatedWithADuplicate[j]) {
                if ((g->lines[i * 2 + 0] == g->lines[j * 2 + 0] && g->lines[i * 2 + 1] == g->lines[j * 2 + 1]) ||
                    (g->lines[i * 2 + 0] == g->lines[j * 2 + 1] && g->lines[i * 2 + 1] == g->lines[j * 2 + 0])) {
                    if (duplicateLineIndexes[j] == -1) {
                        duplicateLineIndexes[j] = i;
                        duplicateCount++;
                    }
                }
            }
        }
    }
    free(lineWasAssociatedWithADuplicate);
    unsigned int newNumLines = g->numLines - duplicateCount;

    unsigned int indexLineOffset = 0;
    // invert duplicate indexes array so duplicates have -1s
    // the rest increment naturally
    // 1  2  3  4  5  6  -1  7  8  -1  9  -1  -1  -1  10  11
    for (int i = 0; i < g->numLines; i++) {
        if (duplicateLineIndexes[i] != -1) {
            duplicateLineIndexes[i] = -1;
            // by how many indexes is the array currently shifting
            // back to cover up the holes of the duplicated indexes
            indexLineOffset++;
        }
        else {
            duplicateLineIndexes[i] = i - indexLineOffset;
        }
    }
    unsigned short *newLines = malloc(sizeof(unsigned short) * newNumLines * 2);
    for (int i = 0; i < g->numLines; i++) {
        if (duplicateLineIndexes[i] != -1) {
            newLines[duplicateLineIndexes[i] * 2 + 0] = g->lines[i * 2 + 0];
            newLines[duplicateLineIndexes[i] * 2 + 1] = g->lines[i * 2 + 1];
        }
    }
    g->numLines = 0;
    free(g->lines);
    g->lines = newLines;
    g->numLines = newNumLines;
    //
    //   END DUPLICATE LINES
    //

    unsigned int indexPointOffset = 0;
    unsigned int newNumPoints = 0;
    // invert duplicate indexes array so duplicates have -1s
    // the rest are their own indexes, in the new collapsed array,
    // which removes all the duplicated indexes completely. looks like:
    // 1  2  3  4  5  6  -1  7  8  -1  9  -1  -1  -1  10  11
    for (int i = 0; i < g->numPoints; i++) {
        if (duplicateIndexes[i] != -1) {
            duplicateIndexes[i] = -1;
            // by how many indexes is the array currently shifting
            // back to cover up the holes of the duplicated indexes
            indexPointOffset++;
        }
        else {
            duplicateIndexes[i] = i - indexPointOffset;
            newNumPoints++;
        }
    }

    flo_t *newPointsArray = malloc(sizeof(flo_t) * newNumPoints * 3);
    for (int i = 0; i < g->numPoints; i++) {
        if (duplicateIndexes[i] != -1) {
            newPointsArray[duplicateIndexes[i] * 3 + 0] = g->points[i * 3 + 0];
            newPointsArray[duplicateIndexes[i] * 3 + 1] = g->points[i * 3 + 1];
            newPointsArray[duplicateIndexes[i] * 3 + 2] = g->points[i * 3 + 2];
        }
    }

    g->numPoints = 0;
    free(g->points);
    g->points = newPointsArray;
    g->numPoints = newNumPoints;

    // finally, update faces and lines with the moved indexes of the shortened point array
    // FACES
    for (int f = 0; f < g->numFaces * 3; f++) {
        if (duplicateIndexes[g->faces[f]] != -1) {
            g->faces[f] = duplicateIndexes[g->faces[f]];
        }
    }
    // LINES
    for (int l = 0; l < g->numLines * 2; l++) {
        if (duplicateIndexes[g->lines[l]] != -1) {
            g->lines[l] = duplicateIndexes[g->lines[l]];
        }
    }
}
