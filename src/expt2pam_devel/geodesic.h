//
//  geodesic
//  geometry builder
//  Copyright (c) 2013 Robby Kraft
//  MIT open source license
//

#ifndef __geodesic__geometry__
#define __geodesic__geometry__

#define M_PI 3.14159265358979323846

typedef double flo_t;
typedef struct geodesicSphere geodesicSphere;

#define TETRAHEDRON_POINT_COUNT 4
#define TETRAHEDRON_LINE_COUNT 6
#define TETRAHEDRON_FACE_COUNT 4

#define OCTAHEDRON_POINT_COUNT 6
#define OCTAHEDRON_LINE_COUNT 12
#define OCTAHEDRON_FACE_COUNT 8

#define HEXAHEDRON_POINT_COUNT 8
#define HEXAHEDRON_LINE_COUNT 12
#define HEXAHEDRON_FACE_COUNT 6
#define HEXAHEDRON_TRIANGLE_FACE_COUNT 12  // OpenGL ES can only render triangles

#define ICOSAHEDRON_POINT_COUNT 12
#define ICOSAHEDRON_LINE_COUNT 30
#define ICOSAHEDRON_FACE_COUNT 20

#define DODECAHEDRON_POINT_COUNT 20
#define DODECAHEDRON_LINE_COUNT 30
#define DODECAHEDRON_FACE_COUNT 12
#define DODECAHEDRON_TRIANGLE_FACE_COUNT 36

// dihedral angle is the angle between two adjacent faces,
//  circling under and running perpendicular to the edge dividing the faces
#define TETRAHEDRON_DIHEDRAL_ANGLE   70.52877936550930863075400066
#define OCTAHEDRON_DIHEDRAL_ANGLE   109.47122063449069136924599934
#define HEXAHEDRON_DIHEDRAL_ANGLE    90
#define ICOSAHEDRON_DIHEDRAL_ANGLE  138.189685104221401934142083269
#define DODECAHEDRON_DIHEDRAL_ANGLE 116.56505117707798935157219372

// inradius, the inscribed sphere, distance from center to midpoint of a face
#define TETRAHEDRON_INRADIUS    0.333333333333333333333333333333
#define OCTAHEDRON_INRADIUS     0.577350269189625764509148780502
#define HEXAHEDRON_INRADIUS     0.5773502691896257645091487805
#define ICOSAHEDRON_INRADIUS    0.794654472291766122955530928331
#define DODECAHEDRON_INRADIUS   0.794654472291766122955530928327

#define TETRAHEDRON_MIDRADIUS   0.577350269189625764509148780501
#define OCTAHEDRON_MIDRADIUS    0.707106781186547524400844362104
#define HEXAHEDRON_MIDRADIUS    0.816496580927726032732428024901
#define ICOSAHEDRON_MIDRADIUS   0.8506508083520399321815404970630
#define DODECAHEDRON_MIDRADIUS  0.9341723589627156964511186235480

#define TETRAHEDRON_SIDE_LENGTH   1.6329931618554520654648560498
#define OCTAHEDRON_SIDE_LENGTH    1.41421356237309504880168872421
#define HEXAHEDRON_SIDE_LENGTH    1.154700538379251529018297561
#define ICOSAHEDRON_SIDE_LENGTH   1.0514622242382672120513381697
#define DODECAHEDRON_SIDE_LENGTH  0.713644179546179863883939686092

#define TETRAHEDRON_VOLUME    0.513200239279667346230354471554
#define OCTAHEDRON_VOLUME     1.33333333333333333333333333333
#define HEXAHEDRON_VOLUME     1.53960071783900203869106341466
#define ICOSAHEDRON_VOLUME    2.53615071012040952564383822238
#define DODECAHEDRON_VOLUME   2.78516386312262296729255491273

extern const double _dodecahedron_points[DODECAHEDRON_POINT_COUNT * 3];

geodesicSphere icosahedronSphere(unsigned int v);
geodesicSphere octahedronSphere(unsigned int v);
geodesicSphere tetrahedronSphere(unsigned int v);
void deleteGeodesicSphere(geodesicSphere *g);

struct geodesicSphere {

    unsigned int numPoints;
    unsigned int numLines;
    unsigned int numFaces;

    flo_t *points;  // count is numPoints * 3
    unsigned short *lines;   // count is numLines * 2    - indices in *points array
    unsigned short *faces;   // count is numFaces * 3    - indices in *points array

    flo_t *pointNormals;  // count is numPoints * 3
    flo_t *lineNormals;   // count is numLines * 3
    flo_t *faceNormals;   // count is numFaces * 3

    unsigned int frequency;

    // legacy data. get creative!
    flo_t *pointsNotSpherized;
    flo_t *pointsDeltaSpherized;  // difference between spherized point and original
    unsigned short *parentFace; // of the original platonic solid faces, from which did the point originate? size of numPoints
};

typedef struct geodesicAnalysis geodesicAnalysis;

geodesicAnalysis classifyLines(geodesicSphere *g);

struct geodesicAnalysis {
    double *lineLengthValues;
    unsigned int numLineLengths;      // the above array size
    unsigned int *lineLengthTypes;    // count is numLines, pointers to indices in lineLengthValues
    unsigned int *lineTypesQuantities;
};

#endif
