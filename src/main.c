#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "constants.h"
#include "marchingcubes.h"
#include "vofi.h"

double Cx, Cy, Cz;

double impl_func_ellipsoid(double xy[]);
double impl_func_sphere(double xy[]);
double impl_func_sphere_chris(double xy[]);
double impl_func_orthocircle(double xy[]);
double impl_func_sinwave1(double xy[]);
double impl_func_plane(double xy[]);
double impl_func_jet(double xy[]);
double impl_func_dodecahedron(double xy[]);
double impl_func_genus2(double xy[]);

int main(int argc, char *argv[])
{
  clock_t start, diff;
  long int nx, ny, nz, msec, *vertCtr, *vert;
  double dx, dy, dz, *colour, *colour_exact, *normalx, *normaly, *normalz;
  struct CONTOUR *contour, *refmesh;
  contour = (struct CONTOUR *) malloc(sizeof(struct CONTOUR));
  refmesh = (struct CONTOUR *) malloc(sizeof(struct CONTOUR));

  /**************** Choose case ********************/
  /* Between: PLANE -- SPHERE -- ELLIPSOID -- ORTHOCIRCLE -- SINWAVE1 --
   * JET -- DODECAHEDRON -- GENUS2 */
  /*************************************************/

  if (argc < 2)
  {
    printf("Usage: curveval <case> <nx> <interpolation_type>\n");
    printf(
        "Case to choose from PLANE -- SPHERE -- ELLIPSOID -- "
        "ORTHOCIRCLE -- SINWAsVE1 -- JET -- DODECAHEDRON \n");
    printf("Interpolation methods to choose from MASON -- LINEAR -- MIDDLE\n");
    return 1;
  }

  char *cas = argv[1];
  nx = ny = nz = atoi(argv[2]);
  char *interp = argv[3];

  printf(
      "+-----------------------------------------------------------------------"
      "-------------------------+\n");
  printf("+    Curvature evaluation test on %i x %i x %i cells for the %s\n", nx, ny, nz, cas);
  printf(
      "+-----------------------------------------------------------------------"
      "-------------------------+\n");

  printf("+    Initialising colour\n");

  start = clock();
  colour = (double *) malloc(nx * ny * nz * sizeof(double));
  normalx = (double *) malloc(nx * ny * nz * sizeof(double));
  normaly = (double *) malloc(nx * ny * nz * sizeof(double));
  normalz = (double *) malloc(nx * ny * nz * sizeof(double));
  vertCtr = (long int *) malloc(nx * ny * nz * sizeof(long int));
  vert = (long int *) malloc(6 * nx * ny * nz * sizeof(long int));
  InitialiseDomain(nx, ny, nz, &dx, &dy, &dz, cas);
  InitialiseColourField(nx, ny, nz, dx, dy, dz, colour, cas);

  diff = clock() - start;
  msec = diff * 1000 / CLOCKS_PER_SEC;
  printf("+      -> dx = %f %f %f\n", dx, dy, dz);

  printf("+      -> Colour initialised in %d s %d ms\n", msec / 1000, msec % 1000);

  printf("+    Generating contour\n");
  long int interptype = 1;
  if (strncasecmp(interp, "MASON", 5) == 0) interptype = 1;
  if (strncasecmp(interp, "LINEAR", 6) == 0) interptype = 2;
  if (strncasecmp(interp, "MIDDLE", 6) == 0) interptype = 3;
  start = clock();
  GenerateContourMC(nx, ny, nz, dx, dy, dz, colour, contour, vertCtr, vert, interptype);
  diff = clock() - start;
  msec = diff * 1000 / CLOCKS_PER_SEC;

  printf("+      -> %i vertices and %i triangles\n", contour->nVertices, contour->nTriangles);
  printf("+      -> Contour generated in %d s %d ms\n", msec / 1000, msec % 1000);

  printf(
      "+-----------------------------------------------------------------------"
      "-------------------------+\n");
  printf("+    Writing contour\n");

  WriteSurfaceMeshSTL(contour, "surface");
  WriteEulerianMeshVTI(nx, ny, nz, dx, dy, dz, colour, "grid");

  printf(
      "+-----------------------------------------------------------------------"
      "-------------------------+\n");

  return 0;
}

int InitialiseDomain(long int nx, long int ny, long int nz, double *dx, double *dy, double *dz, char *cas)
{
  double lx, ly, lz;

  if ((strncasecmp(cas, "SPHERE", 6) == 0) || (strncasecmp(cas, "ELLIPSOID", 9) == 0) || (strncasecmp(cas, "SINWAVE1", 8) == 0) ||
      (strncasecmp(cas, "PLANE", 5) == 0) || (strncasecmp(cas, "JET", 5) == 0) || (strncasecmp(cas, "DODECAHEDRON", 12) == 0))
  {
    lx = ly = lz = 8.0;
    *dx = lx / ((double) nx);
    *dy = ly / ((double) ny);
    *dz = lz / ((double) nz);

    struct timeval tv;
    gettimeofday(&tv, NULL);

    srand(tv.tv_usec);
    Cx = 0.5 * lx;
    Cy = 0.5 * ly;
    Cz = 0.5 * lz;
    printf("+      -> Center is at %.15e, %.15e, %.15e\n", Cx, Cy, Cz);
  }
  else if (strncasecmp(cas, "ORTHOCIRCLE", 11) == 0)
  {
    lx = ly = lz = 3.0;

    *dx = lx / ((double) nx);
    *dy = ly / ((double) ny);
    *dz = lz / ((double) nz);

    struct timeval tv;
    gettimeofday(&tv, NULL);

    srand(tv.tv_usec);
    Cx = 1.5;  // + (0.5 * (double) rand() / (double) (RAND_MAX)) * *dx;
    Cy = 1.5;  // + (0.5 * (double) rand() / (double) (RAND_MAX)) * *dy;
    Cz = 1.5;  // + (0.5 * (double) rand() / (double) (RAND_MAX)) * *dz;
    printf("+      -> Center is at %e, %e, %e\n", Cx, Cy, Cz);
  }
  else if ((strncasecmp(cas, "GENUS2", 6) == 0))
  {
    lx = ly = lz = 4.0;
    *dx = lx / ((double) nx);
    *dy = ly / ((double) ny);
    *dz = lz / ((double) nz);

    struct timeval tv;
    gettimeofday(&tv, NULL);

    srand(tv.tv_usec);
    Cx = 2.0;  // + (0.5 * (double) rand() / (double) (RAND_MAX)) * *dx;
    Cy = 2.0;  // + (0.5 * (double) rand() / (double) (RAND_MAX)) * *dy;
    Cz = 2.0;  // + (0.5 * (double) rand() / (double) (RAND_MAX)) * *dz;
    printf("+      -> Center is at %e, %e, %e\n", Cx, Cy, Cz);
  }

  return 0;
}

int InitialiseColourField(long int nx, long int ny, long int nz, double dx, double dy, double dz, double *colour, char *cas)
{
  long int i, j, k, itrue, ndim0 = 3;
  double Xloc[3], X0[3], fh;

  X0[0] = 0.0;
  X0[1] = 0.0;
  X0[2] = 0.0;

  if (strncasecmp(cas, "SPHERE", 6) == 0)
    fh = vofi_Get_fh(impl_func_sphere, X0, dx, ndim0, itrue);
  else if (strncasecmp(cas, "ELLIPSOID", 9) == 0)
    fh = vofi_Get_fh(impl_func_ellipsoid, X0, dx, ndim0, itrue);
  else if (strncasecmp(cas, "ORTHOCIRCLE", 11) == 0)
    fh = vofi_Get_fh(impl_func_orthocircle, X0, dx, ndim0, itrue);
  else if (strncasecmp(cas, "SINWAVE1", 8) == 0)
    fh = vofi_Get_fh(impl_func_sinwave1, X0, dx, ndim0, itrue);
  else if (strncasecmp(cas, "PLANE", 5) == 0)
    fh = vofi_Get_fh(impl_func_plane, X0, dx, ndim0, itrue);
  else if (strncasecmp(cas, "JET", 3) == 0)
    fh = vofi_Get_fh(impl_func_jet, X0, dx, ndim0, itrue);
  else if (strncasecmp(cas, "DODECAHEDRON", 12) == 0)
    fh = vofi_Get_fh(impl_func_dodecahedron, X0, dx, ndim0, itrue);
  else if (strncasecmp(cas, "GENUS2", 12) == 0)
    fh = vofi_Get_fh(impl_func_genus2, X0, dx, ndim0, itrue);

  for (i = 0; i < nx; i++)
  {
    for (j = 0; j < ny; j++)
    {
      for (k = 0; k < nz; k++)
      {
        Xloc[0] = X0[0] + ((double) i) * dx;
        Xloc[1] = X0[1] + ((double) j) * dx;
        Xloc[2] = X0[2] + ((double) k) * dx;

        if (strncasecmp(cas, "SPHERE", 6) == 0)
          colour[i + j * nx + k * nx * ny] = vofi_Get_cc(impl_func_sphere, Xloc, dx, fh, ndim0);
        else if (strncasecmp(cas, "ELLIPSOID", 9) == 0)
          colour[i + j * nx + k * nx * ny] = vofi_Get_cc(impl_func_ellipsoid, Xloc, dx, fh, ndim0);
        else if (strncasecmp(cas, "ORTHOCIRCLE", 11) == 0)
          colour[i + j * nx + k * nx * ny] = vofi_Get_cc(impl_func_orthocircle, Xloc, dx, fh, ndim0);
        else if (strncasecmp(cas, "SINWAVE1", 8) == 0)
          colour[i + j * nx + k * nx * ny] = vofi_Get_cc(impl_func_sinwave1, Xloc, dx, fh, ndim0);
        else if (strncasecmp(cas, "PLANE", 8) == 0)
          colour[i + j * nx + k * nx * ny] = vofi_Get_cc(impl_func_plane, Xloc, dx, fh, ndim0);
        else if (strncasecmp(cas, "JET", 3) == 0)
          colour[i + j * nx + k * nx * ny] = vofi_Get_cc(impl_func_jet, Xloc, dx, fh, ndim0);
        else if (strncasecmp(cas, "DODECAHEDRON", 12) == 0)
          colour[i + j * nx + k * nx * ny] = vofi_Get_cc(impl_func_dodecahedron, Xloc, dx, fh, ndim0);
        else if (strncasecmp(cas, "GENUS2", 6) == 0)
          colour[i + j * nx + k * nx * ny] = vofi_Get_cc(impl_func_genus2, Xloc, dx, fh, ndim0);
      }
    }
  }

  return 0;
}

int WriteEulerianMeshVTI(long int nx, long int ny, long int nz, double dx, double dy, double dz, double *colour, char *filename)
{
  long int i, j, k;
  FILE *fp;

  char name[500];
  sprintf(name, "%s.vti", filename);

  fp = fopen(name, "w");
  fprintf(fp, "<?xml version=\"1.0\"?>\n");
  fprintf(fp,
          "<VTKFile type=\"ImageData\" version=\"0.1\" "
          "byte_order=\"LittleEndian\">\n");
  fprintf(fp,
          "<ImageData WholeExtent=\"%i %i %i %i %i %i\" Origin=\"%i %i %i\" "
          "Spacing=\"%f %f %f\">\n",
          0, nx, 0, ny, 0, nz, 0, 0, 0, dx, dy, dz);
  fprintf(fp, "<Piece Extent=\"%i %i %i %i %i %i\">\n", 0, nx, 0, ny, 0, nz);
  fprintf(fp, "<CellData>\n");
  fprintf(fp, "<DataArray type=\"Float64\" Name=\"Colour\" format=\"ascii\">\n");
  for (k = 0; k < nz; k++)
    for (j = 0; j < ny; j++)
      for (i = 0; i < nx; i++) fprintf(fp, "%lf\n", colour[i + j * nx + k * nx * ny]);
  fprintf(fp, "</DataArray>\n");
  fprintf(fp, "</CellData>\n");
  fprintf(fp, "</Piece>\n");
  fprintf(fp, "</ImageData>\n");
  fprintf(fp, "</VTKFile>\n");
  fclose(fp);

  return 0;
}

int WriteSurfaceMeshSTL(struct CONTOUR *contour, char *filename)
{
  long int i, j, k;
  FILE *fp;

  /* Output */
  long int outVert, outTri;
  outVert = outTri = 0;

  for (i = 0; i < contour->nTriangles; i++)
    if (contour->triangles[i].cas >= 0) outTri++;

  char header[80] = "MF Extracted Contour";
  char bytecount[2] = "00";
  long int nfacets = outTri;
  float tmpvertexes[9], dummynormal[3];

  char name[500];
  sprintf(name, "%s.stl", filename);
  fp = fopen(name, "wb");
  fwrite(header, sizeof(char), 80, fp);
  fwrite(&nfacets, sizeof(int), 1, fp);

  for (i = 0; i < contour->nTriangles; i++)
  {
    if (contour->triangles[i].cas >= 0)
    {
      fwrite(dummynormal, sizeof(float), 3, fp);
      for (j = 0; j < 3; j++)
      {
        tmpvertexes[3 * j + 0] = (float) contour->vertices[contour->triangles[i].p[j]].x;
        tmpvertexes[3 * j + 1] = (float) contour->vertices[contour->triangles[i].p[j]].y;
        tmpvertexes[3 * j + 2] = (float) contour->vertices[contour->triangles[i].p[j]].z;
      }
      fwrite(tmpvertexes, sizeof(float), 9, fp);
      fwrite(bytecount, sizeof(char), 2, fp);
    }
  }
  fclose(fp);

  return 0;
}

double impl_func_orthocircle(double xy[])
{
  double x, y, z, f0;

  double cx = Cx;
  double cy = Cy;
  double cz = Cz;

  double c1 = 0.075;
  double c2 = 3.0;

  x = xy[0];
  y = xy[1];
  z = xy[2];

  f0 = (MF_SQ(MF_SQ(x - cx) + MF_SQ(y - cy) - 1.0) + MF_SQ(z - cz)) * (MF_SQ(MF_SQ(y - cy) + MF_SQ(z - cz) - 1.0) + MF_SQ(x - cx)) *
           (MF_SQ(MF_SQ(z - cz) + MF_SQ(x - cx) - 1.0) + MF_SQ(y - cy)) -
       MF_SQ(c1) * (1.0 + c2 * (MF_SQ(x - cx) + MF_SQ(y - cy) + MF_SQ(z - cz)));

  return f0;
}

double impl_func_sphere(double xy[])
{
  double x, y, z, f0;

  double r = 2.0;

  x = xy[0];
  y = xy[1];
  z = xy[2];

  f0 = MF_SQ(x - Cx) + MF_SQ(y - Cy) + MF_SQ(z - Cz) - r * r;

  return f0;
}

double impl_func_ellipsoid(double xy[])
{
  double x, y, z, f0;

  double cx = Cx;
  double cy = Cy;
  double cz = Cz;

  double r = 2.0;

  x = xy[0];
  y = xy[1];
  z = xy[2];

  f0 = MF_SQ(x - cx) / MF_SQ(2.0) + MF_SQ(y - cy) / MF_SQ(1.5) + MF_SQ(z - cz) / MF_SQ(1.0) - 1.0;

  return f0;
}

double impl_func_sinwave1(double xy[])
{
  double x, y, z, f0;

  double a = 0.2;
  double b = Cz;

  x = xy[0];
  y = xy[1];
  z = xy[2];

  f0 = a * x * sin(x * 0.5 * MF_PI) * sin(y * 0.5 * MF_PI) + b - z;

  return f0;
}

double impl_func_plane(double xy[])
{
  double x, y, z, f0;

  double a = 1.0;
  double b = 1.0;
  double c = 1.0;

  x = xy[0];
  y = xy[1];
  z = xy[2];

  f0 = a * x + b * y + c * z - (a * 4.0 + b * 4.0 + c * 4.0);

  return f0;
}

double impl_func_jet(double xy[])
{
  double x, y, z, f0;

  double cy = Cy;
  double cz = Cz;
  double ro = 1.0;
  double epso = 0.3;

  x = xy[0];
  y = xy[1];
  z = xy[2];

  f0 = MF_SQ(z - cz) + MF_SQ(y - cy) - MF_SQ(ro * (1.0 + epso * sin(x * MF_PI)));

  return f0;
}

double impl_func_dodecahedron(double xy[])
{
  double x, y, z, f0;
  double cx = Cx;
  double cy = Cy;
  double cz = Cz;
  double ro = 1.0;
  double epso = 0.3;

  x = xy[0];
  y = xy[1];
  z = xy[2];

  f0 = MF_POW(x - cx, 6.0) + MF_POW(y - cy, 6.0) + MF_POW(z - cz, 6.0) +
       20.0 * (MF_POW(x - cx, 4.0) * MF_SQ(y - cy) + MF_POW(y - cy, 4.0) * MF_SQ(x - cx) + MF_POW(z - cz, 4.0) * MF_SQ(x - cx) - 2.0);

  return f0;
}

double impl_func_genus2(double xy[])
{
  double x, y, z, f0;
  double cx = Cx;
  double cy = Cy;
  double cz = Cz;

  x = xy[0];
  y = xy[1];
  z = xy[2];

  f0 = 2.0 * (y - Cy) * (MF_SQ(y - Cy) - 3.0 * MF_SQ(x - Cx)) * (1.0 - MF_SQ(z - Cz)) + MF_SQ(MF_SQ(x - Cx) + MF_SQ(y - Cy)) -
       (9.0 * MF_SQ(z - Cz) - 1.0) * (1.0 - MF_SQ(z - Cz));

  //  f0 = -3.0 + 1.0 * (MF_POW(x - cx, 4.0) + MF_POW(y - cy, 4.0) + MF_POW(z -
  //  cz, 4.0))
  //      - 0.0 * (MF_SQ(x - cx) + MF_SQ(y - cy) + MF_SQ(z - cz));

  return f0;
}
