#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "marchingcubes.h"

/* This code is adapted from Paul Bourke's marching-cube code, available at:
    https://paulbourke.net/geometry/polygonise/
*/

int GenerateContourMC(long int nx, long int ny, long int nz, double dx, double dy, double dz, double *colour, struct CONTOUR *contour, long int *vertCtr,
                      long int *vert, long int interptype)
{
  long int i, j, k, n, m;
  long int edgeindex, cubeindex, cellindex, vertexindex1, vertexindex2, vertexindex3, check1, check2;
  long int vertlist[12], *edgx, *edgy, *edgz;
  double isolevel;
  struct GRIDCELL grid;
  struct STENCIL stencil;

  isolevel = 0.5;

  edgx = (long int *) malloc((nx - 1) * ny * nz * sizeof(long int));
  edgy = (long int *) malloc(nx * (ny - 1) * nz * sizeof(long int));
  edgz = (long int *) malloc(nx * ny * (nz - 1) * sizeof(long int));

  for (i = 0; i < (nx - 1) * ny * nz; i++) edgx[i] = -1;

  for (i = 0; i < nx * (ny - 1) * nz; i++) edgy[i] = -1;

  for (i = 0; i < nx * ny * (nz - 1); i++) edgz[i] = -1;

  for (i = 0; i < nx * ny * nz; i++) vertCtr[i] = 0;

  contour->nTriangles = 0;
  contour->nVertices = 0;

  for (i = 0; i < nx - 1; i++)
  {
    for (j = 0; j < ny - 1; j++)
    {
      for (k = 0; k < nz - 1; k++)
      {
        cellindex = i + j * nx + k * nx * ny;

        /* Determine the index into the edge table which
         tells us which contour->vertices are inside of the surface */

        grid.val[0] = colour[i + j * nx + k * nx * ny];
        grid.val[1] = colour[i + 1 + j * nx + k * nx * ny];
        grid.val[2] = colour[i + 1 + j * nx + (k + 1) * nx * ny];
        grid.val[3] = colour[i + j * nx + (k + 1) * nx * ny];
        grid.val[4] = colour[i + (j + 1) * nx + k * nx * ny];
        grid.val[5] = colour[i + 1 + (j + 1) * nx + k * nx * ny];
        grid.val[6] = colour[i + 1 + (j + 1) * nx + (k + 1) * nx * ny];
        grid.val[7] = colour[i + (j + 1) * nx + (k + 1) * nx * ny];

        grid.p[0].x = grid.p[3].x = grid.p[4].x = grid.p[7].x = ((double) i) * dx + 0.5 * dx;
        grid.p[1].x = grid.p[2].x = grid.p[5].x = grid.p[6].x = ((double) i + 1) * dx + 0.5 * dx;
        grid.p[0].y = grid.p[1].y = grid.p[2].y = grid.p[3].y = ((double) j) * dy + 0.5 * dy;
        grid.p[4].y = grid.p[5].y = grid.p[6].y = grid.p[7].y = ((double) j + 1) * dy + 0.5 * dy;
        grid.p[0].z = grid.p[1].z = grid.p[4].z = grid.p[5].z = ((double) k) * dz + 0.5 * dz;
        grid.p[2].z = grid.p[3].z = grid.p[6].z = grid.p[7].z = ((double) k + 1) * dz + 0.5 * dz;

        cubeindex = 0;
        if (grid.val[0] < isolevel) cubeindex |= 1;
        if (grid.val[1] < isolevel) cubeindex |= 2;
        if (grid.val[2] < isolevel) cubeindex |= 4;
        if (grid.val[3] < isolevel) cubeindex |= 8;
        if (grid.val[4] < isolevel) cubeindex |= 16;
        if (grid.val[5] < isolevel) cubeindex |= 32;
        if (grid.val[6] < isolevel) cubeindex |= 64;
        if (grid.val[7] < isolevel) cubeindex |= 128;

        /* Cube is not entirely in/out of the surface */
        if (edgeTable[cubeindex] != 0)
        {
          /* Find the contour->vertices where the surface intersects the cube
           * if vertice has not been generated before, then generates it */
          if (edgeTable[cubeindex] & 1)
          {
            /* Edge 0 follows axis X with coordinates {i, j, k} */
            edgeindex = i + j * (nx - 1) + k * (nx - 1) * ny;
            if (edgx[edgeindex] == -1)
            {
              DynamicallyVerticeGrow(&contour->vertices, contour->nVertices, 1);
              if (interptype < 2)
                contour->vertices[contour->nVertices] = VertexInterpManson(isolevel, grid.p[0], grid.p[1], grid.val[0], grid.val[1]);
              else if (interptype == 2)
                contour->vertices[contour->nVertices] = VertexInterp(isolevel, grid.p[0], grid.p[1], grid.val[0], grid.val[1]);
              else
                contour->vertices[contour->nVertices] = VertexInterpHalf(isolevel, grid.p[0], grid.p[1], grid.val[0], grid.val[1]);

              if (contour->vertices[contour->nVertices].mu < 0.5)
                contour->vertices[contour->nVertices].parent = i + j * nx + k * nx * ny;
              else
                contour->vertices[contour->nVertices].parent = i + 1 + j * nx + k * nx * ny;

              contour->vertices[contour->nVertices].p[0] = i + j * nx + k * nx * ny;
              contour->vertices[contour->nVertices].p[1] = i + 1 + j * nx + k * nx * ny;

              vert[6 * contour->vertices[contour->nVertices].p[0] + vertCtr[contour->vertices[contour->nVertices].p[0]]] = contour->nVertices;
              vert[6 * contour->vertices[contour->nVertices].p[1] + vertCtr[contour->vertices[contour->nVertices].p[1]]] = contour->nVertices;
              vertCtr[contour->vertices[contour->nVertices].p[0]] += 1;
              vertCtr[contour->vertices[contour->nVertices].p[1]] += 1;

              edgx[edgeindex] = contour->nVertices;
              contour->vertices[contour->nVertices].dir = 0;
              vertlist[0] = contour->nVertices;
              contour->nVertices += 1;
            }
            else
            {
              vertlist[0] = edgx[edgeindex];
            }
          }
          if (edgeTable[cubeindex] & 2)
          {
            /* Edge 1 follows axis Z with coordinates {i + 1, j, k} */
            edgeindex = i + 1 + j * nx + k * nx * ny;
            if (edgz[edgeindex] == -1)
            {
              DynamicallyVerticeGrow(&contour->vertices, contour->nVertices, 1);
              if (interptype < 2)
                contour->vertices[contour->nVertices] = VertexInterpManson(isolevel, grid.p[1], grid.p[2], grid.val[1], grid.val[2]);
              else if (interptype == 2)
                contour->vertices[contour->nVertices] = VertexInterp(isolevel, grid.p[1], grid.p[2], grid.val[1], grid.val[2]);
              else
                contour->vertices[contour->nVertices] = VertexInterpHalf(isolevel, grid.p[1], grid.p[2], grid.val[1], grid.val[2]);
              if (contour->vertices[contour->nVertices].mu < 0.5)
                contour->vertices[contour->nVertices].parent = i + 1 + j * nx + k * nx * ny;
              else
                contour->vertices[contour->nVertices].parent = i + 1 + j * nx + (k + 1) * nx * ny;

              contour->vertices[contour->nVertices].p[0] = i + 1 + j * nx + k * nx * ny;
              contour->vertices[contour->nVertices].p[1] = i + 1 + j * nx + (k + 1) * nx * ny;

              vert[6 * contour->vertices[contour->nVertices].p[0] + vertCtr[contour->vertices[contour->nVertices].p[0]]] = contour->nVertices;
              vert[6 * contour->vertices[contour->nVertices].p[1] + vertCtr[contour->vertices[contour->nVertices].p[1]]] = contour->nVertices;
              vertCtr[contour->vertices[contour->nVertices].p[0]] += 1;
              vertCtr[contour->vertices[contour->nVertices].p[1]] += 1;

              edgz[edgeindex] = contour->nVertices;
              contour->vertices[contour->nVertices].dir = 2;
              vertlist[1] = contour->nVertices;
              contour->nVertices += 1;
            }
            else
            {
              vertlist[1] = edgz[edgeindex];
            }
          }
          if (edgeTable[cubeindex] & 4)
          {
            /* Edge 2 follows axis X with coordinates {i, j, k + 1} */
            edgeindex = i + j * (nx - 1) + (k + 1) * (nx - 1) * ny;
            if (edgx[edgeindex] == -1)
            {
              DynamicallyVerticeGrow(&contour->vertices, contour->nVertices, 1);
              if (interptype < 2)
                contour->vertices[contour->nVertices] = VertexInterpManson(isolevel, grid.p[3], grid.p[2], grid.val[3], grid.val[2]);
              else if (interptype == 2)
                contour->vertices[contour->nVertices] = VertexInterp(isolevel, grid.p[3], grid.p[2], grid.val[3], grid.val[2]);
              else
                contour->vertices[contour->nVertices] = VertexInterpHalf(isolevel, grid.p[3], grid.p[2], grid.val[3], grid.val[2]);
              if (contour->vertices[contour->nVertices].mu < 0.5)
                contour->vertices[contour->nVertices].parent = i + j * nx + (k + 1) * nx * ny;
              else
                contour->vertices[contour->nVertices].parent = i + 1 + j * nx + (k + 1) * nx * ny;

              contour->vertices[contour->nVertices].p[0] = i + j * nx + (k + 1) * nx * ny;
              contour->vertices[contour->nVertices].p[1] = i + 1 + j * nx + (k + 1) * nx * ny;

              vert[6 * contour->vertices[contour->nVertices].p[0] + vertCtr[contour->vertices[contour->nVertices].p[0]]] = contour->nVertices;
              vert[6 * contour->vertices[contour->nVertices].p[1] + vertCtr[contour->vertices[contour->nVertices].p[1]]] = contour->nVertices;
              vertCtr[contour->vertices[contour->nVertices].p[0]] += 1;
              vertCtr[contour->vertices[contour->nVertices].p[1]] += 1;

              edgx[edgeindex] = contour->nVertices;
              contour->vertices[contour->nVertices].dir = 0;
              vertlist[2] = contour->nVertices;
              contour->nVertices += 1;
            }
            else
            {
              vertlist[2] = edgx[edgeindex];
            }
          }
          if (edgeTable[cubeindex] & 8)
          {
            /* Edge 3 follows axis Z with coordinates {i, j, k} */
            edgeindex = i + j * nx + k * nx * ny;
            if (edgz[edgeindex] == -1)
            {
              DynamicallyVerticeGrow(&contour->vertices, contour->nVertices, 1);
              if (interptype < 2)
                contour->vertices[contour->nVertices] = VertexInterpManson(isolevel, grid.p[0], grid.p[3], grid.val[0], grid.val[3]);
              else if (interptype == 2)
                contour->vertices[contour->nVertices] = VertexInterp(isolevel, grid.p[0], grid.p[3], grid.val[0], grid.val[3]);
              else
                contour->vertices[contour->nVertices] = VertexInterpHalf(isolevel, grid.p[0], grid.p[3], grid.val[0], grid.val[3]);
              if (contour->vertices[contour->nVertices].mu < 0.5)
                contour->vertices[contour->nVertices].parent = i + j * nx + k * nx * ny;
              else
                contour->vertices[contour->nVertices].parent = i + j * nx + (k + 1) * nx * ny;

              contour->vertices[contour->nVertices].p[0] = i + j * nx + k * nx * ny;
              contour->vertices[contour->nVertices].p[1] = i + j * nx + (k + 1) * nx * ny;

              vert[6 * contour->vertices[contour->nVertices].p[0] + vertCtr[contour->vertices[contour->nVertices].p[0]]] = contour->nVertices;
              vert[6 * contour->vertices[contour->nVertices].p[1] + vertCtr[contour->vertices[contour->nVertices].p[1]]] = contour->nVertices;
              vertCtr[contour->vertices[contour->nVertices].p[0]] += 1;
              vertCtr[contour->vertices[contour->nVertices].p[1]] += 1;

              edgz[edgeindex] = contour->nVertices;
              contour->vertices[contour->nVertices].dir = 2;
              vertlist[3] = contour->nVertices;
              contour->nVertices += 1;
            }
            else
            {
              vertlist[3] = edgz[edgeindex];
            }
          }
          if (edgeTable[cubeindex] & 16)
          {
            /* Edge 4 follows axis X with coordinates {i, j + 1, k} */
            edgeindex = i + (j + 1) * (nx - 1) + k * (nx - 1) * ny;
            if (edgx[edgeindex] == -1)
            {
              DynamicallyVerticeGrow(&contour->vertices, contour->nVertices, 1);
              if (interptype < 2)
                contour->vertices[contour->nVertices] = VertexInterpManson(isolevel, grid.p[4], grid.p[5], grid.val[4], grid.val[5]);
              else if (interptype == 2)
                contour->vertices[contour->nVertices] = VertexInterp(isolevel, grid.p[4], grid.p[5], grid.val[4], grid.val[5]);
              else
                contour->vertices[contour->nVertices] = VertexInterpHalf(isolevel, grid.p[4], grid.p[5], grid.val[4], grid.val[5]);
              if (contour->vertices[contour->nVertices].mu < 0.5)
                contour->vertices[contour->nVertices].parent = i + (j + 1) * nx + k * nx * ny;
              else
                contour->vertices[contour->nVertices].parent = i + 1 + (j + 1) * nx + k * nx * ny;

              contour->vertices[contour->nVertices].p[0] = i + (j + 1) * nx + k * nx * ny;
              contour->vertices[contour->nVertices].p[1] = i + 1 + (j + 1) * nx + k * nx * ny;

              vert[6 * contour->vertices[contour->nVertices].p[0] + vertCtr[contour->vertices[contour->nVertices].p[0]]] = contour->nVertices;
              vert[6 * contour->vertices[contour->nVertices].p[1] + vertCtr[contour->vertices[contour->nVertices].p[1]]] = contour->nVertices;
              vertCtr[contour->vertices[contour->nVertices].p[0]] += 1;
              vertCtr[contour->vertices[contour->nVertices].p[1]] += 1;

              edgx[edgeindex] = contour->nVertices;
              contour->vertices[contour->nVertices].dir = 0;
              vertlist[4] = contour->nVertices;
              contour->nVertices += 1;
            }
            else
            {
              vertlist[4] = edgx[edgeindex];
            }
          }
          if (edgeTable[cubeindex] & 32)
          {
            /* Edge 5 follows axis Z with coordinates {i + 1, j + 1, k} */
            edgeindex = i + 1 + (j + 1) * nx + k * nx * ny;
            if (edgz[edgeindex] == -1)
            {
              DynamicallyVerticeGrow(&contour->vertices, contour->nVertices, 1);
              if (interptype < 2)
                contour->vertices[contour->nVertices] = VertexInterpManson(isolevel, grid.p[5], grid.p[6], grid.val[5], grid.val[6]);
              else if (interptype == 2)
                contour->vertices[contour->nVertices] = VertexInterp(isolevel, grid.p[5], grid.p[6], grid.val[5], grid.val[6]);
              else
                contour->vertices[contour->nVertices] = VertexInterpHalf(isolevel, grid.p[5], grid.p[6], grid.val[5], grid.val[6]);
              if (contour->vertices[contour->nVertices].mu < 0.5)
                contour->vertices[contour->nVertices].parent = i + 1 + (j + 1) * nx + k * nx * ny;
              else
                contour->vertices[contour->nVertices].parent = i + 1 + (j + 1) * nx + (k + 1) * nx * ny;

              contour->vertices[contour->nVertices].p[0] = i + 1 + (j + 1) * nx + k * nx * ny;
              contour->vertices[contour->nVertices].p[1] = i + 1 + (j + 1) * nx + (k + 1) * nx * ny;

              vert[6 * contour->vertices[contour->nVertices].p[0] + vertCtr[contour->vertices[contour->nVertices].p[0]]] = contour->nVertices;
              vert[6 * contour->vertices[contour->nVertices].p[1] + vertCtr[contour->vertices[contour->nVertices].p[1]]] = contour->nVertices;
              vertCtr[contour->vertices[contour->nVertices].p[0]] += 1;
              vertCtr[contour->vertices[contour->nVertices].p[1]] += 1;

              edgz[edgeindex] = contour->nVertices;
              contour->vertices[contour->nVertices].dir = 2;
              vertlist[5] = contour->nVertices;
              contour->nVertices += 1;
            }
            else
            {
              vertlist[5] = edgz[edgeindex];
            }
          }
          if (edgeTable[cubeindex] & 64)
          {
            /* Edge 6 follows axis X with coordinates {i, j + 1, k + 1} */
            edgeindex = i + (j + 1) * (nx - 1) + (k + 1) * (nx - 1) * ny;
            if (edgx[edgeindex] == -1)
            {
              DynamicallyVerticeGrow(&contour->vertices, contour->nVertices, 1);
              if (interptype < 2)
                contour->vertices[contour->nVertices] = VertexInterpManson(isolevel, grid.p[7], grid.p[6], grid.val[7], grid.val[6]);
              else if (interptype == 2)
                contour->vertices[contour->nVertices] = VertexInterp(isolevel, grid.p[7], grid.p[6], grid.val[7], grid.val[6]);
              else
                contour->vertices[contour->nVertices] = VertexInterpHalf(isolevel, grid.p[7], grid.p[6], grid.val[7], grid.val[6]);
              if (contour->vertices[contour->nVertices].mu < 0.5)
                contour->vertices[contour->nVertices].parent = i + (j + 1) * nx + (k + 1) * nx * ny;
              else
                contour->vertices[contour->nVertices].parent = i + 1 + (j + 1) * nx + (k + 1) * nx * ny;

              contour->vertices[contour->nVertices].p[0] = i + (j + 1) * nx + (k + 1) * nx * ny;
              contour->vertices[contour->nVertices].p[1] = i + 1 + (j + 1) * nx + (k + 1) * nx * ny;

              vert[6 * contour->vertices[contour->nVertices].p[0] + vertCtr[contour->vertices[contour->nVertices].p[0]]] = contour->nVertices;
              vert[6 * contour->vertices[contour->nVertices].p[1] + vertCtr[contour->vertices[contour->nVertices].p[1]]] = contour->nVertices;
              vertCtr[contour->vertices[contour->nVertices].p[0]] += 1;
              vertCtr[contour->vertices[contour->nVertices].p[1]] += 1;

              edgx[edgeindex] = contour->nVertices;
              contour->vertices[contour->nVertices].dir = 0;
              vertlist[6] = contour->nVertices;
              contour->nVertices += 1;
            }
            else
            {
              vertlist[6] = edgx[edgeindex];
            }
          }
          if (edgeTable[cubeindex] & 128)
          {
            /* Edge 7 follows axis Z with coordinates {i, j + 1, k} */
            edgeindex = i + (j + 1) * nx + k * nx * ny;
            if (edgz[edgeindex] == -1)
            {
              DynamicallyVerticeGrow(&contour->vertices, contour->nVertices, 1);
              if (interptype < 2)
                contour->vertices[contour->nVertices] = VertexInterpManson(isolevel, grid.p[4], grid.p[7], grid.val[4], grid.val[7]);
              else if (interptype == 2)
                contour->vertices[contour->nVertices] = VertexInterp(isolevel, grid.p[4], grid.p[7], grid.val[4], grid.val[7]);
              else
                contour->vertices[contour->nVertices] = VertexInterpHalf(isolevel, grid.p[4], grid.p[7], grid.val[4], grid.val[7]);
              if (contour->vertices[contour->nVertices].mu < 0.5)
                contour->vertices[contour->nVertices].parent = i + (j + 1) * nx + k * nx * ny;
              else
                contour->vertices[contour->nVertices].parent = i + (j + 1) * nx + (k + 1) * nx * ny;

              contour->vertices[contour->nVertices].p[0] = i + (j + 1) * nx + k * nx * ny;
              contour->vertices[contour->nVertices].p[1] = i + (j + 1) * nx + (k + 1) * nx * ny;

              vert[6 * contour->vertices[contour->nVertices].p[0] + vertCtr[contour->vertices[contour->nVertices].p[0]]] = contour->nVertices;
              vert[6 * contour->vertices[contour->nVertices].p[1] + vertCtr[contour->vertices[contour->nVertices].p[1]]] = contour->nVertices;
              vertCtr[contour->vertices[contour->nVertices].p[0]] += 1;
              vertCtr[contour->vertices[contour->nVertices].p[1]] += 1;

              edgz[edgeindex] = contour->nVertices;
              contour->vertices[contour->nVertices].dir = 2;
              vertlist[7] = contour->nVertices;
              contour->nVertices += 1;
            }
            else
            {
              vertlist[7] = edgz[edgeindex];
            }
          }
          if (edgeTable[cubeindex] & 256)
          {
            /* Edge 8 follows axis Y with coordinates {i, j, k} */
            edgeindex = i + j * nx + k * nx * (ny - 1);
            if (edgy[edgeindex] == -1)
            {
              DynamicallyVerticeGrow(&contour->vertices, contour->nVertices, 1);
              if (interptype < 2)
                contour->vertices[contour->nVertices] = VertexInterpManson(isolevel, grid.p[0], grid.p[4], grid.val[0], grid.val[4]);
              else if (interptype == 2)
                contour->vertices[contour->nVertices] = VertexInterp(isolevel, grid.p[0], grid.p[4], grid.val[0], grid.val[4]);
              else
                contour->vertices[contour->nVertices] = VertexInterpHalf(isolevel, grid.p[0], grid.p[4], grid.val[0], grid.val[4]);
              if (contour->vertices[contour->nVertices].mu < 0.5)
                contour->vertices[contour->nVertices].parent = i + j * nx + k * nx * ny;
              else
                contour->vertices[contour->nVertices].parent = i + (j + 1) * nx + k * nx * ny;

              contour->vertices[contour->nVertices].p[0] = i + j * nx + k * nx * ny;
              contour->vertices[contour->nVertices].p[1] = i + (j + 1) * nx + k * nx * ny;

              vert[6 * contour->vertices[contour->nVertices].p[0] + vertCtr[contour->vertices[contour->nVertices].p[0]]] = contour->nVertices;
              vert[6 * contour->vertices[contour->nVertices].p[1] + vertCtr[contour->vertices[contour->nVertices].p[1]]] = contour->nVertices;
              vertCtr[contour->vertices[contour->nVertices].p[0]] += 1;
              vertCtr[contour->vertices[contour->nVertices].p[1]] += 1;

              edgy[edgeindex] = contour->nVertices;
              contour->vertices[contour->nVertices].dir = 1;
              vertlist[8] = contour->nVertices;
              contour->nVertices += 1;
            }
            else
            {
              vertlist[8] = edgy[edgeindex];
            }
          }
          if (edgeTable[cubeindex] & 512)
          {
            /* Edge 9 follows axis Y with coordinates {i + 1, j, k} */
            edgeindex = i + 1 + j * nx + k * nx * (ny - 1);
            if (edgy[edgeindex] == -1)
            {
              DynamicallyVerticeGrow(&contour->vertices, contour->nVertices, 1);
              if (interptype < 2)
                contour->vertices[contour->nVertices] = VertexInterpManson(isolevel, grid.p[1], grid.p[5], grid.val[1], grid.val[5]);
              else if (interptype == 2)
                contour->vertices[contour->nVertices] = VertexInterp(isolevel, grid.p[1], grid.p[5], grid.val[1], grid.val[5]);
              else
                contour->vertices[contour->nVertices] = VertexInterpHalf(isolevel, grid.p[1], grid.p[5], grid.val[1], grid.val[5]);
              if (contour->vertices[contour->nVertices].mu < 0.5)
                contour->vertices[contour->nVertices].parent = i + 1 + j * nx + k * nx * ny;
              else
                contour->vertices[contour->nVertices].parent = i + 1 + (j + 1) * nx + k * nx * ny;

              contour->vertices[contour->nVertices].p[0] = i + 1 + j * nx + k * nx * ny;
              contour->vertices[contour->nVertices].p[1] = i + 1 + (j + 1) * nx + k * nx * ny;

              vert[6 * contour->vertices[contour->nVertices].p[0] + vertCtr[contour->vertices[contour->nVertices].p[0]]] = contour->nVertices;
              vert[6 * contour->vertices[contour->nVertices].p[1] + vertCtr[contour->vertices[contour->nVertices].p[1]]] = contour->nVertices;
              vertCtr[contour->vertices[contour->nVertices].p[0]] += 1;
              vertCtr[contour->vertices[contour->nVertices].p[1]] += 1;

              edgy[edgeindex] = contour->nVertices;
              contour->vertices[contour->nVertices].dir = 1;
              vertlist[9] = contour->nVertices;
              contour->nVertices += 1;
            }
            else
            {
              vertlist[9] = edgy[edgeindex];
            }
          }
          if (edgeTable[cubeindex] & 1024)
          {
            /* Edge 10 follows axis Y with coordinates {i + 1, j, k + 1} */
            edgeindex = i + 1 + j * nx + (k + 1) * nx * (ny - 1);
            if (edgy[edgeindex] == -1)
            {
              DynamicallyVerticeGrow(&contour->vertices, contour->nVertices, 1);
              if (interptype < 2)
                contour->vertices[contour->nVertices] = VertexInterpManson(isolevel, grid.p[2], grid.p[6], grid.val[2], grid.val[6]);
              else if (interptype == 2)
                contour->vertices[contour->nVertices] = VertexInterp(isolevel, grid.p[2], grid.p[6], grid.val[2], grid.val[6]);
              else
                contour->vertices[contour->nVertices] = VertexInterpHalf(isolevel, grid.p[2], grid.p[6], grid.val[2], grid.val[6]);
              if (contour->vertices[contour->nVertices].mu < 0.5)
                contour->vertices[contour->nVertices].parent = i + 1 + j * nx + (k + 1) * nx * ny;
              else
                contour->vertices[contour->nVertices].parent = i + 1 + (j + 1) * nx + (k + 1) * nx * ny;

              contour->vertices[contour->nVertices].p[0] = i + 1 + j * nx + (k + 1) * nx * ny;
              contour->vertices[contour->nVertices].p[1] = i + 1 + (j + 1) * nx + (k + 1) * nx * ny;

              vert[6 * contour->vertices[contour->nVertices].p[0] + vertCtr[contour->vertices[contour->nVertices].p[0]]] = contour->nVertices;
              vert[6 * contour->vertices[contour->nVertices].p[1] + vertCtr[contour->vertices[contour->nVertices].p[1]]] = contour->nVertices;
              vertCtr[contour->vertices[contour->nVertices].p[0]] += 1;
              vertCtr[contour->vertices[contour->nVertices].p[1]] += 1;

              edgy[edgeindex] = contour->nVertices;
              contour->vertices[contour->nVertices].dir = 1;
              vertlist[10] = contour->nVertices;
              contour->nVertices += 1;
            }
            else
            {
              vertlist[10] = edgy[edgeindex];
            }
          }
          if (edgeTable[cubeindex] & 2048)
          {
            /* Edge 11 follows axis Y with coordinates {i, j, k + 1} */
            edgeindex = i + j * nx + (k + 1) * nx * (ny - 1);
            if (edgy[edgeindex] == -1)
            {
              DynamicallyVerticeGrow(&contour->vertices, contour->nVertices, 1);
              if (interptype < 2)
                contour->vertices[contour->nVertices] = VertexInterpManson(isolevel, grid.p[3], grid.p[7], grid.val[3], grid.val[7]);
              else if (interptype == 2)
                contour->vertices[contour->nVertices] = VertexInterp(isolevel, grid.p[3], grid.p[7], grid.val[3], grid.val[7]);
              else
                contour->vertices[contour->nVertices] = VertexInterpHalf(isolevel, grid.p[3], grid.p[7], grid.val[3], grid.val[7]);
              if (contour->vertices[contour->nVertices].mu < 0.5)
                contour->vertices[contour->nVertices].parent = i + j * nx + (k + 1) * nx * ny;
              else
                contour->vertices[contour->nVertices].parent = i + (j + 1) * nx + (k + 1) * nx * ny;

              contour->vertices[contour->nVertices].p[0] = i + j * nx + (k + 1) * nx * ny;
              contour->vertices[contour->nVertices].p[1] = i + (j + 1) * nx + (k + 1) * nx * ny;

              vert[6 * contour->vertices[contour->nVertices].p[0] + vertCtr[contour->vertices[contour->nVertices].p[0]]] = contour->nVertices;
              vert[6 * contour->vertices[contour->nVertices].p[1] + vertCtr[contour->vertices[contour->nVertices].p[1]]] = contour->nVertices;
              vertCtr[contour->vertices[contour->nVertices].p[0]] += 1;
              vertCtr[contour->vertices[contour->nVertices].p[1]] += 1;

              edgy[edgeindex] = contour->nVertices;
              contour->vertices[contour->nVertices].dir = 1;
              vertlist[11] = contour->nVertices;
              contour->nVertices += 1;
            }
            else
            {
              vertlist[11] = edgy[edgeindex];
            }
          }

          /* Create the triangle */
          for (n = 0; triTable[cubeindex][n] != -1; n += 3)
          {
            /* Update vertex connectivity */
            vertexindex1 = vertlist[triTable[cubeindex][n]];
            vertexindex2 = vertlist[triTable[cubeindex][n + 1]];
            vertexindex3 = vertlist[triTable[cubeindex][n + 2]];

            contour->vertices[vertexindex1].triangles[contour->vertices[vertexindex1].nTriangles++] = contour->nTriangles;
            contour->vertices[vertexindex2].triangles[contour->vertices[vertexindex2].nTriangles++] = contour->nTriangles;
            contour->vertices[vertexindex3].triangles[contour->vertices[vertexindex3].nTriangles++] = contour->nTriangles;

            DynamicallyTriangleGrow(&contour->triangles, contour->nTriangles, 1);
            contour->triangles[contour->nTriangles].p[0] = vertexindex1;
            contour->triangles[contour->nTriangles].p[1] = vertexindex2;
            contour->triangles[contour->nTriangles].p[2] = vertexindex3;

            contour->triangles[contour->nTriangles].cas = 0;
            contour->triangles[contour->nTriangles].nEdges = 0;
            contour->triangles[contour->nTriangles].nNeigh = 0;
            contour->triangles[contour->nTriangles].split = 0;
            contour->triangles[contour->nTriangles].index = 0;
            contour->nTriangles++;

            check1 = check2 = 0;
            for (m = 0; m < contour->vertices[vertexindex1].nNeigh; m++)
            {
              if (contour->vertices[vertexindex1].neigh[m] == vertexindex2) check1 = 1;
              if (contour->vertices[vertexindex1].neigh[m] == vertexindex3) check2 = 1;
            }
            if (!check1) contour->vertices[vertexindex1].neigh[contour->vertices[vertexindex1].nNeigh++] = vertexindex2;
            if (!check2) contour->vertices[vertexindex1].neigh[contour->vertices[vertexindex1].nNeigh++] = vertexindex3;

            check1 = check2 = 0;
            for (m = 0; m < contour->vertices[vertexindex2].nNeigh; m++)
            {
              if (contour->vertices[vertexindex2].neigh[m] == vertexindex1) check1 = 1;
              if (contour->vertices[vertexindex2].neigh[m] == vertexindex3) check2 = 1;
            }
            if (!check1) contour->vertices[vertexindex2].neigh[contour->vertices[vertexindex2].nNeigh++] = vertexindex1;
            if (!check2) contour->vertices[vertexindex2].neigh[contour->vertices[vertexindex2].nNeigh++] = vertexindex3;

            check1 = check2 = 0;
            for (m = 0; m < contour->vertices[vertexindex3].nNeigh; m++)
            {
              if (contour->vertices[vertexindex3].neigh[m] == vertexindex1) check1 = 1;
              if (contour->vertices[vertexindex3].neigh[m] == vertexindex2) check2 = 1;
            }
            if (!check1) contour->vertices[vertexindex3].neigh[contour->vertices[vertexindex3].nNeigh++] = vertexindex1;
            if (!check2) contour->vertices[vertexindex3].neigh[contour->vertices[vertexindex3].nNeigh++] = vertexindex2;
          }
        }
      }
    }
  }

  free(edgx);
  free(edgy);
  free(edgz);

  return 0;
}

struct XYZ VertexInterpHalf(double isolevel, struct XYZ p1, struct XYZ p2, double valp1, double valp2)
{
  double mu;
  struct XYZ p;

  mu = 0.5;
  p.x = p1.x + mu * (p2.x - p1.x);
  p.y = p1.y + mu * (p2.y - p1.y);
  p.z = p1.z + mu * (p2.z - p1.z);

  p.mu = mu;

  p.nNeigh = 0;
  p.nTriangles = 0;

  return (p);
}

/* Linearly interpolate the position where an isosurface cuts
 an edge between two contour->vertices, each with their own scalar value */

struct XYZ VertexInterp(double isolevel, struct XYZ p1, struct XYZ p2, double valp1, double valp2)
{
  double mu;
  struct XYZ p;

  if (fabs(isolevel - valp1) < 0.00001) return (p1);
  if (fabs(isolevel - valp2) < 0.00001) return (p2);
  if (fabs(valp1 - valp2) < 0.00001) return (p1);

  mu = (isolevel - valp1) / (valp2 - valp1);
  p.x = p1.x + mu * (p2.x - p1.x);
  p.y = p1.y + mu * (p2.y - p1.y);
  p.z = p1.z + mu * (p2.z - p1.z);

  p.mu = mu;

  p.nNeigh = 0;
  p.nTriangles = 0;
  p.cas = 0;

  return (p);
}

struct XYZ VertexInterpManson(double isolevel, struct XYZ p1, struct XYZ p2, double valp1, double valp2)
{
  struct XYZ p;
  long int inv;
  double mu, a1, a2, a1b, a2b;

  if (valp2 > valp1)
  {
    a1 = valp1;
    a2 = valp2;
    a1b = 1.0 - a2;
    a2b = 1.0 - a1;
    inv = 0;
  }
  else
  {
    a1 = valp2;
    a2 = valp1;
    a1b = 1.0 - a2;
    a2b = 1.0 - a1;
    inv = 1;
  }

  if (a2 <= 3.0 * a1)
  {
    if (3.0 * a2 <= a1 + 2.0)
    {
      /* CASE 1*/
      mu = (a1 - 0.5) / (a1 - a2);
      p.cas = 1;
    }
    else
    {
      if (SQ(2.0 * a1b + 2.0 * a2b - 1.0) < 4.0 * a1b * (a1b + a2b))
      {
        /* CASE 4 */
        mu = (2.0 * a2b - 1.0) / (8.0 * a1b + 4.0 * a2b - 8.0 * sqrt(a1b * (a1b + a2b)));
        p.cas = 4;
      }
      else
      {
        /* CASE 2 */
        mu = 1.5 - a1 - a2;
        p.cas = 2;
      }
    }
  }
  else
  {
    if (3.0 * a2 > a1 + 2.0)
    {
      /* CASE 2 */
      mu = 1.5 - a1 - a2;
      p.cas = 2;
    }
    else
    {
      if (SQ(2.0 * a1 + 2.0 * a2 - 1.0) < 4.0 * a1 * (a1 + a2))
      {
        /* CASE 3 */
        mu = 1.0 - ((2.0 * a2 - 1.0) / (8.0 * a1 + 4.0 * a2 - 8.0 * sqrt(a1 * (a1 + a2))));
        p.cas = 3;
      }
      else
      {
        /* CASE 2 */
        mu = 1.5 - a1 - a2;
        p.cas = 2;
      }
    }
  }

  /* Making sure vertex won't coincide with sampling grid point */
  if (mu > 1.0 - VERYSMALL)
    mu = 1.0 - VERYSMALL;
  else if (mu < VERYSMALL)
    mu = VERYSMALL;

  if (inv)
  {
    mu = 1.0 - mu;
  }

  p.x = p1.x + mu * (p2.x - p1.x);
  p.y = p1.y + mu * (p2.y - p1.y);
  p.z = p1.z + mu * (p2.z - p1.z);

  p.mu = mu;

  p.nNeigh = 0;
  p.nTriangles = 0;

  return (p);
}

int DynamicallyTriangleGrow(struct TRIANGLE **A, long int oldsize, long int increase)
{
  struct TRIANGLE *tmp;
  long int newalloc;

  if (increase == 0) return 0;

  if (oldsize == 0 || (oldsize / TRI_PREALLOCATION) < ((oldsize + increase) / TRI_PREALLOCATION))
  {
    newalloc = TRI_PREALLOCATION * (((increase + oldsize) / TRI_PREALLOCATION) + 1);
    tmp = (struct TRIANGLE *) malloc(newalloc * sizeof(struct TRIANGLE));
    if (oldsize > 0)
    {
      memcpy(tmp, *A, oldsize * sizeof(struct TRIANGLE));
      free(*A);
    }
    *A = tmp;
  }

  return 0;
}

int DynamicallyVerticeGrow(struct XYZ **A, long int oldsize, long int increase)
{
  struct XYZ *tmp;
  long int newalloc;

  if (increase == 0) return 0;

  if (oldsize == 0 || (oldsize / TRI_PREALLOCATION) < ((oldsize + increase) / TRI_PREALLOCATION))
  {
    newalloc = TRI_PREALLOCATION * (((increase + oldsize) / TRI_PREALLOCATION) + 1);
    tmp = (struct XYZ *) malloc(newalloc * sizeof(struct XYZ));
    if (oldsize > 0)
    {
      memcpy(tmp, *A, oldsize * sizeof(struct XYZ));
      free(*A);
    }
    *A = tmp;
  }

  return 0;
}

int DynamicallyEdgeGrow(struct EDGE **A, long int oldsize, long int increase)
{
  struct EDGE *tmp;
  long int newalloc;

  if (increase == 0) return 0;

  if (oldsize == 0 || (oldsize / TRI_PREALLOCATION) < ((oldsize + increase) / TRI_PREALLOCATION))
  {
    newalloc = TRI_PREALLOCATION * (((increase + oldsize) / TRI_PREALLOCATION) + 1);
    tmp = (struct EDGE *) malloc(newalloc * sizeof(struct EDGE));
    if (oldsize > 0)
    {
      memcpy(tmp, *A, oldsize * sizeof(struct EDGE));
      free(*A);
    }
    *A = tmp;
  }

  return 0;
}
